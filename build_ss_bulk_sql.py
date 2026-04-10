#!/usr/bin/env python3
import argparse
import re
import sqlite3
from pathlib import Path

import pandas as pd

import app


def create_global_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS expression_raw (
            gse_id TEXT,
            sample_label TEXT,
            gene TEXT,
            value REAL,
            source_file TEXT
        );

        CREATE INDEX IF NOT EXISTS idx_raw_gse ON expression_raw(gse_id);
        CREATE INDEX IF NOT EXISTS idx_raw_gse_sample ON expression_raw(gse_id, sample_label);
        CREATE INDEX IF NOT EXISTS idx_raw_gene ON expression_raw(gene);

        CREATE TABLE IF NOT EXISTS build_log (
            source_file TEXT PRIMARY KEY,
            gse_id TEXT,
            row_count INTEGER,
            loaded_at TEXT,
            note TEXT
        );
        """
    )
    conn.commit()


def extract_gse_from_name(name: str) -> str:
    m = re.search(r"(GSE\d+)", name.upper())
    return m.group(1) if m else ""


def list_all_expression_like_files(root: Path) -> list[Path]:
    files = []
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        lname = p.name.lower()
        if lname.endswith(".tar") or ".raw" in lname:
            continue
        ok_ext = any(lname.endswith(ext) for ext in app.EXPR_EXTENSIONS)
        if not ok_ext:
            continue
        if not re.search(r"GSE\d+", p.name.upper()):
            continue
        score = 0
        if any(h in lname for h in app.EXPR_HINTS):
            score += 2
        if "all" in lname or "matrix" in lname or "merged" in lname:
            score += 1
        if "deg" in lname or "diff" in lname:
            score -= 2
        if score < 1:
            continue
        files.append((score, str(p)))

    files.sort(key=lambda x: (-x[0], x[1]))
    return [Path(x[1]) for x in files]


def process_one_file(conn: sqlite3.Connection, f: Path) -> tuple[int, str]:
    gse_id = extract_gse_from_name(f.name)
    if not gse_id:
        return 0, "no gse id"

    df = app.read_tabular(f)
    df = df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed")]
    if not app.looks_like_expression_matrix(df):
        return 0, "not expression matrix"

    gene_col = df.columns[0]
    sub = df.copy()
    sub[gene_col] = sub[gene_col].astype(str).map(app.fix_mojibake)
    sub = sub[sub[gene_col].str.len() > 0]

    value_cols = []
    for c in sub.columns[1:]:
        numeric_ratio = pd.to_numeric(sub[c], errors="coerce").notna().mean()
        if numeric_ratio > 0.5:
            value_cols.append(c)

    if not value_cols:
        return 0, "no numeric sample columns"

    rows = []
    for c in value_cols:
        vals = pd.to_numeric(sub[c], errors="coerce")
        genes = sub[gene_col]
        label = app.fix_mojibake(str(c))
        for gene, val in zip(genes, vals):
            if pd.isna(val):
                continue
            rows.append((gse_id, label, app.fix_mojibake(str(gene)), float(val), str(f)))

    if rows:
        conn.executemany(
            """
            INSERT INTO expression_raw (gse_id, sample_label, gene, value, source_file)
            VALUES (?, ?, ?, ?, ?)
            """,
            rows,
        )
        conn.execute(
            """
            INSERT OR REPLACE INTO build_log (source_file, gse_id, row_count, loaded_at, note)
            VALUES (?, ?, ?, datetime('now'), ?)
            """,
            (str(f), gse_id, len(rows), "ok"),
        )
        conn.commit()
    return len(rows), "ok"


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a global SQL from SS_Bulk expression-like tables")
    parser.add_argument("--ss-bulk", default=str(app.resolve_ss_bulk_dir()))
    parser.add_argument("--out", default=str(app.resolve_global_expression_db()))
    parser.add_argument("--max-files", type=int, default=0, help="0 means no limit")
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()

    root = Path(args.ss_bulk).expanduser().resolve()
    out = Path(args.out).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(out)
    create_global_schema(conn)

    if args.rebuild:
        conn.execute("DELETE FROM expression_raw")
        conn.execute("DELETE FROM build_log")
        conn.commit()

    files = list_all_expression_like_files(root)
    if args.max_files > 0:
        files = files[: args.max_files]

    print(f"[INFO] candidate files: {len(files)}")
    total_rows = 0
    ok_files = 0

    for i, f in enumerate(files, start=1):
        try:
            n, note = process_one_file(conn, f)
            if n > 0:
                ok_files += 1
                total_rows += n
            print(f"[{i}/{len(files)}] {f.name} -> {n} rows ({note})")
        except Exception as e:
            print(f"[{i}/{len(files)}] {f.name} -> ERROR: {e}")

    print(f"[DONE] db={out}")
    print(f"[DONE] ok_files={ok_files}, total_rows={total_rows}")

    conn.close()


if __name__ == "__main__":
    main()
