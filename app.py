#!/usr/bin/env python3
import argparse
import csv
import json
import os
import re
import sqlite3
from functools import lru_cache
from datetime import datetime
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pandas as pd

BASE_DIR = Path(__file__).resolve().parent
STATIC_DIR = BASE_DIR / "static"
CSV_PATH = BASE_DIR / "data" / "extraction_batch_result.csv"
DB_PATH = BASE_DIR / "geo_samples.db"
SS_BULK_DIR = BASE_DIR.parent / "SS_Bulk"
DEMO_EXPR_DIR = BASE_DIR / "data" / "expression"
GLOBAL_EXPR_DB = BASE_DIR / "data" / "ss_bulk_expression.db"
MASTER_MAPPING_CSV = BASE_DIR.parent / "SS_Bulk" / "gse_gsm_condition_master_mapping.csv"
LABEL_MAPPING_CSV = BASE_DIR.parent / "SS_Bulk" / "gse_expression_label_mapping.csv"
MAPPING_SUMMARY_CSV = BASE_DIR.parent / "SS_Bulk" / "gse_mapping_alignment_summary.csv"
DEMO_MASTER_MAPPING_CSV = BASE_DIR / "data" / "mapping" / "gse_gsm_condition_master_mapping_demo.csv"
DEMO_LABEL_MAPPING_CSV = BASE_DIR / "data" / "mapping" / "gse_expression_label_mapping_demo.csv"
DEMO_MAPPING_SUMMARY_CSV = BASE_DIR / "data" / "mapping" / "gse_mapping_alignment_summary_demo.csv"

EXPR_EXTENSIONS = {".xlsx", ".xls", ".csv", ".tsv", ".txt", ".xlsx.gz", ".xls.gz", ".csv.gz", ".tsv.gz", ".txt.gz"}
EXPR_HINTS = ("fpkm", "rpkm", "tpm", "count", "counts", "expression", "normalized", "matrix")


# ---------- text helpers ----------
def fix_mojibake(text: str) -> str:
    if text is None:
        return ""
    s = str(text)
    if not s:
        return ""

    # common one: Î -> Δ
    replacements = {
        "Î\x94": "Δ",
        "Î": "Δ",
        "â": "'",
        "â": "-",
        "â": "-",
        "Âµ": "µ",
    }
    for a, b in replacements.items():
        s = s.replace(a, b)

    # heuristic latin1->utf8 roundtrip for typical mojibake text
    if re.search(r"[ÃÂÎÐ]", s):
        try:
            candidate = s.encode("latin1", errors="ignore").decode("utf-8", errors="ignore")
            if candidate and ("�" not in candidate):
                s = candidate
        except Exception:
            pass

    return s.strip()


def norm_text(text: str) -> str:
    s = fix_mojibake(text).lower()
    s = re.sub(r"[^a-z0-9]+", "", s)
    return s


def clean_text(text: str) -> str:
    s = fix_mojibake("" if text is None else str(text)).strip()
    if s.lower() in {"nan", "none", "null"}:
        return ""
    return s


def text_tokens(text: str) -> list:
    s = clean_text(text).lower()
    return [x for x in re.split(r"[^a-z0-9]+", s) if x]


def expand_tokens(tokens: list) -> set:
    t = set(tokens)
    if "wt" in t:
        t.update({"wild", "type"})
    if "q" in t:
        t.update({"q", "quiescence"})
    if "ds" in t:
        t.update({"ds", "diauxic"})
    return t


def resolve_csv_path() -> Path:
    env_path = os.getenv("GEO_CSV_PATH", "").strip()
    if env_path:
        p = Path(env_path).expanduser().resolve()
        if p.exists():
            return p

    if CSV_PATH.exists():
        return CSV_PATH

    legacy = BASE_DIR.parent / "extraction_batch_result.csv"
    if legacy.exists():
        return legacy

    return CSV_PATH


def resolve_ss_bulk_dir() -> Path:
    env_path = os.getenv("SS_BULK_DIR", "").strip()
    if env_path:
        p = Path(env_path).expanduser().resolve()
        if p.exists():
            return p
    if SS_BULK_DIR.exists():
        return SS_BULK_DIR
    return DEMO_EXPR_DIR


def resolve_global_expression_db() -> Path:
    env_path = os.getenv("GLOBAL_EXPR_DB", "").strip()
    if env_path:
        return Path(env_path).expanduser().resolve()
    return GLOBAL_EXPR_DB


def resolve_master_mapping_csv() -> Path:
    env_path = os.getenv("MASTER_MAPPING_CSV", "").strip()
    if env_path:
        return Path(env_path).expanduser().resolve()
    if MASTER_MAPPING_CSV.exists():
        return MASTER_MAPPING_CSV
    return DEMO_MASTER_MAPPING_CSV


def resolve_label_mapping_csv() -> Path:
    env_path = os.getenv("LABEL_MAPPING_CSV", "").strip()
    if env_path:
        return Path(env_path).expanduser().resolve()
    if LABEL_MAPPING_CSV.exists():
        return LABEL_MAPPING_CSV
    return DEMO_LABEL_MAPPING_CSV


def resolve_mapping_summary_csv() -> Path:
    env_path = os.getenv("MAPPING_SUMMARY_CSV", "").strip()
    if env_path:
        return Path(env_path).expanduser().resolve()
    if MAPPING_SUMMARY_CSV.exists():
        return MAPPING_SUMMARY_CSV
    return DEMO_MAPPING_SUMMARY_CSV


@lru_cache(maxsize=1)
def load_master_mapping_rows() -> list:
    path = resolve_master_mapping_csv()
    if not path.exists():
        return []
    try:
        df = pd.read_csv(path, dtype=str).fillna("")
    except Exception:
        return []

    needed = {"gse_id", "gsm_id", "sample_title", "condition_text"}
    if not needed.issubset(set(df.columns)):
        return []

    rows = []
    for _, r in df.iterrows():
        gse_id = clean_text(r.get("gse_id", "")).upper()
        gsm_id = clean_text(r.get("gsm_id", "")).upper()
        if not gse_id or not gsm_id:
            continue
        sample_title = clean_text(r.get("sample_title", ""))
        condition_text = clean_text(r.get("condition_text", ""))
        rows.append(
            {
                "gse_id": gse_id,
                "gsm_id": gsm_id,
                "sample_title": sample_title,
                "sample_title_norm": norm_text(sample_title),
                "condition_text": condition_text,
            }
        )
    return rows


@lru_cache(maxsize=1)
def load_label_mapping_dict() -> dict:
    path = resolve_label_mapping_csv()
    if not path.exists():
        return {}
    try:
        df = pd.read_csv(path, dtype=str).fillna("")
    except Exception:
        return {}

    needed = {"gse_id", "sample_label", "match_status", "matched_gsm_ids"}
    if not needed.issubset(set(df.columns)):
        return {}

    out = {}
    for _, r in df.iterrows():
        status = clean_text(r.get("match_status", ""))
        if not status.startswith("MATCH_"):
            continue
        gse_id = clean_text(r.get("gse_id", "")).upper()
        label = clean_text(r.get("sample_label", ""))
        if not gse_id or not label:
            continue
        key = (gse_id, norm_text(label))
        gsms = [x.strip().upper() for x in clean_text(r.get("matched_gsm_ids", "")).split(";") if x.strip()]
        if not gsms:
            continue
        old = out.get(key, [])
        out[key] = sorted(set(old + gsms))
    return out


# ---------- schema ----------
def create_tables(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS mapping_status (
            gse_id TEXT PRIMARY KEY,
            mapping_status TEXT,
            label_count INTEGER,
            sample_like_count INTEGER,
            mapped_any INTEGER,
            unmatched INTEGER
        );

        CREATE INDEX IF NOT EXISTS idx_mapping_status ON mapping_status(mapping_status);

        CREATE TABLE IF NOT EXISTS mapping_samples (
            gse_id TEXT,
            gsm_id TEXT,
            sample_title TEXT,
            condition_text TEXT,
            mapping_source TEXT,
            PRIMARY KEY (gse_id, gsm_id)
        );

        CREATE INDEX IF NOT EXISTS idx_mapping_gse ON mapping_samples(gse_id);
        CREATE INDEX IF NOT EXISTS idx_mapping_gsm ON mapping_samples(gsm_id);

        CREATE TABLE IF NOT EXISTS samples (
            gse_id TEXT,
            gsm_id TEXT,
            sample_title TEXT,
            status TEXT,
            raw_characteristics TEXT,
            parse_quality TEXT,
            source_name TEXT,
            growth_condition TEXT,
            cell_line TEXT,
            group_name TEXT,
            genotype TEXT,
            treatment TEXT,
            drug TEXT,
            concentration TEXT,
            time_value TEXT,
            other_info TEXT,
            error_msg TEXT,
            extra_fields_json TEXT,
            PRIMARY KEY (gse_id, gsm_id)
        );

        CREATE INDEX IF NOT EXISTS idx_samples_gse ON samples(gse_id);
        CREATE INDEX IF NOT EXISTS idx_samples_gsm ON samples(gsm_id);
        CREATE INDEX IF NOT EXISTS idx_samples_title ON samples(sample_title);

        CREATE TABLE IF NOT EXISTS sample_features (
            gse_id TEXT,
            gsm_id TEXT,
            feature_key TEXT,
            feature_value TEXT,
            source TEXT,
            PRIMARY KEY (gse_id, gsm_id, feature_key, feature_value)
        );

        CREATE INDEX IF NOT EXISTS idx_features_gse_gsm ON sample_features(gse_id, gsm_id);

        CREATE TABLE IF NOT EXISTS expression_values (
            gse_id TEXT,
            gsm_id TEXT,
            sample_label TEXT,
            sample_title TEXT,
            condition_text TEXT,
            gene TEXT,
            value REAL,
            source_file TEXT
        );

        CREATE INDEX IF NOT EXISTS idx_expr_gse ON expression_values(gse_id);
        CREATE INDEX IF NOT EXISTS idx_expr_gse_gsm ON expression_values(gse_id, gsm_id);
        CREATE INDEX IF NOT EXISTS idx_expr_gene ON expression_values(gene);

        CREATE TABLE IF NOT EXISTS expression_load_log (
            gse_id TEXT PRIMARY KEY,
            loaded_at TEXT,
            file_count INTEGER,
            row_count INTEGER,
            notes TEXT
        );
        """
    )
    conn.commit()


def parse_raw_characteristics(raw_text: str) -> dict:
    data = {}
    if not raw_text:
        return data
    for item in str(raw_text).split("|"):
        s = item.strip()
        if not s:
            continue
        if ":" in s:
            k, v = s.split(":", 1)
            key = fix_mojibake(k).strip().lower().replace(" ", "_")
            val = fix_mojibake(v).strip()
            if key:
                data[key] = val
    return data


def to_feature_pairs(row: dict) -> list:
    pairs = []
    standardized = [
        ("source_name", row.get("source_name", "")),
        ("growth_condition", row.get("growth_condition", "")),
        ("cell_line", row.get("cell_line", "")),
        ("group", row.get("group", "")),
        ("genotype", row.get("genotype", "")),
        ("treatment", row.get("treatment", "")),
        ("drug", row.get("drug", "")),
        ("concentration", row.get("concentration", "")),
        ("time", row.get("time", "")),
        ("other_info", row.get("other_info", "")),
    ]
    for key, val in standardized:
        text = fix_mojibake(val)
        if text:
            pairs.append((key, text, "standardized"))

    from_raw = parse_raw_characteristics(row.get("Raw_Characteristics", ""))
    for key, val in from_raw.items():
        pairs.append((key, val, "raw_characteristics"))

    extra_json_text = row.get("extra_fields_json", "") or ""
    if extra_json_text:
        try:
            obj = json.loads(extra_json_text)
            for key, val in obj.items():
                text = fix_mojibake("" if val is None else str(val))
                if text:
                    pairs.append((str(key), text, "extra_json"))
        except Exception:
            pass

    dedup = []
    seen = set()
    for p in pairs:
        if p not in seen:
            seen.add(p)
            dedup.append(p)
    return dedup


# ---------- import sample metadata ----------
def import_csv_to_db(csv_path: Path, db_path: Path) -> None:
    conn = sqlite3.connect(db_path)
    conn.executescript(
        """
        DROP TABLE IF EXISTS mapping_status;
        DROP TABLE IF EXISTS mapping_samples;
        DROP TABLE IF EXISTS sample_features;
        DROP TABLE IF EXISTS samples;
        DROP TABLE IF EXISTS expression_values;
        DROP TABLE IF EXISTS expression_load_log;
        """
    )
    create_tables(conn)

    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        sample_rows = []
        feature_rows = []

        for row in reader:
            gse_id = (row.get("GSE_ID", "") or "").strip().upper()
            gsm_id = (row.get("GSM_ID", "") or "").strip().upper()
            sample_title = fix_mojibake(row.get("Sample_Title", ""))
            raw = fix_mojibake(row.get("Raw_Characteristics", ""))

            sample_rows.append(
                (
                    gse_id,
                    gsm_id,
                    sample_title,
                    row.get("Status", ""),
                    raw,
                    row.get("Parse_Quality", ""),
                    fix_mojibake(row.get("source_name", "")),
                    fix_mojibake(row.get("growth_condition", "")),
                    fix_mojibake(row.get("cell_line", "")),
                    fix_mojibake(row.get("group", "")),
                    fix_mojibake(row.get("genotype", "")),
                    fix_mojibake(row.get("treatment", "")),
                    fix_mojibake(row.get("drug", "")),
                    fix_mojibake(row.get("concentration", "")),
                    fix_mojibake(row.get("time", "")),
                    fix_mojibake(row.get("other_info", "")),
                    fix_mojibake(row.get("Error_Msg", "")),
                    row.get("extra_fields_json", ""),
                )
            )

            if gse_id and gsm_id:
                for key, val, source in to_feature_pairs(row):
                    feature_rows.append((gse_id, gsm_id, key, val, source))

    conn.executemany(
        """
        INSERT INTO samples (
            gse_id, gsm_id, sample_title, status, raw_characteristics, parse_quality,
            source_name, growth_condition, cell_line, group_name, genotype, treatment,
            drug, concentration, time_value, other_info, error_msg, extra_fields_json
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        sample_rows,
    )

    conn.executemany(
        """
        INSERT OR IGNORE INTO sample_features (gse_id, gsm_id, feature_key, feature_value, source)
        VALUES (?, ?, ?, ?, ?)
        """,
        feature_rows,
    )

    # Build a DB-native master mapping table for all GSE/GSM pairs.
    merged = {}
    for r in load_master_mapping_rows():
        key = (r["gse_id"], r["gsm_id"])
        merged[key] = {
            "gse_id": r["gse_id"],
            "gsm_id": r["gsm_id"],
            "sample_title": clean_text(r.get("sample_title", "")),
            "condition_text": clean_text(r.get("condition_text", "")),
            "mapping_source": "MASTER_MAPPING",
        }

    local_rows = conn.execute(
        """
        SELECT gse_id, gsm_id, sample_title,
               COALESCE(treatment,''), COALESCE(genotype,''), COALESCE(raw_characteristics,'')
        FROM samples
        WHERE gse_id <> '' AND gsm_id <> ''
        """
    ).fetchall()
    for gse_id, gsm_id, sample_title, treatment, genotype, raw in local_rows:
        key = (gse_id, gsm_id)
        cond_local = " | ".join([x for x in [treatment, genotype, raw] if clean_text(x)]).strip(" |")
        cur = merged.get(
            key,
            {
                "gse_id": gse_id,
                "gsm_id": gsm_id,
                "sample_title": "",
                "condition_text": "",
                "mapping_source": "LOCAL_META",
            },
        )
        if clean_text(sample_title):
            cur["sample_title"] = clean_text(sample_title)
        if clean_text(cond_local):
            cur["condition_text"] = clean_text(cond_local)
        if cur.get("mapping_source") != "MASTER_MAPPING":
            cur["mapping_source"] = "LOCAL_META"
        merged[key] = cur

    mapping_rows = [
        (
            x["gse_id"],
            x["gsm_id"],
            x["sample_title"],
            x["condition_text"],
            x["mapping_source"],
        )
        for x in merged.values()
        if x["gse_id"] and x["gsm_id"]
    ]
    conn.executemany(
        """
        INSERT OR REPLACE INTO mapping_samples (gse_id, gsm_id, sample_title, condition_text, mapping_source)
        VALUES (?, ?, ?, ?, ?)
        """,
        mapping_rows,
    )

    # optional summary table for UI status badge
    summary_path = resolve_mapping_summary_csv()
    if summary_path.exists():
        try:
            sm = pd.read_csv(summary_path, dtype=str).fillna("")
            if {"gse_id", "mapping_status"}.issubset(set(sm.columns)):
                status_rows = []
                for _, r in sm.iterrows():
                    status_rows.append(
                        (
                            clean_text(r.get("gse_id", "")).upper(),
                            clean_text(r.get("mapping_status", "")),
                            int(clean_text(r.get("label_count", "0")) or 0),
                            int(clean_text(r.get("sample_like_count", "0")) or 0),
                            int(clean_text(r.get("mapped_any", "0")) or 0),
                            int(clean_text(r.get("unmatched", "0")) or 0),
                        )
                    )
                conn.executemany(
                    """
                    INSERT OR REPLACE INTO mapping_status (
                        gse_id, mapping_status, label_count, sample_like_count, mapped_any, unmatched
                    ) VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    status_rows,
                )
        except Exception:
            pass

    conn.commit()
    # Reclaim disk space after repeated rebuilds; otherwise SQLite file can bloat.
    conn.execute("VACUUM")
    conn.close()


# ---------- sample search ----------
def search_samples(conn: sqlite3.Connection, params: dict) -> list:
    gse = params.get("gse", "").strip().upper()
    gsm = params.get("gsm", "").strip().upper()
    q = fix_mojibake(params.get("q", "").strip()).lower()
    limit = int(params.get("limit", "100") or 100)
    limit = max(1, min(limit, 500))

    where, args = [], []
    if gse:
        where.append("m.gse_id = ?")
        args.append(gse)
    if gsm:
        where.append("m.gsm_id = ?")
        args.append(gsm)
    if q:
        where.append(
            "(" \
            "LOWER(COALESCE(NULLIF(s.sample_title,''), m.sample_title)) LIKE ? OR " \
            "LOWER(COALESCE(s.raw_characteristics, m.condition_text, '')) LIKE ? OR " \
            "LOWER(COALESCE(s.treatment, '')) LIKE ? OR " \
            "LOWER(COALESCE(s.genotype, '')) LIKE ? OR " \
            "LOWER(m.gse_id) LIKE ? OR " \
            "LOWER(m.gsm_id) LIKE ?" \
            ")"
        )
        like = f"%{q}%"
        args.extend([like, like, like, like, like, like])

    where_sql = " AND ".join(where) if where else "1=1"
    sql = f"""
        SELECT m.gse_id, m.gsm_id,
               COALESCE(NULLIF(s.sample_title,''), m.sample_title) AS sample_title,
               COALESCE(s.treatment, '') AS treatment,
               COALESCE(s.genotype, '') AS genotype,
               COALESCE(NULLIF(s.raw_characteristics,''), m.condition_text, '') AS raw_characteristics,
               COUNT(f.feature_key) AS feature_count,
               COALESCE(ms.mapping_status, 'UNKNOWN') AS mapping_status
        FROM mapping_samples m
        LEFT JOIN samples s ON s.gse_id = m.gse_id AND s.gsm_id = m.gsm_id
        LEFT JOIN sample_features f ON f.gse_id = m.gse_id AND f.gsm_id = m.gsm_id
        LEFT JOIN mapping_status ms ON ms.gse_id = m.gse_id
        WHERE {where_sql}
        GROUP BY m.gse_id, m.gsm_id
        ORDER BY m.gse_id, m.gsm_id
        LIMIT ?
    """
    rows = conn.execute(sql, (*args, limit)).fetchall()
    return [
        {
            "gse_id": r[0],
            "gsm_id": r[1],
            "sample_title": r[2],
            "treatment": r[3],
            "genotype": r[4],
            "raw_characteristics": r[5],
            "feature_count": r[6],
            "mapping_status": r[7],
        }
        for r in rows
    ]


def get_sample_detail(conn: sqlite3.Connection, gsm_id: str) -> dict:
    gsm_id = (gsm_id or "").strip().upper()
    row = conn.execute(
        """
        SELECT gse_id, gsm_id, sample_title, raw_characteristics,
               source_name, growth_condition, cell_line, group_name, genotype, treatment,
               drug, concentration, time_value, other_info
        FROM samples
        WHERE gsm_id = ?
        ORDER BY CASE WHEN gse_id = '' THEN 1 ELSE 0 END, LENGTH(raw_characteristics) DESC
        LIMIT 1
        """,
        (gsm_id,),
    ).fetchone()
    if not row:
        return {}

    features = conn.execute(
        """
        SELECT feature_key, feature_value, source
        FROM sample_features
        WHERE gse_id = ? AND gsm_id = ?
        ORDER BY feature_key, feature_value
        LIMIT 1000
        """,
        (row[0], row[1]),
    ).fetchall()

    return {
        "sample": {
            "gse_id": row[0], "gsm_id": row[1], "sample_title": row[2],
            "raw_characteristics": row[3], "source_name": row[4], "growth_condition": row[5],
            "cell_line": row[6], "group": row[7], "genotype": row[8], "treatment": row[9],
            "drug": row[10], "concentration": row[11], "time": row[12], "other_info": row[13],
        },
        "gene_like_rows": [{"gene": f[0], "value": f[1], "source": f[2]} for f in features],
    }


# ---------- gse overview ----------
def get_gse_basic_table(conn: sqlite3.Connection, gse_id: str) -> list:
    gse_id = (gse_id or "").strip().upper()
    if not gse_id:
        return []

    samples = get_gse_samples(conn, gse_id)
    return [
        {
            "gse_id": x["gse_id"],
            "gsm_id": x["gsm_id"],
            "sample_name": x["sample_title"],
            "condition": x["condition_text"],
        }
        for x in samples
    ]


def get_gse_mapping_status(conn: sqlite3.Connection, gse_id: str) -> dict:
    gse_id = (gse_id or "").strip().upper()
    if not gse_id:
        return {"gse_id": "", "mapping_status": "UNKNOWN"}
    row = conn.execute(
        """
        SELECT gse_id, mapping_status, label_count, sample_like_count, mapped_any, unmatched
        FROM mapping_status
        WHERE gse_id = ?
        """,
        (gse_id,),
    ).fetchone()
    if not row:
        return {"gse_id": gse_id, "mapping_status": "UNKNOWN"}
    return {
        "gse_id": row[0],
        "mapping_status": row[1],
        "label_count": row[2],
        "sample_like_count": row[3],
        "mapped_any": row[4],
        "unmatched": row[5],
    }


# ---------- expression matrix ingest ----------
def get_gse_samples(conn: sqlite3.Connection, gse_id: str) -> list:
    gse_id = (gse_id or "").strip().upper()
    if not gse_id:
        return []

    rows = conn.execute(
        """
        SELECT m.gse_id, m.gsm_id,
               COALESCE(NULLIF(s.sample_title,''), m.sample_title) AS sample_title,
               COALESCE(s.treatment,''), COALESCE(s.genotype,''), COALESCE(s.raw_characteristics,''),
               COALESCE(m.condition_text,'')
        FROM mapping_samples m
        LEFT JOIN samples s ON s.gse_id = m.gse_id AND s.gsm_id = m.gsm_id
        WHERE m.gse_id = ?
        ORDER BY m.gsm_id
        """,
        (gse_id,),
    ).fetchall()
    if rows:
        out = []
        for r in rows:
            cond_local = " | ".join([x for x in [r[3], r[4], r[5]] if clean_text(x)]).strip(" |")
            cond = cond_local if cond_local else clean_text(r[6])
            sample_title = clean_text(r[2])
            out.append(
                {
                    "gse_id": r[0],
                    "gsm_id": r[1],
                    "sample_title": sample_title,
                    "sample_title_norm": norm_text(sample_title),
                    "condition_text": cond,
                }
            )
        return out

    # fallback: if mapping_samples table is unavailable, use old local+csv merge path
    rows = conn.execute(
        """
        SELECT gse_id, gsm_id, sample_title,
               COALESCE(treatment,''), COALESCE(genotype,''), COALESCE(raw_characteristics,'')
        FROM samples
        WHERE gse_id = ?
        ORDER BY gsm_id
        """,
        (gse_id,),
    ).fetchall()
    local_samples = []
    for r in rows:
        cond = " | ".join([x for x in [r[3], r[4], r[5]] if x]).strip(" |")
        local_samples.append(
            {
                "gse_id": r[0],
                "gsm_id": r[1],
                "sample_title": clean_text(r[2]),
                "sample_title_norm": norm_text(clean_text(r[2])),
                "condition_text": cond,
            }
        )

    master_rows = [x for x in load_master_mapping_rows() if x["gse_id"] == gse_id]

    merged = {}
    for s in master_rows:
        merged[s["gsm_id"]] = dict(s)
    for s in local_samples:
        cur = merged.get(s["gsm_id"])
        if not cur:
            merged[s["gsm_id"]] = dict(s)
            continue
        # local metadata has higher priority; fill only missing master fields
        if not clean_text(cur.get("sample_title", "")) and clean_text(s.get("sample_title", "")):
            cur["sample_title"] = clean_text(s["sample_title"])
            cur["sample_title_norm"] = norm_text(cur["sample_title"])
        if not clean_text(cur.get("condition_text", "")) and clean_text(s.get("condition_text", "")):
            cur["condition_text"] = clean_text(s["condition_text"])
        merged[s["gsm_id"]] = cur

    out = list(merged.values())
    out.sort(key=lambda x: x.get("gsm_id", ""))
    return out


def list_candidate_expression_files_from_root(gse_id: str, root: Path) -> list:
    if not root.exists():
        return []
    gse_upper = gse_id.upper()
    candidates = []
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        name = p.name
        lname = name.lower()
        if gse_upper not in name.upper():
            continue
        if lname.endswith(".tar") or ".raw" in lname:
            continue
        ok_ext = any(lname.endswith(ext) for ext in EXPR_EXTENSIONS)
        if not ok_ext:
            continue
        score = 0
        if any(h in lname for h in EXPR_HINTS):
            score += 2
        if "all" in lname or "merged" in lname or "matrix" in lname:
            score += 1
        if "deg" in lname or "diff" in lname:
            score -= 2
        candidates.append((score, str(p)))

    candidates.sort(key=lambda x: (-x[0], x[1]))
    return [Path(x[1]) for x in candidates[:20]]


def list_candidate_expression_files(gse_id: str) -> list:
    primary = resolve_ss_bulk_dir()
    files = list_candidate_expression_files_from_root(gse_id, primary)
    if files:
        return files

    # cloud fallback: bundled demo expression files in repo
    if primary != DEMO_EXPR_DIR:
        files = list_candidate_expression_files_from_root(gse_id, DEMO_EXPR_DIR)
        if files:
            return files
    return []


def read_tabular(path: Path) -> pd.DataFrame:
    lower = path.name.lower()
    if lower.endswith(".xlsx") or lower.endswith(".xls"):
        return pd.read_excel(path)

    if lower.endswith(".gz"):
        try:
            df = pd.read_csv(path, sep="\t", compression="gzip")
            if df.shape[1] <= 1:
                df = pd.read_csv(path, sep=",", compression="gzip")
            return df
        except Exception:
            return pd.read_csv(path, sep=",", compression="gzip")

    try:
        df = pd.read_csv(path, sep="\t")
        if df.shape[1] <= 1:
            df = pd.read_csv(path, sep=",")
        return df
    except Exception:
        return pd.read_csv(path, sep=",")


def looks_like_expression_matrix(df: pd.DataFrame) -> bool:
    if df is None or df.empty or df.shape[1] < 2:
        return False
    first_col_name = str(df.columns[0]).lower()
    if any(k in first_col_name for k in ["gene", "symbol", "orf", "locus", "id"]):
        return True

    sample = df.iloc[:50, 0].astype(str)
    alpha_ratio = (sample.str.len() > 0).mean()
    numeric_ratio = pd.to_numeric(df.iloc[:50, 1], errors="coerce").notna().mean()
    return (alpha_ratio > 0.8) and (numeric_ratio > 0.6)


def map_column_to_samples(gse_id: str, col_name: str, gse_samples: list) -> list:
    col_name = str(col_name).strip()
    m = re.search(r"(GSM\d+)", col_name, flags=re.I)
    sample_by_gsm = {s["gsm_id"]: s for s in gse_samples}
    if m:
        gsm = m.group(1).upper()
        s = sample_by_gsm.get(gsm)
        if s:
            return [(gsm, s["sample_title"], s["condition_text"])]
        return [(gsm, "", "")]

    c_norm = norm_text(col_name)
    if not c_norm:
        return []

    # first try precomputed label mapping (handles labels like WT_Log -> multiple GSM replicates)
    label_map = load_label_mapping_dict()
    by_label = label_map.get((gse_id.upper(), c_norm), [])
    if by_label:
        out = []
        for gsm in by_label:
            s = sample_by_gsm.get(gsm)
            if s:
                out.append((gsm, s["sample_title"], s["condition_text"]))
            else:
                out.append((gsm, "", ""))
        return out

    # exact/contains match on sample title
    for s in gse_samples:
        t_norm = s["sample_title_norm"]
        if not t_norm:
            continue
        if c_norm == t_norm or c_norm in t_norm or t_norm in c_norm:
            return [(s["gsm_id"], s["sample_title"], s["condition_text"])]

    # token-based fallback for non-GSM labels
    label_tokens = expand_tokens(text_tokens(col_name))
    scored = []
    ll = col_name.lower().replace("-", "_")
    for s in gse_samples:
        st = clean_text(s.get("sample_title", ""))
        ct = clean_text(s.get("condition_text", ""))
        st_tokens = expand_tokens(text_tokens(st))
        ct_tokens = expand_tokens(text_tokens(ct))
        score = len(label_tokens & st_tokens) + len(label_tokens & ct_tokens)

        stl = st.lower()
        if "wt" in ll and ("wild type" in stl or "wild-type" in stl):
            score += 2
        if "rpd3" in ll and "rpd3" in stl:
            score += 2
        if "_log" in ll and "log" in stl:
            score += 2
        if "_ds" in ll and ("ds" in stl or "diauxic" in stl):
            score += 2
        if ll.endswith("_q") and (" q " in f" {stl} " or "quiescence" in stl):
            score += 2

        if score > 0:
            scored.append((score, s))

    if scored:
        best = max(x[0] for x in scored)
        top = [s for score, s in scored if score == best]
        top.sort(key=lambda x: x.get("gsm_id", ""))
        return [(x["gsm_id"], x["sample_title"], x["condition_text"]) for x in top]

    return []


def load_expression_from_global_db(conn: sqlite3.Connection, gse_id: str, gse_samples: list) -> dict:
    db_path = resolve_global_expression_db()
    if not db_path.exists():
        return {"ok": False, "message": f"global expression db not found: {db_path}"}

    # global DB can be large; allow wait time for transient write locks
    gconn = sqlite3.connect(db_path, timeout=30)
    try:
        exists = gconn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='expression_raw'"
        ).fetchone()
        if not exists:
            return {"ok": False, "message": "global expression db has no table expression_raw"}

        rows = gconn.execute(
            """
            SELECT gse_id, sample_label, gene, value, source_file
            FROM expression_raw
            WHERE gse_id = ?
            """,
            (gse_id,),
        ).fetchall()
        if not rows:
            return {"ok": False, "message": f"{gse_id} not found in global expression db"}

        inserts = []
        for gse, sample_label, gene, value, source_file in rows:
            matches = map_column_to_samples(gse_id, sample_label, gse_samples)
            if not matches:
                matches = [("", "", "")]
            for gsm_id, sample_title, condition_text in matches:
                inserts.append(
                    (
                        gse,
                        gsm_id,
                        sample_label,
                        sample_title,
                        condition_text,
                        gene,
                        float(value),
                        source_file,
                    )
                )

        conn.executemany(
            """
            INSERT INTO expression_values (
                gse_id, gsm_id, sample_label, sample_title, condition_text,
                gene, value, source_file
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            inserts,
        )
        conn.commit()
        return {"ok": True, "message": "loaded from global expression db", "row_count": len(inserts)}
    finally:
        gconn.close()


def ensure_expression_loaded(conn: sqlite3.Connection, gse_id: str) -> dict:
    gse_id = (gse_id or "").strip().upper()
    if not re.fullmatch(r"GSE\d+", gse_id):
        return {"ok": False, "message": "invalid gse"}

    already = conn.execute("SELECT gse_id, loaded_at, file_count, row_count FROM expression_load_log WHERE gse_id = ?", (gse_id,)).fetchone()
    if already:
        return {
            "ok": True,
            "message": "already loaded",
            "gse_id": already[0],
            "loaded_at": already[1],
            "file_count": already[2],
            "row_count": already[3],
        }

    gse_samples = get_gse_samples(conn, gse_id)
    if not gse_samples:
        return {"ok": False, "message": f"{gse_id} not found in sample metadata"}

    # 1) prefer prebuilt global SQL from all SS_Bulk expression tables
    global_loaded = load_expression_from_global_db(conn, gse_id, gse_samples)
    if global_loaded.get("ok"):
        row_count = int(global_loaded.get("row_count", 0))
        conn.execute(
            "INSERT OR REPLACE INTO expression_load_log (gse_id, loaded_at, file_count, row_count, notes) VALUES (?, ?, ?, ?, ?)",
            (gse_id, datetime.now().strftime("%Y-%m-%d %H:%M:%S"), 0, row_count, "global_db"),
        )
        conn.commit()
        return {
            "ok": True,
            "message": "global_db",
            "gse_id": gse_id,
            "file_count": 0,
            "row_count": row_count,
        }

    # 2) fallback to raw file scan
    files = list_candidate_expression_files(gse_id)
    if not files:
        return {
            "ok": False,
            "message": f"no expression-like files found in {resolve_ss_bulk_dir()} or {DEMO_EXPR_DIR}",
        }

    total_rows = 0
    used_files = 0

    for f in files:
        try:
            df = read_tabular(f)
            # clean columns
            df = df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed")]
            if not looks_like_expression_matrix(df):
                continue

            gene_col = df.columns[0]
            sub = df.copy()
            sub[gene_col] = sub[gene_col].astype(str).map(fix_mojibake)
            sub = sub[sub[gene_col].str.len() > 0]

            value_cols = []
            for c in sub.columns[1:]:
                numeric_ratio = pd.to_numeric(sub[c], errors="coerce").notna().mean()
                if numeric_ratio > 0.5:
                    value_cols.append(c)

            if not value_cols:
                continue

            rows_to_insert = []
            for c in value_cols:
                matches = map_column_to_samples(gse_id, str(c), gse_samples)
                if not matches:
                    matches = [("", "", "")]
                vals = pd.to_numeric(sub[c], errors="coerce")
                genes = sub[gene_col]
                for gene, val in zip(genes, vals):
                    if pd.isna(val):
                        continue
                    for gsm_id, sample_title, condition_text in matches:
                        rows_to_insert.append(
                            (
                                gse_id,
                                gsm_id,
                                fix_mojibake(str(c)),
                                sample_title,
                                condition_text,
                                fix_mojibake(str(gene)),
                                float(val),
                                str(f),
                            )
                        )

            if rows_to_insert:
                conn.executemany(
                    """
                    INSERT INTO expression_values (
                        gse_id, gsm_id, sample_label, sample_title, condition_text,
                        gene, value, source_file
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    rows_to_insert,
                )
                total_rows += len(rows_to_insert)
                used_files += 1
                conn.commit()

        except Exception:
            continue

    note = "ok" if total_rows > 0 else "no valid matrix parsed"
    conn.execute(
        "INSERT OR REPLACE INTO expression_load_log (gse_id, loaded_at, file_count, row_count, notes) VALUES (?, ?, ?, ?, ?)",
        (gse_id, datetime.now().strftime("%Y-%m-%d %H:%M:%S"), used_files, total_rows, note),
    )
    conn.commit()

    return {
        "ok": total_rows > 0,
        "message": note,
        "gse_id": gse_id,
        "file_count": used_files,
        "row_count": total_rows,
    }


def search_expression(conn: sqlite3.Connection, params: dict) -> list:
    gse = (params.get("gse", "") or "").strip().upper()
    if not gse:
        return []

    gsm = (params.get("gsm", "") or "").strip().upper()
    sample = fix_mojibake(params.get("sample", "").strip()).lower()
    condition = fix_mojibake(params.get("condition", "").strip()).lower()
    gene = fix_mojibake(params.get("gene", "").strip()).lower()
    limit = int(params.get("limit", "200") or 200)
    limit = max(1, min(limit, 5000))

    where = ["e.gse_id = ?"]
    args = [gse]

    if gsm:
        where.append("e.gsm_id = ?")
        args.append(gsm)
    if sample:
        where.append("(LOWER(e.sample_label) LIKE ? OR LOWER(e.sample_title) LIKE ?)")
        like = f"%{sample}%"
        args.extend([like, like])
    if condition:
        where.append("LOWER(e.condition_text) LIKE ?")
        args.append(f"%{condition}%")
    if gene:
        where.append("LOWER(e.gene) LIKE ?")
        args.append(f"%{gene}%")

    sql = f"""
        SELECT e.gse_id, e.gsm_id, e.sample_label, e.sample_title,
               e.condition_text, e.gene, e.value, e.source_file
        FROM expression_values e
        WHERE {' AND '.join(where)}
        ORDER BY e.sample_label, e.gene
        LIMIT ?
    """
    rows = conn.execute(sql, (*args, limit)).fetchall()
    return [
        {
            "gse_id": r[0],
            "gsm_id": r[1],
            "sample_label": r[2],
            "sample_title": r[3],
            "condition_text": r[4],
            "gene": r[5],
            "value": r[6],
            "source_file": r[7],
        }
        for r in rows
    ]


# ---------- stats ----------
def get_stats(conn: sqlite3.Connection) -> dict:
    n_samples = conn.execute("SELECT COUNT(*) FROM mapping_samples").fetchone()[0]
    n_gse = conn.execute("SELECT COUNT(DISTINCT gse_id) FROM mapping_samples WHERE gse_id <> ''").fetchone()[0]
    n_features = conn.execute("SELECT COUNT(*) FROM sample_features").fetchone()[0]
    n_expr = conn.execute("SELECT COUNT(*) FROM expression_values").fetchone()[0]
    return {
        "samples": n_samples,
        "gses": n_gse,
        "gene_like_rows": n_features,
        "expression_rows": n_expr,
        "db_updated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }


# ---------- HTTP API ----------
class AppHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        self.conn = kwargs.pop("conn")
        super().__init__(*args, directory=str(STATIC_DIR), **kwargs)

    def _send_json(self, payload: dict, status: int = 200) -> None:
        data = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        parsed = urlparse(self.path)

        if parsed.path == "/api/stats":
            self._send_json(get_stats(self.conn))
            return

        if parsed.path == "/api/search":
            params = {k: v[0] for k, v in parse_qs(parsed.query).items()}
            rows = search_samples(self.conn, params)
            self._send_json({"rows": rows, "count": len(rows)})
            return

        if parsed.path == "/api/expression_prepare":
            params = {k: v[0] for k, v in parse_qs(parsed.query).items()}
            gse = params.get("gse", "")
            result = ensure_expression_loaded(self.conn, gse)
            self._send_json(result)
            return

        if parsed.path == "/api/expression_search":
            params = {k: v[0] for k, v in parse_qs(parsed.query).items()}
            gse = (params.get("gse", "") or "").strip().upper()
            if gse:
                ensure_expression_loaded(self.conn, gse)
            rows = search_expression(self.conn, params)
            self._send_json({"rows": rows, "count": len(rows)})
            return

        if parsed.path.startswith("/api/sample/"):
            gsm = parsed.path.rsplit("/", 1)[-1].upper().strip()
            payload = get_sample_detail(self.conn, gsm)
            self._send_json(payload if payload else {"error": "not found"}, 200 if payload else 404)
            return

        if parsed.path in {"/", "/index.html"}:
            self.path = "/index.html"
        return super().do_GET()


# ---------- runtime ----------
def run_server(port: int) -> None:
    csv_path = resolve_csv_path()
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    import_csv_to_db(csv_path, DB_PATH)
    conn = sqlite3.connect(DB_PATH, check_same_thread=False)

    def handler(*args, **kwargs):
        AppHandler(*args, conn=conn, **kwargs)

    server = ThreadingHTTPServer(("0.0.0.0", port), handler)
    print(f"[OK] Server running at http://127.0.0.1:{port}")
    print(f"[OK] DB file: {DB_PATH}")
    print("[INFO] Press Ctrl+C to stop")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        conn.close()
        server.server_close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GSE/GSM prototype query web app")
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    run_server(args.port)
