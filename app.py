#!/usr/bin/env python3
import argparse
import csv
import json
import os
import sqlite3
from datetime import datetime
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

BASE_DIR = Path(__file__).resolve().parent
STATIC_DIR = BASE_DIR / "static"
CSV_PATH = BASE_DIR / "data" / "extraction_batch_result.csv"
DB_PATH = BASE_DIR / "geo_samples.db"


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


def create_tables(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
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
            gsm_id TEXT,
            feature_key TEXT,
            feature_value TEXT,
            source TEXT,
            PRIMARY KEY (gsm_id, feature_key, feature_value)
        );

        CREATE INDEX IF NOT EXISTS idx_features_gsm ON sample_features(gsm_id);
        CREATE INDEX IF NOT EXISTS idx_features_key ON sample_features(feature_key);
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
            key = k.strip().lower().replace(" ", "_")
            val = v.strip()
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
        if val:
            pairs.append((key, str(val), "standardized"))

    from_raw = parse_raw_characteristics(row.get("Raw_Characteristics", ""))
    for key, val in from_raw.items():
        pairs.append((key, val, "raw_characteristics"))

    extra_json_text = row.get("extra_fields_json", "") or ""
    if extra_json_text:
        try:
            obj = json.loads(extra_json_text)
            for key, val in obj.items():
                text = "" if val is None else str(val)
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


def import_csv_to_db(csv_path: Path, db_path: Path) -> None:
    conn = sqlite3.connect(db_path)
    conn.executescript(
        """
        DROP TABLE IF EXISTS sample_features;
        DROP TABLE IF EXISTS samples;
        """
    )
    create_tables(conn)

    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        sample_rows = []
        feature_rows = []

        for row in reader:
            sample_rows.append(
                (
                    row.get("GSE_ID", ""),
                    row.get("GSM_ID", ""),
                    row.get("Sample_Title", ""),
                    row.get("Status", ""),
                    row.get("Raw_Characteristics", ""),
                    row.get("Parse_Quality", ""),
                    row.get("source_name", ""),
                    row.get("growth_condition", ""),
                    row.get("cell_line", ""),
                    row.get("group", ""),
                    row.get("genotype", ""),
                    row.get("treatment", ""),
                    row.get("drug", ""),
                    row.get("concentration", ""),
                    row.get("time", ""),
                    row.get("other_info", ""),
                    row.get("Error_Msg", ""),
                    row.get("extra_fields_json", ""),
                )
            )

            gsm_id = row.get("GSM_ID", "")
            if gsm_id:
                for key, val, source in to_feature_pairs(row):
                    feature_rows.append((gsm_id, key, val, source))

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
        INSERT OR IGNORE INTO sample_features (gsm_id, feature_key, feature_value, source)
        VALUES (?, ?, ?, ?)
        """,
        feature_rows,
    )

    conn.commit()
    conn.close()


def search_samples(conn: sqlite3.Connection, params: dict) -> list:
    gse = params.get("gse", "").strip().upper()
    gsm = params.get("gsm", "").strip().upper()
    q = params.get("q", "").strip().lower()
    limit = int(params.get("limit", "100") or 100)
    limit = max(1, min(limit, 500))

    where = []
    args = []

    if gse:
        where.append("s.gse_id = ?")
        args.append(gse)
    if gsm:
        where.append("s.gsm_id = ?")
        args.append(gsm)
    if q:
        where.append(
            "(" \
            "LOWER(s.sample_title) LIKE ? OR " \
            "LOWER(s.raw_characteristics) LIKE ? OR " \
            "LOWER(s.treatment) LIKE ? OR " \
            "LOWER(s.genotype) LIKE ? OR " \
            "LOWER(s.gse_id) LIKE ? OR " \
            "LOWER(s.gsm_id) LIKE ?" \
            ")"
        )
        like = f"%{q}%"
        args.extend([like, like, like, like, like, like])

    where_sql = " AND ".join(where) if where else "1=1"
    sql = f"""
        SELECT
            s.gse_id,
            s.gsm_id,
            s.sample_title,
            s.treatment,
            s.genotype,
            s.parse_quality,
            s.status,
            s.raw_characteristics,
            COUNT(f.feature_key) AS feature_count
        FROM samples s
        LEFT JOIN sample_features f ON f.gsm_id = s.gsm_id
        WHERE {where_sql}
        GROUP BY s.gse_id, s.gsm_id
        ORDER BY s.gse_id, s.gsm_id
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
            "parse_quality": r[5],
            "status": r[6],
            "raw_characteristics": r[7],
            "feature_count": r[8],
        }
        for r in rows
    ]


def get_sample_detail(conn: sqlite3.Connection, gsm_id: str) -> dict:
    row = conn.execute(
        """
        SELECT gse_id, gsm_id, sample_title, status, parse_quality, raw_characteristics,
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
        WHERE gsm_id = ?
        ORDER BY feature_key, feature_value
        LIMIT 1000
        """,
        (gsm_id,),
    ).fetchall()

    return {
        "sample": {
            "gse_id": row[0],
            "gsm_id": row[1],
            "sample_title": row[2],
            "status": row[3],
            "parse_quality": row[4],
            "raw_characteristics": row[5],
            "source_name": row[6],
            "growth_condition": row[7],
            "cell_line": row[8],
            "group": row[9],
            "genotype": row[10],
            "treatment": row[11],
            "drug": row[12],
            "concentration": row[13],
            "time": row[14],
            "other_info": row[15],
        },
        "gene_like_rows": [
            {"gene": f[0], "value": f[1], "source": f[2]} for f in features
        ],
    }


def get_stats(conn: sqlite3.Connection) -> dict:
    n_samples = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    n_gse = conn.execute("SELECT COUNT(DISTINCT gse_id) FROM samples WHERE gse_id <> ''").fetchone()[0]
    n_features = conn.execute("SELECT COUNT(*) FROM sample_features").fetchone()[0]
    return {
        "samples": n_samples,
        "gses": n_gse,
        "gene_like_rows": n_features,
        "db_updated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }


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

        if parsed.path.startswith("/api/sample/"):
            gsm = parsed.path.rsplit("/", 1)[-1].upper().strip()
            payload = get_sample_detail(self.conn, gsm)
            if not payload:
                self._send_json({"error": "not found"}, 404)
            else:
                self._send_json(payload)
            return

        if parsed.path in {"/", "/index.html"}:
            self.path = "/index.html"
        return super().do_GET()


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
