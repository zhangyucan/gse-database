#!/usr/bin/env python3
import argparse
import sqlite3
import subprocess
from pathlib import Path

import app


def run_cmd(cmd: list[str]) -> None:
    print("[RUN]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="One-click rebuild: mapping + browser DB + global expression DB")
    parser.add_argument("--skip-global-sql", action="store_true", help="Skip rebuilding ss_bulk_expression.db")
    parser.add_argument("--global-max-files", type=int, default=0, help="Optional max files for global SQL builder")
    parser.add_argument("--global-rebuild", action="store_true", help="Force truncate and rebuild global SQL")
    args = parser.parse_args()

    base = Path(__file__).resolve().parent
    mapping_builder = base / "build_mapping_tables.py"
    global_builder = base / "build_ss_bulk_sql.py"
    browser_db = base / "geo_samples.db"

    # 1) refresh mapping outputs + status summary
    run_cmd(["python3", str(mapping_builder)])

    # 2) rebuild browser DB (samples + mapping_samples + mapping_status)
    csv_path = app.resolve_csv_path()
    print("[RUN] import_csv_to_db ->", browser_db)
    app.import_csv_to_db(csv_path, browser_db)

    conn = sqlite3.connect(browser_db)
    stats = app.get_stats(conn)
    conn.close()
    print("[DONE] browser db stats:", stats)

    # 3) optional global expression SQL refresh
    if not args.skip_global_sql:
        cmd = [
            "python3",
            str(global_builder),
            "--ss-bulk",
            str(app.resolve_ss_bulk_dir()),
            "--out",
            str(app.resolve_global_expression_db()),
        ]
        if args.global_max_files > 0:
            cmd.extend(["--max-files", str(args.global_max_files)])
        if args.global_rebuild:
            cmd.append("--rebuild")
        run_cmd(cmd)

    print("[ALL DONE]")


if __name__ == "__main__":
    main()
