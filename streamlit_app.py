#!/usr/bin/env python3
import sqlite3
from pathlib import Path

import pandas as pd
import streamlit as st

import app

BASE_DIR = Path(__file__).resolve().parent
DB_PATH = BASE_DIR / "geo_samples.db"


@st.cache_resource(show_spinner=False)
def get_conn() -> sqlite3.Connection:
    return sqlite3.connect(DB_PATH, check_same_thread=False)


def ensure_database() -> None:
    csv_path = app.resolve_csv_path()
    if not csv_path.exists():
        st.error(f"CSV not found: {csv_path}")
        st.stop()

    needs_rebuild = (not DB_PATH.exists()) or (DB_PATH.stat().st_mtime < csv_path.stat().st_mtime)
    if needs_rebuild:
        app.import_csv_to_db(csv_path, DB_PATH)


def run() -> None:
    st.set_page_config(page_title="GEO Query Prototype", layout="wide")
    st.title("GSE · GSM 条件数据库")
    st.caption("基于现有 extraction_batch_result.csv 的 Streamlit 查询页")

    ensure_database()
    conn = get_conn()
    stats = app.get_stats(conn)

    c1, c2, c3 = st.columns(3)
    c1.metric("样本数", stats["samples"])
    c2.metric("GSE 数", stats["gses"])
    c3.metric("Gene/Value 行数", stats["gene_like_rows"])

    with st.sidebar:
        st.header("查询条件")
        gse = st.text_input("GSE", placeholder="GSE250283").strip().upper()
        gsm = st.text_input("GSM", placeholder="GSM7976778").strip().upper()
        q = st.text_input("关键词", placeholder="DMM00002 / CAD / treatment").strip()
        limit = st.number_input("结果条数", min_value=1, max_value=500, value=100, step=10)

    params = {"gse": gse, "gsm": gsm, "q": q, "limit": str(limit)}
    rows = app.search_samples(conn, params)

    st.subheader(f"查询结果 ({len(rows)})")
    if not rows:
        st.info("没有匹配结果")
        return

    df = pd.DataFrame(rows)
    display_cols = [
        "gse_id",
        "gsm_id",
        "sample_title",
        "treatment",
        "genotype",
        "feature_count",
        "parse_quality",
        "status",
    ]
    st.dataframe(df[display_cols], use_container_width=True, hide_index=True)

    default_gsm = rows[0]["gsm_id"]
    gsm_pick = st.selectbox("查看样本详情", options=[r["gsm_id"] for r in rows], index=0)
    if not gsm_pick:
        gsm_pick = default_gsm

    detail = app.get_sample_detail(conn, gsm_pick)
    if not detail:
        st.warning("未找到样本详情")
        return

    st.subheader(f"样本详情: {gsm_pick}")
    sample_df = pd.DataFrame([detail["sample"]])
    st.dataframe(sample_df, use_container_width=True, hide_index=True)

    st.subheader("Gene/Value 结构视图（字段映射版）")
    feature_df = pd.DataFrame(detail["gene_like_rows"])
    if feature_df.empty:
        st.info("暂无明细")
    else:
        st.dataframe(feature_df, use_container_width=True, hide_index=True)


if __name__ == "__main__":
    run()
