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
    st.set_page_config(page_title="GSE Gene-Value Query", layout="wide")
    st.title("GSE Gene-Value 查询")
    st.caption("给定 GSE 后，可按 GSM / 样本名 / 条件 / Gene 查询表达值")

    ensure_database()
    conn = get_conn()
    stats = app.get_stats(conn)

    c1, c2, c3 = st.columns(3)
    c1.metric("样本数", stats["samples"])
    c2.metric("GSE 数", stats["gses"])
    c3.metric("已载入表达行", stats["expression_rows"])

    with st.sidebar:
        st.header("筛选")
        gse = st.text_input("GSE（必填）", placeholder="GSE250283").strip().upper()
        gsm = st.text_input("GSM", placeholder="GSM7976778").strip().upper()
        sample = st.text_input("样本名", placeholder="例如 WT_Log 或 DMM00002").strip()
        condition = st.text_input("条件关键词", placeholder="例如 CAD / NaCl / heat shock").strip()
        gene = st.text_input("Gene", placeholder="例如 AAC1").strip()
        limit = st.number_input("结果上限", min_value=10, max_value=5000, value=500, step=50)
        do_search = st.button("查询 Gene-Value", type="primary")

    if not gse:
        st.info("请先输入 GSE")
        return

    if do_search:
        with st.spinner(f"正在准备 {gse} 的表达矩阵..."):
            prep = app.ensure_expression_loaded(conn, gse)

        if not prep.get("ok"):
            st.error(prep.get("message", "表达矩阵载入失败"))
            return

        st.success(f"{gse} 表达矩阵已就绪：{prep.get('row_count', 0)} 行，文件 {prep.get('file_count', 0)} 个")

    params = {
        "gse": gse,
        "gsm": gsm,
        "sample": sample,
        "condition": condition,
        "gene": gene,
        "limit": str(limit),
    }

    rows = app.search_expression(conn, params)

    st.subheader(f"Gene-Value 结果 ({len(rows)})")
    if not rows:
        st.warning("当前条件下没有结果。请先点一次“查询 Gene-Value”进行该 GSE 的矩阵加载。")
        return

    df = pd.DataFrame(rows)
    show_cols = [
        "gse_id",
        "gsm_id",
        "sample_label",
        "sample_title",
        "condition_text",
        "gene",
        "value",
    ]
    st.dataframe(df[show_cols], use_container_width=True, hide_index=True)


if __name__ == "__main__":
    run()
