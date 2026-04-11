#!/usr/bin/env python3
import sqlite3
from pathlib import Path

import pandas as pd
import streamlit as st

import app

BASE_DIR = Path(__file__).resolve().parent
DB_PATH = BASE_DIR / "geo_samples.db"
SOURCE_INDEX_PATH = BASE_DIR.parent / "SS_Bulk" / "gse_table_source_index.csv"


@st.cache_resource(show_spinner=False)
def get_conn() -> sqlite3.Connection:
    return sqlite3.connect(DB_PATH, check_same_thread=False)


@st.cache_data(show_spinner=False)
def load_source_index() -> pd.DataFrame:
    if not SOURCE_INDEX_PATH.exists():
        return pd.DataFrame(columns=["gse_id", "table_source", "chosen_table_file", "final_status"])
    try:
        df = pd.read_csv(SOURCE_INDEX_PATH)
    except Exception:
        return pd.DataFrame(columns=["gse_id", "table_source", "chosen_table_file", "final_status"])
    if "gse_id" not in df.columns:
        return pd.DataFrame(columns=["gse_id", "table_source", "chosen_table_file", "final_status"])
    df["gse_id"] = df["gse_id"].astype(str).str.upper().str.strip()
    if "table_source" not in df.columns:
        df["table_source"] = "NONE"
    if "chosen_table_file" not in df.columns:
        df["chosen_table_file"] = ""
    if "final_status" not in df.columns:
        df["final_status"] = ""
    return df


def ensure_database() -> None:
    csv_path = app.resolve_csv_path()
    if not csv_path.exists():
        st.error(f"CSV not found: {csv_path}")
        st.stop()

    dep_paths = [
        csv_path,
        app.resolve_master_mapping_csv(),
        app.resolve_label_mapping_csv(),
        app.resolve_mapping_summary_csv(),
    ]
    newest_dep_mtime = 0.0
    for p in dep_paths:
        if p.exists():
            newest_dep_mtime = max(newest_dep_mtime, p.stat().st_mtime)

    needs_rebuild = (not DB_PATH.exists()) or (DB_PATH.stat().st_mtime < newest_dep_mtime)
    if needs_rebuild:
        app.import_csv_to_db(csv_path, DB_PATH)


def run() -> None:
    st.set_page_config(page_title="GSE Gene-Value Query", layout="wide")
    st.title("GSE Gene-Value 查询")
    st.caption("给定一个 GSE，先看 GSM 对应基本表，再锁定 GSM 查看整张 Gene-Value 表")

    ensure_database()
    conn = get_conn()
    stats = app.get_stats(conn)
    source_df = load_source_index()
    source_map = {}
    chosen_map = {}
    if not source_df.empty:
        source_map = dict(zip(source_df["gse_id"], source_df["table_source"]))
        chosen_map = dict(zip(source_df["gse_id"], source_df["chosen_table_file"]))

    c1, c2, c3 = st.columns(3)
    c1.metric("样本数", stats["samples"])
    c2.metric("GSE 数", stats["gses"])
    c3.metric("已载入表达行", stats["expression_rows"])

    with st.sidebar:
        st.header("筛选")
        gse = st.text_input("GSE", placeholder="GSE250283").strip().upper()
        gsm = st.text_input("GSM", placeholder="GSM7976778").strip().upper()
        sample = st.text_input("样本名", placeholder="AR13_1 或 DMM00002").strip()
        condition = st.text_input("条件关键词", placeholder="CAD / NaCl / heat shock").strip()
        source_filter = st.selectbox("来源", options=["ALL", "MANUAL", "AUTO", "NONE"], index=0)
        do_search = st.button("查询", type="primary")

    if not gse:
        st.subheader("样本总览大表")
        overview_limit = st.number_input("总览行数", min_value=100, max_value=5000, value=1000, step=100)
        overview_params = {"limit": str(overview_limit)}
        if gsm:
            overview_params["gsm"] = gsm
        if sample:
            overview_params["q"] = sample
        overview_rows = app.search_samples(conn, overview_params)
        if overview_rows:
            odf = pd.DataFrame(overview_rows)
            odf["table_source"] = odf["gse_id"].map(source_map).fillna("NONE")
            if source_filter != "ALL":
                odf = odf[odf["table_source"] == source_filter]
            odf = odf.copy()
            show_cols = [
                "gse_id",
                "table_source",
                "mapping_status",
                "gsm_id",
                "sample_title",
                "treatment",
                "genotype",
                "raw_characteristics",
                "feature_count",
            ]
            st.dataframe(odf[show_cols], use_container_width=True, height=420, hide_index=True)
        else:
            st.info("总览没有匹配结果")
        st.info("输入 GSE 后，将用该 GSE 基本情况表替代总览大表。")
        return

    st.subheader(f"{gse} 基本情况表")
    gse_source = source_map.get(gse, "NONE")
    chosen_file = chosen_map.get(gse, "")
    gse_map_status = app.get_gse_mapping_status(conn, gse)
    st.caption(f"来源: `{gse_source}`")
    st.caption(
        f"映射状态: `{gse_map_status.get('mapping_status','UNKNOWN')}` "
        f"(mapped={gse_map_status.get('mapped_any', 0)}, unmatched={gse_map_status.get('unmatched', 0)})"
    )
    if chosen_file:
        st.caption(f"优先表: `{chosen_file}`")

    gse_rows = app.get_gse_basic_table(conn, gse)
    if gse_rows:
        gse_df = pd.DataFrame(gse_rows)
        st.dataframe(
            gse_df[["gse_id", "gsm_id", "sample_name", "condition"]],
            use_container_width=True,
            height=320,
            hide_index=True,
        )

        options = [f"{x['gsm_id']} | {x['sample_name']}" for x in gse_rows]
        picked = st.selectbox("锁定一个 GSM（推荐）", options=options, index=0)
        picked_gsm = picked.split("|", 1)[0].strip().upper()
        if not gsm:
            gsm = picked_gsm
    else:
        st.warning(f"{gse} 在样本元信息中没有记录。")
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
        "gene": "",
        "limit": "10000",
    }

    limit = st.number_input("结果上限", min_value=100, max_value=50000, value=10000, step=100)
    params["limit"] = str(limit)
    rows = app.search_expression(conn, params)

    if gsm:
        sample_name = ""
        sample_condition = ""
        for x in gse_rows:
            if x["gsm_id"] == gsm:
                sample_name = x["sample_name"]
                sample_condition = x["condition"]
                break
        st.subheader(f"{gsm}|{sample_name}|{sample_condition}  Gene-Value 表")
    else:
        st.subheader(f"Gene-Value 结果 ({len(rows)})")

    if not rows:
        st.warning("当前条件下没有结果。请先点一次“查询 Gene-Value”进行该 GSE 的矩阵加载。")
        return

    df = pd.DataFrame(rows)
    if gsm:
        gv = df[["gene", "value"]].dropna().copy()
        gv = gv.sort_values("gene")
        # 锁定 GSM 后，按 gene 去重，保留一个值用于样本级展示
        gv = gv.drop_duplicates(subset=["gene"], keep="first")
        st.dataframe(gv, use_container_width=True, height=520, hide_index=True)
    else:
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
