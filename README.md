# GSE Gene-Value Query (Streamlit)

## 你当前遇到的报错

`no expression-like files found in /mount/src/SS_Bulk`

原因是 Streamlit Cloud 没有你本地的 `/data1/zyc/lu2026/SS_Bulk` 目录。

## 已修复

1. 查询会优先读取全局表达 SQL（`data/ss_bulk_expression.db`）。
2. 如果全局 SQL 不存在，再尝试扫描 `SS_Bulk`。
3. 扫描失败时会回退到项目内 `data/expression` demo 文件。

## 把 SS_Bulk 全量做成一个 SQL

在本机执行：

```bash
cd /data1/zyc/lu2026/web_prototype
python3 build_ss_bulk_sql.py \
  --ss-bulk /data1/zyc/lu2026/SS_Bulk \
  --out /data1/zyc/lu2026/web_prototype/data/ss_bulk_expression.db \
  --rebuild
```

说明：
- 这个 SQL 会把所有识别为表达矩阵的表格统一写入 `expression_raw`。
- 表结构字段：`gse_id, sample_label, gene, value, source_file`。

## 查询逻辑（与你需求一致）

1. 输入 `GSE`。
2. 先看该 `GSE` 的基本情况表（GSM + 样本名 + 条件）。
3. 锁定 `GSM` 后，直接展示该样本的整张 `Gene-Value` 表。

## 本地运行

```bash
cd /data1/zyc/lu2026/web_prototype
python3 -m pip install -r requirements.txt
streamlit run streamlit_app.py
```
