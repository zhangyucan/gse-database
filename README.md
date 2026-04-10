# GSE Gene-Value Query (Streamlit)

## 目标

给定一个 `GSE`，可以用 `GSM / 样本名 / 条件 / Gene` 检索对应的 gene-value 信息。

## 数据来源

1. 样本元信息：`data/extraction_batch_result.csv`
2. 表达矩阵：`/data1/zyc/lu2026/SS_Bulk` 下与 GSE 对应的表格文件（如 `xlsx/csv/tsv/txt`）

## 已实现

1. 不再在前端展示 `parse_quality` 和 `status`。
2. 自动修复样本名乱码（如 `cth1Îcth2Î` -> `cth1Δcth2Δ`）。
3. 新增表达矩阵加载与检索：
   - 首次查询某个 GSE 时，自动扫描 `SS_Bulk` 中该 GSE 的表达矩阵文件。
   - 解析第一列 gene 和后续样本列，入库到 `expression_values`。
   - 支持按 `GSM/样本名/条件/Gene` 查询 gene-value。

## 本地运行

```bash
cd /data1/zyc/lu2026/web_prototype
python3 -m pip install -r requirements.txt
streamlit run streamlit_app.py
```

## Streamlit Cloud 部署

- Repository: `zhangyucan/gse-database`
- Branch: `main`
- Main file path: `streamlit_app.py`

