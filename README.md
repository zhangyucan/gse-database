# GSE Gene-Value Query (Demo + Full)

## Demo 版（推荐给同学试用）

仓库默认是轻量 demo 数据，只含 3 个 GSE：
- `GSE237235`
- `GSE67149`
- `GSE110818`

直接运行：

```bash
cd /data1/zyc/lu2026/web_prototype
python3 -m pip install -r requirements.txt
python3 -m streamlit run streamlit_app.py --server.port 8501 --server.address 0.0.0.0
```

浏览器打开：`http://127.0.0.1:8501`

## Demo 验证示例

- `GSE67149` + `GSM1640265`
- `GSE237235` + `GSM7597583`

## 全量模式（你本机）

如果你要接本地全量数据（`/data1/zyc/lu2026/SS_Bulk`），可运行：

```bash
cd /data1/zyc/lu2026/web_prototype
python3 rebuild_all.py --global-rebuild
```

这会重建：
- 映射表（GSE/GSM/条件）
- 浏览库 `geo_samples.db`
- 全量表达库 `data/ss_bulk_expression.db`

## 环境变量（可选）

- `GEO_CSV_PATH`：覆盖样本元信息 CSV
- `SS_BULK_DIR`：覆盖表达文件根目录
- `GLOBAL_EXPR_DB`：覆盖全局表达库路径
- `MASTER_MAPPING_CSV`：覆盖主映射表
- `LABEL_MAPPING_CSV`：覆盖标签映射表
- `MAPPING_SUMMARY_CSV`：覆盖映射状态汇总

