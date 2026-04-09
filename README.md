# GEO Query Prototype (Web + Streamlit)

这个项目是你当前 `GSE-GSM-条件` 数据的第一版查询站，支持：

1. Web API + 前端查询页（`app.py` + `static/`）
2. Streamlit 查询页（`streamlit_app.py`）

数据源使用：`data/extraction_batch_result.csv`

## 本地运行

```bash
cd /data1/zyc/lu2026/web_prototype
python3 app.py --port 8765
```

打开：`http://127.0.0.1:8765`

## 本地运行 Streamlit

```bash
cd /data1/zyc/lu2026/web_prototype
python3 -m pip install -r requirements.txt
streamlit run streamlit_app.py
```

## GitHub 上传

```bash
cd /data1/zyc/lu2026/web_prototype
git init
git add .
git commit -m "feat: initial GEO query prototype with streamlit"
git branch -M main
git remote add origin <你的仓库URL>
git push -u origin main
```

## Streamlit Cloud 关联

1. 打开 `https://share.streamlit.io/`
2. 选择你的 GitHub 仓库
3. `Main file path` 填：`streamlit_app.py`
4. Python 依赖自动读取 `requirements.txt`
5. 点击 Deploy

## 说明

- `Gene/Value` 视图目前是基于样本条件字段映射，不是完整表达矩阵。
- 下一步可把真实表达矩阵按 `GSE+GSM` 合并进来，替换为真正的 gene-level value。
