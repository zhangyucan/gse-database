#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import pandas as pd


NON_SAMPLE_LABELS = {
    "start",
    "end",
    "length",
    "tags",
    "maxheight",
    "summit",
    "percent",
    "accumulativevalue",
    "exonlen",
    "introlen",
    "depth",
    "coverage",
    "fpkm",
    "totalreads",
    "sensereads",
    "antisensereads",
    "intronreads",
    "allreadsfpkm",
    "geneid",
    "efflength",
    "estcounts",
    "tpm",
    "meannormcounts",
    "log2foldchange",
    "seoflog2foldchange",
    "waldstatistic",
    "pvalue",
    "fdr",
    "logfc",
    "logcpm",
    "replicate",
}


def clean(text) -> str:
    s = str(text).strip()
    if s.lower() in {"nan", "none", "null"}:
        return ""
    return s


def norm(text: str) -> str:
    return "".join(ch for ch in clean(text).lower() if ch.isalnum())


def tokens(text: str) -> list:
    s = clean(text).lower()
    out = []
    cur = []
    for ch in s:
        if ch.isalnum():
            cur.append(ch)
        else:
            if cur:
                out.append("".join(cur))
                cur = []
    if cur:
        out.append("".join(cur))
    return out


def expand_tokens(token_list: list) -> set:
    t = set(token_list)
    if "wt" in t:
        t.update({"wild", "type"})
    if "q" in t:
        t.update({"q", "quiescence"})
    if "ds" in t:
        t.update({"ds", "diauxic"})
    return t


def label_is_meta_field(label: str) -> bool:
    n = norm(label)
    if n in NON_SAMPLE_LABELS:
        return True
    if n.endswith("pvalue") or n.endswith("fdr") or n.endswith("logfc"):
        return True
    if re.fullmatch(r"a+b+pvalue", n) or re.fullmatch(r"a+b+fdr", n) or re.fullmatch(r"a+b+logfc", n):
        return True
    return False


def load_whitelist_rules(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(
            columns=["gse_id", "label_pattern", "target_gsm_ids", "enabled", "is_regex", "note"]
        )
    df = pd.read_csv(path, dtype=str).fillna("")
    for c in ["gse_id", "label_pattern", "target_gsm_ids", "enabled", "is_regex", "note"]:
        if c not in df.columns:
            df[c] = ""
    df["gse_id"] = df["gse_id"].map(clean).str.upper()
    df["label_pattern"] = df["label_pattern"].map(clean)
    df["target_gsm_ids"] = df["target_gsm_ids"].map(clean)
    df["enabled"] = df["enabled"].map(clean).replace("", "1")
    df["is_regex"] = df["is_regex"].map(clean).replace("", "0")
    df["note"] = df["note"].map(clean)
    return df[["gse_id", "label_pattern", "target_gsm_ids", "enabled", "is_regex", "note"]]


def find_whitelist_match(gse_id: str, label: str, rules: pd.DataFrame) -> list:
    if rules.empty:
        return []
    sub = rules[(rules["gse_id"] == gse_id) & (rules["enabled"] == "1")]
    if sub.empty:
        return []
    nlabel = norm(label)
    for _, r in sub.iterrows():
        pat = r["label_pattern"]
        if not pat:
            continue
        if r["is_regex"] == "1":
            if re.search(pat, label, flags=re.I):
                return [x.strip().upper() for x in r["target_gsm_ids"].split(";") if x.strip()]
        else:
            if nlabel == norm(pat):
                return [x.strip().upper() for x in r["target_gsm_ids"].split(";") if x.strip()]
    return []


def build_master_mapping(local_csv: Path, geo_csv: Path, out_master: Path) -> pd.DataFrame:
    local = pd.read_csv(local_csv, dtype=str).fillna("")
    for c in local.columns:
        local[c] = local[c].map(clean)
    local["gse_id"] = local["GSE_ID"].str.upper().str.strip()
    local["gsm_id"] = local["GSM_ID"].str.upper().str.strip()
    local["sample_title_local"] = local["Sample_Title"].map(clean)
    for c in ["treatment", "group", "genotype", "growth_condition"]:
        if c not in local.columns:
            local[c] = ""
        local[c] = local[c].map(clean)
    local["condition_local"] = local[["treatment", "group", "genotype", "growth_condition"]].apply(
        lambda r: " | ".join([x for x in r.tolist() if x]),
        axis=1,
    )
    local_base = local[["gse_id", "gsm_id", "sample_title_local", "condition_local"]].copy()
    local_base = local_base[(local_base["gse_id"] != "") & (local_base["gsm_id"] != "")]

    geo = pd.read_csv(geo_csv, dtype=str).fillna("")
    for c in geo.columns:
        geo[c] = geo[c].map(clean)
    geo["gse_id"] = geo["gse_id"].str.upper().str.strip()
    geo["gsm_id"] = geo["gsm_id"].str.upper().str.strip()
    geo = geo[["gse_id", "gsm_id", "sample_title_geo", "condition_inferred"]].rename(
        columns={"condition_inferred": "condition_geo"}
    )

    master = local_base.merge(geo, on=["gse_id", "gsm_id"], how="outer")
    for c in ["sample_title_local", "condition_local", "sample_title_geo", "condition_geo"]:
        master[c] = master[c].map(clean)
    master["mapping_source"] = master.apply(
        lambda r: "LOCAL_META"
        if clean(r["sample_title_local"])
        else ("GEO_FETCH" if clean(r["sample_title_geo"]) else "UNKNOWN"),
        axis=1,
    )
    master["sample_title"] = master.apply(
        lambda r: clean(r["sample_title_local"]) if clean(r["sample_title_local"]) else clean(r["sample_title_geo"]),
        axis=1,
    )
    master["condition_text"] = master.apply(
        lambda r: clean(r["condition_local"]) if clean(r["condition_local"]) else clean(r["condition_geo"]),
        axis=1,
    )
    master = master[
        [
            "gse_id",
            "gsm_id",
            "sample_title",
            "condition_text",
            "mapping_source",
            "sample_title_local",
            "condition_local",
            "sample_title_geo",
            "condition_geo",
        ]
    ].drop_duplicates(["gse_id", "gsm_id"])
    master = master.sort_values(["gse_id", "gsm_id"])
    out_master.parent.mkdir(parents=True, exist_ok=True)
    master.to_csv(out_master, index=False)
    return master


def build_label_mapping(
    labels_csv: Path,
    master: pd.DataFrame,
    whitelist_rules: pd.DataFrame,
    out_label_map: Path,
    out_summary: Path,
):
    labels = pd.read_csv(labels_csv, dtype=str).fillna("")
    for c in labels.columns:
        labels[c] = labels[c].map(clean)
    labels["gse_id"] = labels["gse_id"].str.upper().str.strip()
    labels["gsm_in_label"] = labels["gsm_in_label"].str.upper().str.strip()

    by_gse = {g: sub.to_dict("records") for g, sub in master.groupby("gse_id")}
    mapped = []
    for _, r in labels.iterrows():
        gse_id = r["gse_id"]
        label = r["sample_label"]
        gsm_in = r["gsm_in_label"]
        candidates = by_gse.get(gse_id, [])
        status = "UNMATCHED"
        matched = []
        reason = ""

        if gsm_in:
            hit = [x for x in candidates if x["gsm_id"] == gsm_in]
            if hit:
                matched = hit
                status = "MATCH_BY_GSM"
                reason = "GSM in label"
            else:
                status = "GSM_NOT_IN_MAPPING"
                reason = "GSM not found in master"
        else:
            if label_is_meta_field(label):
                status = "IGNORED_META_FIELD"
                reason = "meta/stat field"
            else:
                wl_gsms = find_whitelist_match(gse_id, label, whitelist_rules)
                if wl_gsms:
                    matched = [x for x in candidates if x["gsm_id"] in set(wl_gsms)]
                    if matched:
                        status = "MATCH_BY_WHITELIST"
                        reason = "whitelist"

            if status == "UNMATCHED":
                lt = expand_tokens(tokens(label))
                ln = norm(label)
                ll = label.lower().replace("-", "_")
                scored = []
                for x in candidates:
                    st = clean(x.get("sample_title", ""))
                    ct = clean(x.get("condition_text", ""))
                    stl = st.lower()
                    tt = expand_tokens(tokens(st))
                    ctok = expand_tokens(tokens(ct))
                    sn = norm(st)
                    score = len(lt & tt) + len(lt & ctok)
                    if ln and sn and (ln in sn or sn in ln):
                        score += 2
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
                        scored.append((score, x))
                if scored:
                    best = max(s for s, _ in scored)
                    matched = [x for s, x in scored if s == best]
                    status = "MATCH_BY_HEURISTIC"
                    reason = f"score={best}"
                else:
                    status = "UNMATCHED_NO_GSM"
                    reason = "no overlap"

        mapped.append(
            {
                "gse_id": gse_id,
                "source_file": r["source_file"],
                "sample_label": label,
                "gsm_in_label": gsm_in,
                "match_status": status,
                "matched_gsm_ids": ";".join(sorted({x["gsm_id"] for x in matched})),
                "matched_sample_titles": ";".join(
                    sorted({clean(x.get("sample_title", "")) for x in matched if clean(x.get("sample_title", ""))})
                ),
                "matched_conditions": ";".join(
                    sorted({clean(x.get("condition_text", "")) for x in matched if clean(x.get("condition_text", ""))})
                ),
                "matched_source": ";".join(
                    sorted({clean(x.get("mapping_source", "")) for x in matched if clean(x.get("mapping_source", ""))})
                ),
                "match_count": len(set(x["gsm_id"] for x in matched)),
                "match_reason": reason,
            }
        )

    map_df = pd.DataFrame(mapped)
    out_label_map.parent.mkdir(parents=True, exist_ok=True)
    map_df.to_csv(out_label_map, index=False)

    sum_df = map_df.groupby("gse_id", as_index=False).agg(
        label_count=("sample_label", "count"),
        sample_like_count=("match_status", lambda s: int((s != "IGNORED_META_FIELD").sum())),
        by_gsm=("match_status", lambda s: int((s == "MATCH_BY_GSM").sum())),
        by_whitelist=("match_status", lambda s: int((s == "MATCH_BY_WHITELIST").sum())),
        by_heuristic=("match_status", lambda s: int((s == "MATCH_BY_HEURISTIC").sum())),
        ignored_meta=("match_status", lambda s: int((s == "IGNORED_META_FIELD").sum())),
        unmatched=("match_status", lambda s: int((s.str.startswith("UNMATCHED")).sum())),
        gsm_not_in_mapping=("match_status", lambda s: int((s == "GSM_NOT_IN_MAPPING").sum())),
    )
    sum_df["mapped_any"] = sum_df["by_gsm"] + sum_df["by_whitelist"] + sum_df["by_heuristic"]
    sum_df["mapping_status"] = sum_df.apply(
        lambda r: "NONE"
        if int(r["sample_like_count"]) == 0 or int(r["mapped_any"]) == 0
        else ("FULL" if int(r["mapped_any"]) >= int(r["sample_like_count"]) else "PARTIAL"),
        axis=1,
    )
    sum_df.to_csv(out_summary, index=False)
    return map_df, sum_df


def update_whitelist_template(
    whitelist_path: Path,
    map_df: pd.DataFrame,
    summary_df: pd.DataFrame,
) -> pd.DataFrame:
    rules = load_whitelist_rules(whitelist_path)
    existing_keys = set(
        (clean(r["gse_id"]).upper(), norm(r["label_pattern"]))
        for _, r in rules.iterrows()
        if clean(r["gse_id"]) and clean(r["label_pattern"])
    )

    need_gse = summary_df[summary_df["mapping_status"] == "NONE"]["gse_id"].tolist()
    add_rows = []
    for gse_id in need_gse:
        sub = map_df[(map_df["gse_id"] == gse_id) & (map_df["match_status"] == "UNMATCHED_NO_GSM")]
        for label in sorted(set(sub["sample_label"].tolist())):
            key = (gse_id, norm(label))
            if key in existing_keys:
                continue
            add_rows.append(
                {
                    "gse_id": gse_id,
                    "label_pattern": label,
                    "target_gsm_ids": "",
                    "enabled": "0",
                    "is_regex": "0",
                    "note": "AUTO_TODO",
                }
            )
            existing_keys.add(key)

    if add_rows:
        rules = pd.concat([rules, pd.DataFrame(add_rows)], ignore_index=True)
        rules = rules.sort_values(["gse_id", "enabled", "label_pattern"], ascending=[True, False, True])
        whitelist_path.parent.mkdir(parents=True, exist_ok=True)
        rules.to_csv(whitelist_path, index=False)
    return rules


def main():
    parser = argparse.ArgumentParser(description="Build GSE-GSM-condition mapping tables with whitelist rules")
    parser.add_argument("--local-csv", default="/data1/zyc/lu2026/web_prototype/data/extraction_batch_result.csv")
    parser.add_argument("--geo-csv", default="/data1/zyc/lu2026/SS_Bulk/gse_gsm_condition_from_geo.csv")
    parser.add_argument("--labels-csv", default="/data1/zyc/lu2026/SS_Bulk/gse_sample_labels_extracted.csv")
    parser.add_argument("--whitelist", default="/data1/zyc/lu2026/SS_Bulk/mapping_whitelist_rules.csv")
    parser.add_argument("--out-master", default="/data1/zyc/lu2026/SS_Bulk/gse_gsm_condition_master_mapping.csv")
    parser.add_argument("--out-label-map", default="/data1/zyc/lu2026/SS_Bulk/gse_expression_label_mapping.csv")
    parser.add_argument("--out-summary", default="/data1/zyc/lu2026/SS_Bulk/gse_mapping_alignment_summary.csv")
    args = parser.parse_args()

    local_csv = Path(args.local_csv)
    geo_csv = Path(args.geo_csv)
    labels_csv = Path(args.labels_csv)
    whitelist_path = Path(args.whitelist)
    out_master = Path(args.out_master)
    out_label_map = Path(args.out_label_map)
    out_summary = Path(args.out_summary)

    whitelist = load_whitelist_rules(whitelist_path)
    master = build_master_mapping(local_csv, geo_csv, out_master)
    map_df, summary_df = build_label_mapping(labels_csv, master, whitelist, out_label_map, out_summary)
    rules_after = update_whitelist_template(whitelist_path, map_df, summary_df)

    print(f"[DONE] master: {out_master} rows={len(master)} gse={master['gse_id'].nunique()}")
    print(
        f"[DONE] label_map: {out_label_map} rows={len(map_df)} mapped="
        f"{int((map_df['match_status'].str.startswith('MATCH_')).sum())}"
    )
    print(
        f"[DONE] summary: {out_summary} gse={len(summary_df)} "
        f"FULL={int((summary_df['mapping_status']=='FULL').sum())} "
        f"PARTIAL={int((summary_df['mapping_status']=='PARTIAL').sum())} "
        f"NONE={int((summary_df['mapping_status']=='NONE').sum())}"
    )
    print(f"[DONE] whitelist: {whitelist_path} rows={len(rules_after)}")


if __name__ == "__main__":
    main()
