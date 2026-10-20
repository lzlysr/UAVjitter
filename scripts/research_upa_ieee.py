"""查询Crossref，保存可审计的检索记录；命中不等于满足UPA纳入条件。"""
import argparse
import concurrent.futures
import datetime
import json
import pathlib
import re
import urllib.parse
import urllib.request

QUERIES = [
    "planar array adaptive beamforming covariance reconstruction",
    "uniform rectangular array interference covariance reconstruction",
    "two dimensional angular integration interference noise beamforming",
    "uniform planar robust beamforming covariance",
    "planar interference covariance sparse reconstruction",
    "planar beamforming spatial spectrum matching",
]


def query(term):
    """官方元数据接口；独立记录失败，不用空结果伪装检索成功。"""
    parameters = urllib.parse.urlencode({"query": term, "rows": 100,
        "filter": "from-pub-date:2019-01-01,until-pub-date:2026-10-20"})
    url = "https://api.crossref.org/works?" + parameters
    request = urllib.request.Request(url, headers={"User-Agent": "UPA-INCM-LiteratureAudit/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = json.load(response)["message"]
        return {"query": term, "url": url, "status": "ok", "total": payload["total-results"], "items": payload["items"]}
    except Exception as error:
        return {"query": term, "url": url, "status": "failed", "error": str(error), "items": []}


def main():
    parser = argparse.ArgumentParser(description="UPA/IEEE文献候选元数据检索，不作全文结论")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        searches = list(executor.map(query, QUERIES))
    merged = {}
    for search in searches:
        for item in search["items"]:
            doi = item.get("DOI", "").strip().lower()
            if not doi:
                continue
            if doi not in merged:
                merged[doi] = {"metadata": item, "queries": []}
            merged[doi]["queries"].append(search["query"])
    candidates = []
    for doi, item in merged.items():
        metadata = item["metadata"]
        title = " ".join(metadata.get("title", []))
        journal = " ".join(metadata.get("container-title", []))
        abstract = re.sub("<[^>]+>", " ", metadata.get("abstract", ""))
        text = (title + " " + abstract).lower()
        if doi.startswith("10.1109/") and "ieee access" not in journal.lower() and "beamform" in text:
            score = sum(weight for pattern, weight in [("planar", 5), ("rectangular", 5), ("reconstruct", 3),
                ("interference", 2), ("covariance", 2), ("two-dimensional", 1)] if pattern in text)
            candidates.append({"doi": doi, "title": title, "journal": journal, "score": score,
                "abstract": abstract, "queries": item["queries"]})
    candidates.sort(key=lambda row: (-row["score"], row["title"]))
    output = pathlib.Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    audit = {"created_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "source": "Crossref", "searches": searches, "unique_records": len(merged), "candidates": candidates}
    (output / "crossref_search_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"search_status": [{"query": q["query"], "status": q["status"]} for q in searches],
        "unique_records": len(merged), "candidates": [{k: r[k] for k in ("doi", "title", "journal", "score")} for r in candidates]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
