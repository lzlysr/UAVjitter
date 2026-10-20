"""汇总本轮 Tavily 原始检索；网址去重不等于论文去重，更不等于合格文献。"""
import concurrent.futures
import csv
import json
import re
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "data/literature/tavily_upa_ieee_2026-10-20"
DOIS = ["10.1109/PAST49659.2022.9975052", "10.1109/JSEN.2020.3030043",
        "10.1109/TVT.2023.3288611", "10.1109/ICSIDP47821.2019.9173458"]


def metadata(doi):
    """通过官方 Crossref API 交叉核对选定文献题录，不将题录当作全文证据。"""
    try:
        with urllib.request.urlopen("https://api.crossref.org/works/" + doi, timeout=25) as response:
            record = json.load(response)["message"]
        keys = ["DOI", "title", "author", "container-title", "published", "volume", "issue", "page", "URL"]
        return {key: record.get(key) for key in keys}
    except Exception as error:
        return {"DOI": doi, "error": str(error)}


def main():
    requests = []
    hits = []
    for file in sorted(ROOT.glob("*.json")):
        record = json.loads(file.read_text(encoding="utf-8"))
        if "result" not in record:
            continue
        text = "\n".join(block.get("text", "") for block in record["result"].get("content", []))
        limited = bool(re.search("hourly keyless Tavily limit|Retry after:", text, re.I))
        requests.append({"record": file.name, "mode": record["mode"], "query": record["query"],
                         "status": "limited" if limited else "returned",
                         "isError": record["result"].get("isError", False)})
        if limited:
            continue
        for chunk in re.split(r"(?=^Title: )", text, flags=re.M):
            title = re.search(r"^Title: (.+)$", chunk, re.M)
            url = re.search(r"^URL: (.+)$", chunk, re.M)
            if not title or not url:
                continue
            snippet = chunk.split("Raw Content:")[0]
            hits.append({"title": title[1], "url": url[1].rstrip("/"),
                         "query": record["query"], "record": file.name,
                         "snippet": snippet})
    unique = {hit["url"]: hit for hit in hits}
    summary = {"requests": requests, "returnedRequestCount": sum(x["status"] == "returned" for x in requests),
               "returnedSearchCount": sum(x["status"] == "returned" and x["mode"] != "extract" for x in requests),
               "returnedExtractCount": sum(x["status"] == "returned" and x["mode"] == "extract" for x in requests),
               "limitedRequestCount": sum(x["status"] == "limited" for x in requests),
               "hitCountIncludingExtract": len(hits), "uniqueUrlCount": len(unique),
               "warning": "计数包含误命中、重复论文、非论文网页，不能当作合格论文数量。"}
    (ROOT / "audit_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    with (ROOT / "unique_urls.tsv").open("w", encoding="utf-8-sig", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=["title", "url", "query", "record", "snippet"], delimiter="\t")
        writer.writeheader()
        writer.writerows(unique.values())
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        verified = list(pool.map(metadata, DOIS))
    (ROOT / "crossref_verified.json").write_text(json.dumps(verified, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "requests"}, ensure_ascii=False, indent=2))
    print(json.dumps(verified, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
