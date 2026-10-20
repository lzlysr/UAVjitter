"""只读检索本地 Zotero 的 UAV jitter 分类，保存可复核的全文索引与清单。

不访问公网、不改动 Zotero，不将关键词命中自动判定为论文采用了 UPA。
全文缓存属于本机研究资料，由 data/literature 忽略规则排除在 Git 之外。
"""
import json
import re
import urllib.error
import urllib.request
from pathlib import Path

BASE = "http://127.0.0.1:23119/api/users/0/"
COLLECTION = "TPNHRYML"
OUT = Path(__file__).resolve().parents[1] / "data/literature/zotero_uav_jitter_2026-10-20"


def get(route, raw=False):
    """只使用 Zotero 本地只读 API；超时和未索引状态保留到审计清单。"""
    request = urllib.request.Request(BASE + route, headers={"Zotero-API-Version": "3"})
    with urllib.request.urlopen(request, timeout=20) as response:
        return response.read().decode("utf-8").strip() if raw else json.load(response)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    items = get(f"collections/{COLLECTION}/items/top?limit=100")
    if len(items) == 100:
        raise RuntimeError("分类可能超过一页，须增加分页，不能将当前结果视为全部条目。")
    records = []
    for item in items:
        row = {"key": item["key"], "data": item["data"], "attachments": []}
        for child in get(f"items/{item['key']}/children"):
            data = child["data"]
            if data.get("itemType") != "attachment":
                continue
            attachment = {"key": child["key"], "filename": data.get("filename"),
                          "contentType": data.get("contentType")}
            try:
                fulltext = get(f"items/{child['key']}/fulltext")
                text = fulltext.get("content", "")
                attachment["textLength"] = len(text)
                attachment["indexedPages"] = fulltext.get("indexedPages")
                attachment["totalPages"] = fulltext.get("totalPages")
                attachment["fileUrl"] = get(f"items/{child['key']}/file/view/url", raw=True)
                (OUT / f"{item['key']}_{child['key']}.txt").write_text(text, encoding="utf-8")
                pattern = r".{0,100}(?:uniform\s+(?:planar|rectangular|linear)|\bUPA\b|\bURA\b|\bULA\b|Frobenius|nonnegative|non-negative|least.squares|covariance fitting).{0,130}"
                attachment["hits"] = re.findall(pattern, re.sub(r"\s+", " ", text), re.I)
            except (urllib.error.URLError, ValueError) as error:
                attachment["error"] = str(error)
            row["attachments"].append(attachment)
        records.append(row)
        print(item["key"], item["data"].get("title", ""))
        for attachment in row["attachments"]:
            print("  ", attachment["key"], attachment.get("textLength", 0),
                  attachment.get("error", ""))
            for hit in attachment.get("hits", [])[:8]:
                print("   ", hit)
    (OUT / "inventory.json").write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8")
    print("分类条目数:", len(records), "清单:", OUT / "inventory.json")


if __name__ == "__main__":
    main()
