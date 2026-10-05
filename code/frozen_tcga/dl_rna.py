#!/usr/bin/env python3
"""Download TCGA-BRCA RNA-seq V2 RSEM from cBioPortal API in batches."""
import json, sys, time, urllib.request

BASE = "https://www.cbioportal.org/api"
OUT = sys.argv[1] if len(sys.argv) > 1 else "brca_rna_rsem.tsv"

def get(path):
    req = urllib.request.Request(BASE + path, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)

def post(path, payload):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(BASE + path, data=data,
                                 headers={"Content-Type": "application/json",
                                          "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.load(r)

print("fetching gene list...", flush=True)
genes = get("/genes?pageSize=40000&projection=ID")
# keep real genes (positive entrez), drop phosphoprotein pseudo-entries
genes = [g for g in genes if g["entrezGeneId"] > 0]
entrez = [g["entrezGeneId"] for g in genes]
id2sym = {g["entrezGeneId"]: g["hugoGeneSymbol"] for g in genes}
print(f"{len(entrez)} genes", flush=True)

BATCH = 800
# sampleId -> {entrez: value}
mat = {}
order = []
for i in range(0, len(entrez), BATCH):
    chunk = entrez[i:i + BATCH]
    for attempt in range(3):
        try:
            rows = post("/molecular-data/fetch", {
                "entrezGeneIds": chunk,
                "molecularProfileIds": ["brca_tcga_rna_seq_v2_mrna"],
                "sampleListId": "brca_tcga_rna_seq_v2_mrna"})
            break
        except Exception as e:
            print(f"  batch {i//BATCH} attempt {attempt+1} failed: {e}", flush=True)
            time.sleep(5)
    else:
        print(f"FAILED batch starting at {i}", flush=True); sys.exit(1)
    for row in rows:
        sid = row["sampleId"]; eid = row["entrezGeneId"]
        if sid not in mat:
            mat[sid] = {}
            order.append(sid)
        try:
            mat[sid][eid] = float(row["value"])
        except (TypeError, ValueError):
            pass
    print(f"  batch {i//BATCH + 1}/{(len(entrez)+BATCH-1)//BATCH}: {len(order)} samples", flush=True)

print(f"writing {OUT}: {len(order)} samples x {len(entrez)} genes", flush=True)
with open(OUT, "w") as f:
    f.write("sample\t" + "\t".join(id2sym[e] for e in entrez) + "\n")
    for sid in order:
        d = mat[sid]
        f.write(sid + "\t" + "\t".join(
            ("%.4g" % d[e]) if e in d else "NA" for e in entrez) + "\n")
print("done", flush=True)
