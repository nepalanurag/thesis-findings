#!/usr/bin/env python3
"""Download TCGA-BRCA clinical data (ER/PR/HER2 IHC etc.) from cBioPortal API."""
import json, sys, urllib.request

BASE = "https://www.cbioportal.org/api"

def get(path):
    req = urllib.request.Request(BASE + path, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)

def post(path, payload, tries=4):
    import time
    data = json.dumps(payload).encode()
    for t in range(tries):
        try:
            req = urllib.request.Request(BASE + path, data=data,
                                         headers={"Content-Type": "application/json",
                                                  "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return json.load(r)
        except Exception as e:
            print(f"  retry {t+1}/{tries}: {e}", flush=True)
            time.sleep(3)
    raise RuntimeError("fetch failed")

print("patients...", flush=True)
patients = get("/studies/brca_tcga/patients?pageSize=2000&projection=ID")
pids = [p["patientId"] for p in patients]
print(len(pids), "patients", flush=True)

print("patient clinical (batched)...", flush=True)
pat, samp = [], []
B = 200
for i in range(0, len(pids), B):
    chunk = pids[i:i + B]
    pat += post("/clinical-data/fetch", {
        "clinicalDataType": "PATIENT", "studyId": "brca_tcga",
        "patientIds": chunk, "projection": "SUMMARY"})
    samp += post("/clinical-data/fetch", {
        "clinicalDataType": "SAMPLE", "studyId": "brca_tcga",
        "patientIds": chunk, "projection": "SUMMARY"})
    print(f"  {i//B + 1}/{(len(pids)+B-1)//B}", flush=True)

def to_tsv(rows, out):
    keys = sorted({r["clinicalAttributeId"] for r in rows})
    idx = {}
    for r in rows:
        k = (r.get("patientId"), r.get("sampleId"))
        idx.setdefault(k, {})[r["clinicalAttributeId"]] = r["value"]
    with open(out, "w") as f:
        f.write("patientId\tsampleId\t" + "\t".join(keys) + "\n")
        for (pid, sid), d in sorted(idx.items(), key=lambda x: (x[0][0], str(x[0][1]))):
            f.write("\t".join([pid or "", sid or ""] +
                               [d.get(k, "") for k in keys]) + "\n")
    print(f"wrote {out}: {len(idx)} rows x {len(keys)} attrs", flush=True)

to_tsv(pat, "brca_clinical_patient.tsv")
to_tsv(samp, "brca_clinical_sample.tsv")
