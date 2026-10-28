#!/usr/bin/env python3
"""
Facts about a 2-D run, read from its log.

Every figure tool used to hard-code the storm ("100-yr 6-hr event"), the run
length ("9-hour transient run") and, in plot_profiles.py, the whole inflow
table. Run on a 24-hour case those captions are simply false, and the
NM15-141 figures shipped in Report_for_C10 carried a 6-hr caption and a
pre-culvert-10 inflow table on 24-hour post-culvert-10 results. Reading the
facts from the log means a figure cannot disagree with the run that made it.

    from run_facts import facts
    f = facts("out_nm141_1m_6hr_carved_2m_postC10")
    f["storm"]      -> "100-yr 6-hr"
    f["sim_hours"]  -> 9.0
    f["inflows"]    -> [("Reach_147  south", 5414), ("Reach_122  west", 4664)]

The log lives either INSIDE the run folder (run.log, the post-C10 layout) or
beside it (<run>.log, the original layout); both are found.
"""
import os, re

_EDGE = re.compile(r"^inflow\s+(\S+?),\s*(.*?)\s+(SOUTH|WEST|EAST|NORTH)\s+edge.*?\(\s*([\d,]+)\s*cfs\)", re.M)


def find_log(run):
    for p in (os.path.join(run, "run.log"), run.rstrip("/") + ".log"):
        if os.path.exists(p):
            return p
    return None


def facts(run):
    f = dict(storm=None, sim_hours=None, inflows=[], inflow_cfs=None,
             inflow_acft=None, rain_acft=None, rain_in=None, total_acft=None,
             hms_cfs=None, hms_acft=None, terrain=None, log=None)
    p = find_log(run)
    if p is None:
        return f
    f["log"] = p
    t = open(p, errors="ignore").read()

    m = re.search(r"^storm\s+(\d+)-hour\s+100-yr", t, re.M)
    if m:
        f["storm"] = "100-yr %s-hr" % m.group(1)
    m = re.search(r"^duration\s+([\d.]+)\s*s", t, re.M)
    if m:
        f["sim_hours"] = float(m.group(1)) / 3600.0
    # A log with no storm line predates that line; the simulated length is 1.5x
    # the storm in every case built here (9 h for 6-hr, 36 h for 24-hr).
    if f["storm"] is None and f["sim_hours"]:
        f["storm"] = "100-yr %d-hr" % (6 if f["sim_hours"] < 20 else 24)

    for name, desc, edge, cfs in _EDGE.findall(t):
        short = re.sub(r"\s+", " ", desc).strip()
        short = re.sub(r"^(main wash|west|east|north|south)\b.*?\((\w+)\)", r"\1", short)
        f["inflows"].append(("%s  %s" % (name, edge.lower()), int(cfs.replace(",", ""))))
    m = re.search(r"---\s*\d+ inflows, peaks sum to\s+([\d,]+)\s*cfs, volume\s+([\d,]+)\s*ac-ft", t)
    if m:
        f["inflow_cfs"] = int(m.group(1).replace(",", ""))
        f["inflow_acft"] = float(m.group(2).replace(",", ""))
    m = re.search(r"HMS Outlet:\s*([\d,]+)\s*cfs,\s*([\d,]+)\s*ac-ft", t)
    if m:
        f["hms_cfs"] = int(m.group(1).replace(",", ""))
        f["hms_acft"] = float(m.group(2).replace(",", ""))
    m = re.search(r"^\s+([\d.]+)\s*in\s*\(([\d.]+)\s*ac-ft\)", t, re.M)
    if m:
        f["rain_in"] = float(m.group(1)); f["rain_acft"] = float(m.group(2))
    if f["inflow_acft"] is not None and f["rain_acft"] is not None:
        f["total_acft"] = f["inflow_acft"] + f["rain_acft"]
    m = re.search(r"^terrain\s+(\S+)", t, re.M)
    if m:
        f["terrain"] = m.group(1)
    return f


def storm_caption(run, default="100-yr 6-hr"):
    return (facts(run).get("storm") or default) + " event"


if __name__ == "__main__":
    import sys, json
    for r in sys.argv[1:]:
        print(r); print(json.dumps(facts(r), indent=2)); print()
