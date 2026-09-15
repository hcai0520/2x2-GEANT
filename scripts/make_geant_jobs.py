#!/usr/bin/env python3

import csv
import glob
import sys
from pathlib import Path


EMITTERS = {
    111: "pi0",
    221: "eta",
    331: "etap",
    113: "rho0",
    223: "omega",
    333: "phi",
    443: "jpsi",
    0: "dy",
}


mass_file = Path(sys.argv[1])
requested_emitters = [int(value) for value in sys.argv[2:]]

for emitter_pdg in requested_emitters:
    if emitter_pdg not in EMITTERS:
        raise SystemExit(f"Unknown emitter: {emitter_pdg}")

masses = []
for line in mass_file.read_text().splitlines():
    line = line.strip()
    if line and not line.startswith("#"):
        float(line)
        masses.append(line)

rows = []
for emitter_pdg in requested_emitters:
    emitter = EMITTERS[emitter_pdg]
    for mass_text in masses:
        mass_tag = f"{float(mass_text):.6f}".replace(".", "p")
        if emitter_pdg == 0:
            pattern = (
                "/pscratch/sd/h/hogcai/mcp-dy-production/raw/"
                f"emitter_dy/mass_{mass_tag}/*.root"
            )
        else:
            pattern = (
                "/pscratch/sd/h/hogcai/mcp-meson-production-50pt/raw/"
                f"emitter_{emitter}/mass_{mass_tag}/*.root"
            )

        if glob.glob(pattern):
            rows.append({
                "job_id": len(rows),
                "emitter": emitter_pdg,
                "mcp_mass_GeV": mass_text,
            })

with open("job.csv", "w", newline="") as output:
    writer = csv.DictWriter(
        output, fieldnames=["job_id", "emitter", "mcp_mass_GeV"]
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} jobs to job.csv")
