#!/usr/bin/env python3

import csv
import re
from pathlib import Path

import numpy as np
import uproot


LOG_DIR = Path("output")
SUMMARY_CSV = LOG_DIR / "aggregate_summary.csv"
SIGNAL_CSV = LOG_DIR / "selected_event_signals.csv"
CHARGE_PER_MEV_KE = 21.78


def read_completed_logs():
    results = {}

    for log in LOG_DIR.glob("slurm-mcp_geant_*.out"):
        values = {}
        for line in log.read_text(errors="ignore").splitlines():
            if "=" in line:
                key, value = line.split("=", 1)
                values[key.strip()] = value.strip()

        required = {
            "emitter", "product", "mass_GeV", "events",
            "input_ROOT_files", "accepted_PYTHIA_MCPs",
            "sampled_MCPs", "selected_events", "output",
        }
        if not required.issubset(values):
            continue

        match = re.search(r"_(\d+)\.out$", log.name)
        if match is None:
            continue

        output = Path(values["output"])
        if not output.is_file():
            continue

        sampled = int(values["sampled_MCPs"])
        selected = int(values["selected_events"])
        row = {
            "job_id": int(match.group(1)),
            "emitter_pdg": int(values["emitter"]),
            "product": values["product"],
            "mcp_mass_GeV": float(values["mass_GeV"]),
            "input_ROOT_files": int(values["input_ROOT_files"]),
            "accepted_PYTHIA_MCPs": int(values["accepted_PYTHIA_MCPs"]),
            "sampled_MCPs": sampled,
            "selected_events": selected,
            "geant_acceptance_fraction": selected / sampled,
            "output": str(output),
        }

        key = (row["emitter_pdg"], row["mcp_mass_GeV"])
        if key not in results or log.stat().st_mtime > results[key][0]:
            results[key] = (log.stat().st_mtime, row)

    rows = [value[1] for value in results.values()]
    return sorted(rows, key=lambda row: row["job_id"])


def cluster_labels(x, y, z):
    xyz = np.column_stack((x, y, z))
    labels = np.full(len(xyz), -1, dtype=np.int32)
    cluster_id = 0

    for seed in range(len(xyz)):
        if labels[seed] != -1:
            continue

        labels[seed] = cluster_id
        pending = [seed]

        while pending:
            current = pending.pop()
            distance2 = np.sum((xyz - xyz[current]) ** 2, axis=1)
            neighbors = np.flatnonzero((labels == -1) & (distance2 <= 1.0))
            labels[neighbors] = cluster_id
            pending.extend(neighbors.tolist())

        cluster_id += 1

    return labels


def extract_signals(summary_rows):
    signal_rows = []

    for job in summary_rows:
        with uproot.open(job["output"]) as root_file:
            hits = root_file["PixelEDep"].arrays(
                ["EventID", "EnergyDep_Primary", "x", "y", "z"],
                library="np",
            )

        for event_id in np.unique(hits["EventID"]):
            event = hits["EventID"] == event_id
            energy = hits["EnergyDep_Primary"][event]
            x = hits["x"][event]
            y = hits["y"][event]
            z = hits["z"][event]
            labels = cluster_labels(x, y, z)

            if len(np.unique(labels)) != 2:
                continue

            clusters = []
            for cluster_id in range(2):
                keep = labels == cluster_id
                clusters.append({
                    "NPixels": int(np.count_nonzero(keep)),
                    "Energy_MeV": float(energy[keep].sum()),
                    "Charge_ke": float(CHARGE_PER_MEV_KE * energy[keep].sum()),
                    "x_cm": float(x[keep].mean()),
                    "y_cm": float(y[keep].mean()),
                    "z_cm": float(z[keep].mean()),
                })

            delta = np.array([
                clusters[1]["x_cm"] - clusters[0]["x_cm"],
                clusters[1]["y_cm"] - clusters[0]["y_cm"],
                clusters[1]["z_cm"] - clusters[0]["z_cm"],
            ])
            distance = float(np.linalg.norm(delta))
            theta_z = float(np.degrees(np.arccos(abs(delta[2]) / distance)))

            signal_rows.append({
                "job_id": job["job_id"],
                "emitter_pdg": job["emitter_pdg"],
                "product": job["product"],
                "mcp_mass_GeV": job["mcp_mass_GeV"],
                "EventID": int(event_id),
                "total_primary_energy_MeV": clusters[0]["Energy_MeV"] + clusters[1]["Energy_MeV"],
                "total_primary_charge_ke": clusters[0]["Charge_ke"] + clusters[1]["Charge_ke"],
                "cluster1_NPixels": clusters[0]["NPixels"],
                "cluster1_Energy_MeV": clusters[0]["Energy_MeV"],
                "cluster1_Charge_ke": clusters[0]["Charge_ke"],
                "cluster1_x_cm": clusters[0]["x_cm"],
                "cluster1_y_cm": clusters[0]["y_cm"],
                "cluster1_z_cm": clusters[0]["z_cm"],
                "cluster2_NPixels": clusters[1]["NPixels"],
                "cluster2_Energy_MeV": clusters[1]["Energy_MeV"],
                "cluster2_Charge_ke": clusters[1]["Charge_ke"],
                "cluster2_x_cm": clusters[1]["x_cm"],
                "cluster2_y_cm": clusters[1]["y_cm"],
                "cluster2_z_cm": clusters[1]["z_cm"],
                "Distance_cm": distance,
                "Theta_z_deg": theta_z,
                "output": job["output"],
            })

    return signal_rows


def write_csv(path, rows, fieldnames):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


summary_rows = read_completed_logs()
if not summary_rows:
    raise SystemExit("No completed Geant4 results found")

write_csv(SUMMARY_CSV, summary_rows, list(summary_rows[0]))

signal_rows = extract_signals(summary_rows)
signal_fields = [
    "job_id", "emitter_pdg", "product", "mcp_mass_GeV", "EventID",
    "total_primary_energy_MeV", "total_primary_charge_ke",
    "cluster1_NPixels", "cluster1_Energy_MeV", "cluster1_Charge_ke",
    "cluster1_x_cm", "cluster1_y_cm", "cluster1_z_cm",
    "cluster2_NPixels", "cluster2_Energy_MeV", "cluster2_Charge_ke",
    "cluster2_x_cm", "cluster2_y_cm", "cluster2_z_cm",
    "Distance_cm", "Theta_z_deg", "output",
]
write_csv(SIGNAL_CSV, signal_rows, signal_fields)

print(f"Wrote {len(summary_rows)} job results to {SUMMARY_CSV}")
print(f"Wrote {len(signal_rows)} selected events to {SIGNAL_CSV}")
