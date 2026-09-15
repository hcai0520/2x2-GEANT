#!/usr/bin/env python3

import glob
import os
import sys
from pathlib import Path

import numpy as np
import uproot


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

N_SAMPLE = int(os.environ.get("N_SAMPLE", 10_000_000))
BINS_SAMPLE = 50
SAMPLE_SEED = 1331
L_DETECTOR_M = 1040.0


mass = float(sys.argv[1])
emitter_pdg = int(sys.argv[2])
output = Path(sys.argv[3])
mass_tag = f"{mass:.6f}".replace(".", "p")
emitter = EMITTERS[emitter_pdg]

if emitter_pdg == 0:
    files = sorted(glob.glob(
        f"/pscratch/sd/h/hogcai/mcp-dy-production/raw/"
        f"emitter_dy/mass_{mass_tag}/*.root"
    ))
    accepted_branch = "accepted"
else:
    files = sorted(glob.glob(
        f"/pscratch/sd/h/hogcai/mcp-meson-production-50pt/raw/"
        f"emitter_{emitter}/mass_{mass_tag}/*.root"
    ))
    accepted_branch = "passed_geometry"

if not files:
    raise SystemExit(f"No ROOT files found for emitter={emitter}, mass={mass:g} GeV")

d = uproot.concatenate(
    [f"{path}:mcp_spectra" for path in files],
    ["mcp_pdg", "E_GeV", "theta_x_rad", "theta_y_rad", accepted_branch],
    library="np",
)

keep = (
    (d[accepted_branch] == 1)
    & (d["E_GeV"] > mass)
    & np.isfinite(d["E_GeV"])
    & np.isfinite(d["theta_x_rad"])
    & np.isfinite(d["theta_y_rad"])
    & (np.hypot(d["theta_x_rad"], d["theta_y_rad"]) > 0)
)

energy_original = d["E_GeV"][keep]
theta_x_original = d["theta_x_rad"][keep]
theta_y_original = d["theta_y_rad"][keep]
mcp_pdg_original = d["mcp_pdg"][keep]

if len(energy_original) == 0:
    raise SystemExit(f"No accepted MCPs for emitter={emitter}, mass={mass:g} GeV")

theta_original = np.arctan(np.hypot(theta_x_original, theta_y_original)) * 1e3
phi_original = np.arctan2(theta_y_original, theta_x_original)
theta_min, theta_max = theta_original.min(), theta_original.max()
u_original = np.log(energy_original - mass)
v_original = np.log(theta_original)
u_edges = np.linspace(u_original.min(), u_original.max(), BINS_SAMPLE + 1)
v_edges = np.linspace(v_original.min(), v_original.max(), BINS_SAMPLE + 1)
phi_edges = np.linspace(-np.pi, np.pi, BINS_SAMPLE + 1)

H, _ = np.histogramdd(
    (u_original, v_original, phi_original),
    bins=(u_edges, v_edges, phi_edges),
)
probability = H.ravel() / H.sum()

rng = np.random.default_rng(SAMPLE_SEED)
energy_sampled = np.empty(N_SAMPLE)
theta_x_sampled = np.empty(N_SAMPLE)
theta_y_sampled = np.empty(N_SAMPLE)
filled = 0

while filled < N_SAMPLE:
    selected = rng.choice(H.size, size=min(1_000_000, N_SAMPLE - filled), p=probability)
    iu, iv, ip = np.unravel_index(selected, H.shape)
    energy_new = mass + np.exp(rng.uniform(u_edges[iu], u_edges[iu + 1]))
    theta_new = np.exp(rng.uniform(v_edges[iv], v_edges[iv + 1]))
    phi_new = rng.uniform(phi_edges[ip], phi_edges[ip + 1])

    slope = np.tan(theta_new * 1e-3)
    tx, ty = slope * np.cos(phi_new), slope * np.sin(phi_new)
    x, y = L_DETECTOR_M * tx, L_DETECTOR_M * ty
    accept = (
        (theta_new >= theta_min) & (theta_new <= theta_max)
        & (((x >= -0.65) & (x <= -0.05)) | ((x >= 0.05) & (x <= 0.65)))
        & (y >= -0.70) & (y <= 0.70)
    )
    stop = filled + int(accept.sum())
    energy_sampled[filled:stop] = energy_new[accept]
    theta_x_sampled[filled:stop] = tx[accept]
    theta_y_sampled[filled:stop] = ty[accept]
    filled = stop
    del selected, iu, iv, ip, energy_new, theta_new, phi_new, slope, tx, ty, x, y, accept

mcp_pdg_sampled = rng.choice(mcp_pdg_original, size=N_SAMPLE)
del H, probability, d

output.parent.mkdir(parents=True, exist_ok=True)

chunk_size = 1_000_000
with output.open("w") as stream:
    stream.write(
        "event_index,mcp_pdg,px_GeV,py_GeV,pz_GeV,E_GeV,"
        "x_at_detector_m,y_at_detector_m\n"
    )
    for start in range(0, N_SAMPLE, chunk_size):
        stop = min(start + chunk_size, N_SAMPLE)
        energy = energy_sampled[start:stop]
        theta_x = theta_x_sampled[start:stop]
        theta_y = theta_y_sampled[start:stop]
        momentum = np.sqrt(energy**2 - mass**2)
        direction_norm = np.sqrt(1.0 + theta_x**2 + theta_y**2)

        chunk = np.column_stack((
            np.arange(start, stop),
            mcp_pdg_sampled[start:stop],
            momentum * theta_x / direction_norm,
            momentum * theta_y / direction_norm,
            momentum / direction_norm,
            energy,
            L_DETECTOR_M * theta_x,
            L_DETECTOR_M * theta_y,
        ))
        np.savetxt(
            stream,
            chunk,
            delimiter=",",
            fmt=["%d", "%d"] + ["%.17g"] * 6,
        )

print(f"emitter={emitter}")
print(f"mass_GeV={mass:g}")
print(f"input_ROOT_files={len(files)}")
print(f"accepted_PYTHIA_MCPs={len(energy_original)}")
print(f"sampled_MCPs={N_SAMPLE}")
print(f"output={output}")
