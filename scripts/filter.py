#!/usr/bin/env python3

import sys
import time
from pathlib import Path

import numpy as np
import uproot


INPUT = Path(sys.argv[1])
OUTPUT = Path(sys.argv[2])
THRESHOLD_DIR = Path("/global/homes/h/hogcai/data")

BRANCHES = [
    "EventID", "PixelID", "VolName",
    "EnergyDep", "EnergyDep_Primary",
    "EnergyDep_Primary_only", "EnergyDep_Secondary_only",
    "x", "y", "z",
    "x_start", "y_start", "z_start",
    "x_end", "y_end", "z_end",
    "IOGroup", "Plane", "PixelY", "PixelZ", "ThresholdLocalKey",
]

INT_BRANCHES = {
    "EventID", "PixelID", "VolName", "IOGroup", "Plane",
    "PixelY", "PixelZ", "ThresholdLocalKey",
}

GRIDS = [(140, 280), (140, 280), (160, 320), (140, 280)]
TABLES = []

for module, (nx, ny) in enumerate(GRIDS):
    table = np.load(THRESHOLD_DIR / f"thresholds_module{module}.npz")
    keys = table["keys"].astype(np.int64)
    order = np.argsort(keys)
    TABLES.append((
        keys[order],
        table["values"][order].astype(float) / 1000.0,
        float(np.asarray(table["default"]).flat[0]) / 1000.0,
        nx,
        ny,
        int(keys.min()),
    ))


def threshold_ke(io_group, plane, pixel_y, pixel_z):
    module = (io_group.astype(int) - 1) // 2
    result = np.full(module.shape, np.inf)

    for m in range(4):
        mask = module == m
        if not np.any(mask):
            continue

        keys, values, default, nx, ny, offset = TABLES[m]
        key = offset + pixel_z[mask] + nx * (pixel_y[mask] + ny * plane[mask])
        position = np.searchsorted(keys, key)
        found = (
            (position < len(keys))
            & (keys[np.minimum(position, len(keys) - 1)] == key)
        )
        result[mask] = default
        destination = np.flatnonzero(mask)[found]
        result.flat[destination] = values[position[found]]

    return result


def concatenate(first, second):
    if first is None:
        return second
    return {
        name: np.concatenate((first[name], second[name]))
        for name in second
    }


def count_clusters(x, y, z):
    xyz = np.column_stack((x, y, z))
    labels = np.full(len(xyz), -1, dtype=np.int32)
    n_clusters = 0

    for seed in range(len(xyz)):
        if labels[seed] != -1:
            continue

        if n_clusters == 2:
            return 3

        labels[seed] = n_clusters
        pending = [seed]

        while pending:
            current = pending.pop()
            distance2 = np.sum((xyz - xyz[current])**2, axis=1)
            neighbors = np.flatnonzero((labels == -1) & (distance2 <= 1.0))
            labels[neighbors] = n_clusters
            pending.extend(neighbors.tolist())

        n_clusters += 1

    return n_clusters


def select_entries(data):
    if len(data["EventID"]) == 0:
        return data["entry"], 0

    event_ids, starts = np.unique(data["EventID"], return_index=True)
    stops = np.r_[starts[1:], len(data["EventID"])]
    keep = np.zeros(len(data["EventID"]), dtype=bool)
    selected_events = 0

    for start, stop in zip(starts, stops):
        if stop - start >= 2 and count_clusters(
            data["x"][start:stop],
            data["y"][start:stop],
            data["z"][start:stop],
        ) == 2:
            keep[start:stop] = True
            selected_events += 1

    return data["entry"][keep], selected_events


def read_rows(tree, names, entries, basket_offsets):
    """Read only ROOT baskets containing the requested, ordered row indices."""
    result = {}
    for name in names:
        branch = tree[name]
        offsets = basket_offsets[name]
        baskets = np.unique(np.searchsorted(offsets, entries, side="right") - 1)
        groups = np.split(baskets, np.flatnonzero(np.diff(baskets) > 1) + 1)
        pieces = []
        for group in groups:
            start, stop = int(offsets[group[0]]), int(offsets[group[-1] + 1])
            left, right = np.searchsorted(entries, [start, stop])
            values = branch.array(entry_start=start, entry_stop=stop, library="np")
            pieces.append(values[entries[left:right] - start])
        result[name] = np.concatenate(pieces)
    return result


def passing_rows(tree, energy, entry_start, basket_offsets, threshold_floor):
    # A row below the lowest threshold (including defaults) cannot pass anywhere.
    candidate = np.flatnonzero(
        np.isfinite(energy) & (energy > 0) & (21.78 * energy > threshold_floor)
    )
    if len(candidate) == 0:
        return None

    entries = candidate + entry_start
    channels = read_rows(
        tree, ["IOGroup", "Plane", "PixelY", "PixelZ"], entries, basket_offsets
    )
    threshold = threshold_ke(
        channels["IOGroup"], channels["Plane"], channels["PixelY"], channels["PixelZ"]
    )
    entries = entries[21.78 * energy[candidate] > threshold]
    if len(entries) == 0:
        return None

    data = read_rows(tree, ["EventID", "x", "y", "z"], entries, basket_offsets)
    fiducial = (
        (data["x"] >= -63.931) & (data["x"] <= 63.931)
        & (data["y"] >= -51.85) & (data["y"] <= 51.85)
        & (
            ((data["z"] >= -54.32) & (data["z"] <= -12.68))
            | ((data["z"] >= 12.68) & (data["z"] <= 54.32))
        )
    )
    if not np.any(fiducial):
        return None
    data["entry"] = entries
    return {name: values[fiducial] for name, values in data.items()}


OUTPUT.parent.mkdir(parents=True, exist_ok=True)
branch_types = {
    name: np.dtype("int32" if name in INT_BRANCHES else "float64")
    for name in BRANCHES
}

selected_event_count = 0
carry = None
selected_entries = []
threshold_floor = min(
    float(np.nanmin(np.r_[values, default]))
    for _, values, default, *_ in TABLES
)
started = last_progress = time.monotonic()

with uproot.open(INPUT, array_cache=None) as source, uproot.recreate(OUTPUT) as destination:
    source_tree = source["PixelEDep"]
    output_tree = destination.mktree("PixelEDep", branch_types)
    basket_offsets = {
        name: np.asarray(source_tree[name].entry_offsets, dtype=np.int64)
        for name in BRANCHES
    }
    print(f"Scanning {int(source_tree.num_entries):,} pixel rows", flush=True)

    for chunk, report in source_tree.iterate(
        ["EnergyDep_Primary"], step_size=1_000_000, library="np", report=True
    ):
        filtered = passing_rows(
            source_tree, chunk["EnergyDep_Primary"], report.tree_entry_start,
            basket_offsets, threshold_floor,
        )
        if filtered is not None:
            combined = concatenate(carry, filtered)
            # Retain the last surviving event until a later event is encountered.
            complete = combined["EventID"] != combined["EventID"][-1]
            complete_data = {name: values[complete] for name, values in combined.items()}
            carry = {name: values[~complete] for name, values in combined.items()}
            entries, count = select_entries(complete_data)
            if len(entries):
                selected_entries.append(entries)
            selected_event_count += count

        now = time.monotonic()
        if now - last_progress >= 30 or report.tree_entry_stop == source_tree.num_entries:
            print(
                f"Scanned {report.tree_entry_stop:,}/{int(source_tree.num_entries):,} rows; "
                f"selected events so far: {selected_event_count}; elapsed {now - started:.0f}s",
                flush=True,
            )
            last_progress = now

    if carry is not None:
        entries, count = select_entries(carry)
        if len(entries):
            selected_entries.append(entries)
        selected_event_count += count

    if selected_entries:
        entries = np.concatenate(selected_entries)
        print(f"Copying {len(entries):,} selected pixel rows with all output fields", flush=True)
        for start in range(0, len(entries), 100_000):
            output_tree.extend(read_rows(
                source_tree, BRANCHES, entries[start:start + 100_000], basket_offsets
            ))

print(f"selected_events={selected_event_count}")
print(f"output={OUTPUT}")
