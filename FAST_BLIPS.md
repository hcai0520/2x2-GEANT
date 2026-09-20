# Weighted two-blip mode

This mode transports the same input MCP through the existing geometry and physics,
then generates one weighted two-point response on that MCP's transported path.
It bypasses pixel bookkeeping and writes `mcp_fast.root`. It does not force two
microscopic energy deposits. No build, simulation, or notebook execution was
performed while making this change.

## Configuration

Keep the existing `MCP_INPUT_CSV`, `MCP_MASS_GEV`, and run macro. Set
`MCP_FAST_BLIPS=1` before launching the rebuilt executable. Unset it, or set it to
`0`, for the original `PixelEDep` mode and `mcp_02.root` output. Configuration is
read once per process. A rebuild is required because C++ sources were added.
The existing CMake source glob includes the new file after reconfiguration.

The updated `scripts/submit_geant_array.slurm` defaults to `MCP_FAST_BLIPS=1`.
It skips the pixel filter and saves `mcp_fast_<product>_<mass>.root` in
`/pscratch/sd/h/hogcai/mcp-geant4/output-fast` unless `OUTPUT_DIR` is supplied.
Set `MCP_FAST_BLIPS=0` explicitly to use the original batch workflow. The job
table, sample count, and Slurm resources are unchanged. Relative profile paths
are resolved from the project directory before changing to the job directory.

`filter.py` remains available for original pixel output. Do not apply it to
`FastBlips`; thresholding and the definition of a reconstructed blip are already
included in the calibration. Geant4 production cuts are unchanged: those are
transport settings, not the analysis threshold.

## Default calibration

The default uses the saved `blip_counts.ipynb` result:

```
N0 = 49,894,650
N1 = 105,251
reference active length = 83.28 cm
mu_reference = N1 / N0 = 0.0021094646
lambda = mu_reference / 83.28 cm
```

This is a constant-intensity approximation using the measured mean, not the
bin-by-bin measured spatial curve. For each path of active length L,
`mu(L) = lambda * L`. Thus the integrated mean increases with active length.
The full-length two-blip probability is approximately `2.22e-6`, before any
distance or angular selection. There is no extra threshold efficiency factor.

To use the measured cumulative curve instead, set `MCP_BLIP_PROFILE` to a CSV
with two numeric columns, no header, and optional `#` comments:

```
# cumulative_active_length_cm, probability_at_least_one_blip
0,0
# Add the measured points here, in increasing length order.
```

At least two points are required. The first must be `(0,0)`; lengths must strictly
increase and probabilities must be nondecreasing in `[0,1)`. The code converts
`P_at_least_one(L)` to `mu(L) = -log(1 - P_at_least_one(L))`, then interpolates
**mu** linearly between points. It samples each point by inverting
`mu(s) = U * mu(L)`. It does not use the first-blip CDF directly as the CDF for
each member of the pair. Flat portions have zero sampling probability.

The profile is indexed by cumulative distance along the active path from its
first fiducial entry, not by physical z or by normalized fractional depth. It
must cover the entire active length of every simulated track. Longer paths
cause an explicit error; the code neither stretches nor extrapolates the table.
Angled tracks can exceed the 83.28 cm axial reference. No measured numeric
profile was exported or supplied with this change.

## Path and sampling

Only primary MCP steps whose pre-step volume is `Prisms_M0` through `Prisms_M3`
contribute. Each step chord is clipped to the existing fiducial bounds, in cm:

```
x: [-63.931, 63.931]
y: [-51.85, 51.85]
z: [-54.32, -12.68] or [12.68, 54.32]
```

World, module gaps, and material outside these bounds add no active length.
The underlying transport still crosses those regions normally. The cumulative
coordinate starts at zero and does not advance through gaps. Geometric chord
lengths are used consistently for sampling and position interpolation, rather
than mixing Geant4's corrected step length with geometric distances.

Two independent active coordinates are sampled conditional on exactly two
points, ordered along the path, and mapped back to global positions on the
same piecewise path. Each candidate has

```
PairWeight = exp(-mu(L)) * mu(L)^2 / 2
```

There is no rejection based on that small probability: it is an output weight.
There is also no smearing, distance cut, or resampling of nearby points. The
20 cm cut belongs in downstream analysis. Applying DBSCAN to these two points
would introduce a different model; no new cluster-merging efficiency is applied.

## Output and normalization

`FastBlips` has one row per event reaching the end-event action, including events
without a usable active path. Its fields are:

| Fields | Meaning |
|---|---|
| `RunID`, `EventID`, `TrackID` | Event and primary identity; one MCP per event |
| `Status` | 1: candidate; 0: zero active length or zero integrated rate; 2: aborted |
| `ActiveSegments` | Number of retained fiducial step chords |
| `ActiveLength_cm`, `MeanBlips`, `PairWeight` | Active length, integrated rate, exactly-two probability |
| `s1_cm`, `s2_cm` | Ordered cumulative active coordinates |
| `x1_cm`, `y1_cm`, `z1_cm`, `x2_cm`, `y2_cm`, `z2_cm` | Global point positions |
| `Distance_cm` | Euclidean separation, without a cut |

Unavailable positions and distances are NaN; unusable and aborted events have
zero pair weight. An aborted event's length is partial and its mean is not
evaluated. `RunSummary` records the run ID and Geant4's event count. Compare this
with the number of event rows and check for aborted events before normalization;
an incomplete run must not silently be treated as a complete input sample.
Subsequent runs in the same process use `mcp_fast_run<ID>.root`.

For a complete run, divide the sum of passing `PairWeight` values by **all
incident MCP events**, including zero-path events. Do not divide by just the
candidate count. In physical event weights, multiply the original production
and geometry normalization per simulated incident MCP by `PairWeight` once.
The existing pixel aggregator and signal notebook do not read this new schema;
they were not changed.

## Model limits and review status

This implements the independent Poisson response approximation. The calibration
was measured at a fixed incident energy; applying it to other energies and masses
assumes the same effective rate. The code does not infer an unmeasured energy or
charge dependence. A supplied length profile carries the same transfer assumption.

The points lie on Geant4 step chords. They are not reconstructed centroids and
do not include the centroid displacement caused by a secondary deposit. Nor does
an unconditional transported path reproduce the path distribution conditioned
on two visible interactions. Centroid/angular agreement, clustering effects,
step-size dependence, and speedup require validation against original-mode
results. The new mode still transports secondary particles. Sampling also uses
the Geant4 random stream, so subsequent event histories need not match original
mode event by event even with the same seed.

Changes were reviewed statically only. Compilation and runtime behavior have
not been tested.
