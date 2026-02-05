# Profiling and Optimization Checklist

Use this checklist when validating performance of simulation code.

## Setup
- Record CPU model and Python version.
- Use consistent input parameters between runs.

## Profiling
- Profile with `cProfile` or `py-spy` to find hot loops.
- Measure total runtime and memory peak.

## Roofline / Efficiency Notes
- This is a Python numerical loop; expect it to be memory and interpreter bound.
- If performance is a concern, consider vectorization or numba.

## Reporting
- Note top 5 functions by time and any obvious bottlenecks.
- Capture before/after timings when optimizing.

## Latest Baseline (2026-02-05)
- Total runtime: ~31.25s (cProfile)
- Top costs: list appends (~4.43s), numpy array conversions (~2.81s)
- Notes: performance dominated by Python-level loops and append-heavy inner steps.
