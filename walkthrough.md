# Validation Walkthrough

Use this log to record validation runs, profiling results, and notable findings.

## Template
- Date:
- Script:
- Parameters:
- Runtime:
- Peak memory:
- Profiling summary:
- Findings:
- Actions taken:

---
- Date: 2026-02-05
- Script: cloud_physics_project.py
- Parameters: default constants
- Runtime: ~26s before crash
- Peak memory: not measured
- Profiling summary: cProfile aborted due to matplotlib error (plt.xlim called with 3 args). No profiling output produced.
- Findings: script crashes during plotting; see review notes.
- Actions taken: none (review only)

---
- Date: 2026-02-05
- Script: cloud_physics_project.py
- Parameters: default constants
- Runtime: ~30s
- Peak memory: not measured
- Profiling summary: not run in this pass (plotting fix validation only)
- Findings: script runs successfully; plots generated without matplotlib errors.
- Actions taken: fixed plt.xlim/plt.ylim arg usage

---
- Date: 2026-02-05
- Script: cloud_physics_project.py
- Parameters: default constants
- Runtime: ~81s
- Peak memory: not measured
- Profiling summary: not run; script crashed before completion
- Findings: radius growth change caused ascent to exceed cloud; loop overflowed array bounds (i reached 12,000,000). IndexError in descent loop.
- Actions taken: none (needs guard or larger arrays)

---
- Date: 2026-02-05
- Script: cloud_physics_project.py
- Parameters: default constants
- Runtime: ~119s
- Peak memory: not measured
- Profiling summary: not run in this pass (list-based fix validation)
- Findings: script runs to completion after switching to dynamic lists; droplet exceeds cloud but no array overflow.
- Actions taken: replaced preallocated arrays with lists; recomputed E using searchsorted; plots still generated.

---
- Date: 2026-02-05
- Script: cloud_physics_project.py
- Parameters: default constants
- Runtime: ~22s
- Peak memory: not measured
- Profiling summary: not run in this pass
- Findings: restored original E behavior; outputs match prior baseline and plots generated.
- Actions taken: reverted E update to ascent-only with original j lookup; kept dynamic lists.

---
- Date: 2026-02-05
- Script: cloud_physics_project.py
- Parameters: default constants
- Runtime: ~33s (cProfile)
- Peak memory: not measured
- Profiling summary: total 31.25s; top time in module execution; major overhead in list.append (~4.43s) and numpy.array conversions (~2.81s).
- Findings: performance dominated by looped appends and array conversion at end.
- Actions taken: removed unused pandas import; fixed duplicate plot filename to cloud-height-radius.jpg.

---
