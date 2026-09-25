# Secondary pilot runtime amendment — recorded before resuming

The original frozen run plan and training engine remain unchanged. This amendment concerns one documented interruption and time accounting; it does not change eligibility, cohorts, seeds, order, objectives, architecture, optimizer settings or the 2,000-update endpoints.

## Observation and preserved record

Eight planned runs completed. The ninth, human ATP8/base with seed 29, stopped at update 268 after the computer entered lid-triggered Modern Standby. On resumption, the existing time guard correctly preserved an incomplete checkpoint. The controller recorded 20,466.25 seconds of wall time, including standby, and stopped. This exceeded four elapsed clock hours because execution was suspended; it was not four hours of uninterrupted training.

The original execution/results and partial-run files are archived with hashes in `reports/secondary/interruption/` and `runs/secondary-interruption/`. No prior measurement, checkpoint time, source hash or protocol is rewritten. The original pretraining commit remains `70941dc89abd40803a153b9b407a8d3a8f779b1d`.

## Conservative accounting correction

The owner authorized approximately four hours of local training. Windows Kernel-Power event records 52920 and 52931 report 13,390.759541 and 2,329.413931 seconds of hardware deep-idle residency, respectively. Both entire standby intervals lie within the interrupted worker's lifetime. Microsoft describes hardware DRIPS as measured physical residency of the system-on-chip in its lowest power state ([Microsoft documentation](https://learn.microsoft.com/en-us/windows-hardware/design/device-experiences/modern-standby-sleepstudy-common-problem-examples)).

Exclude only those **15,720.173472 seconds** of documented hardware deep idle. Keep every other second of the original execution charged, including all uncertain and non-deep-idle standby time. This leaves **4,746.076528 seconds charged** (79.10 minutes) and **9,653.923472 seconds remaining** (160.90 minutes) under the original 14,400-second limit. The existing 120-second shutdown reserve remains in force. This is conservative budget accounting, not a measurement of exact active compute time.

This correction is applied once, with a receipt and immutable evidence hashes. A repeated application must fail. The unchanged runner then resumes the saved optimizer, random state and sampling position and skips already completed jobs. Future interruptions receive the original conservative accounting; no further credit is automatic. No operating-system power setting is changed.

## Interpretation and reporting

The eight first-seed results and one partial second-seed result were already observed when this amendment was written. Continuation is to finish the previously fixed matrix, not to select favorable settings or extend update budgets. The extra update-268 diagnostic arose from interruption; it is retained.

Raw worker and controller wall times remain available. The final report distinguishes them from budget-charged time. The interrupted seed-29 base run is excluded from hardware-speed comparisons because its timing includes suspension and cannot be made directly comparable by this conservative correction. Its completed model results remain eligible for the originally planned seed-paired scientific contrasts.

Prepare and verify the archive, publish this amendment before further training, then apply it once and run the original scheduler:

```powershell
.venv\Scripts\python.exe scripts/continue_secondary_after_standby.py prepare
.venv\Scripts\python.exe scripts/continue_secondary_after_standby.py apply
.venv\Scripts\python.exe scripts/run_secondary.py
```

The machine-readable amendment records the exact evidence, archive and script hashes. This document supplements the historical protocol; it does not pretend this timing change was registered before the first eight outcomes.
