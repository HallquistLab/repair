# Optimizer-default clarification

Recorded during the frozen robust rerun, before full-grid performance inspection.
The numerical fitting strategy and simulation design are unchanged.

The development CSV label `nlopt_tight` and the initial protocol's phrase
"tighter nloptwrap" refer to explicitly setting xtol_abs=ftol_abs=1e-8 and raising
maxeval to 200000. Inspection of installed lme4 2.0.6 shows that 1e-8 is already
its default stopping tolerance (default maxeval=100000). Thus this comparison
increased the evaluation budget; it did not tighten the installed default
stopping tolerances. It left all 120 selected failures unresolved. It does not
show what an even tighter nloptwrap tolerance would do.

BOBYQA here is minqa's BOBYQA, with rhobeg=.002, rhoend=1e-8 and maxfun=200000.
The default nloptwrap algorithm is NLOPT_LN_BOBYQA; switching implementations,
initialization and numerical controls is the meaningful distinction. The robust
recipe remains nloptwrap initialization followed by likelihood-checked BOBYQA
with bounded deterministic retries. No convergence check is relaxed.

The initial design/source snapshots are intentionally preserved byte-for-byte
for the active run's hash verification. This note corrects their descriptive
terminology; the executed controls and diagnostic CSVs are explicit and unchanged.
