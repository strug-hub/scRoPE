# Joint-LN availability runtime harness

`run_joint_ln_availability_runtime.R` is a bounded, synthetic software timing
harness for the shadow-only candidate classifier and one-state-evaluation basin
probe. It exercises the fast, immediate-veto, and narrow-review paths.

It does not fit real data, estimate inferential calibration, or measure the
whole joint-LN optimizer. Total optimizer/polish evaluation counts remain
unavailable by design; the harness reports the probe's exact counts. Its base-R
memory measure is an R heap proxy, not process peak RSS.

After installing the exact package candidate, run:

```sh
R_LIBS=/tmp/scrope-rc-library \
  Rscript tests/validation/run_joint_ln_availability_runtime.R \
  /tmp/scrope-joint-ln-availability-runtime.rds \
  50
```

The repetition count is bounded to 20–200. The output RDS records the package
version, a path-free MD5 manifest of the installed package tree, the complete
resolution policy/threshold object, per-path median/p90/p99 elapsed time,
review fraction/time, probe evaluation counts, the R heap proxy, and a deep
round-trip check.
