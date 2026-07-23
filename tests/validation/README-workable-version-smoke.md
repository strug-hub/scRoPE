# Unified-route workable-version smoke test

This deliberately small end-to-end test exercises the installed package on the
bundled sample data, whose average cells per donor is realistic for a software
plumbing check. It fits two genes through legacy LN, joint-LN, and PGMM and one
gene through the corrected first-order HL route. It checks the shared result
schema, route provenance, method-specific test availability, saved-result
round-tripping, and elapsed-time recording.

It is not a calibration study and must not be used to support type-I error,
coverage, power, FDR, or method-comparison claims.

After installing the exact release-candidate source into a clean library, pass
the installed package path explicitly. The script refuses to run against a
different installation:

```sh
R_LIBS=/tmp/scrope-rc-library \
  Rscript tests/validation/run_workable_version_smoke.R \
  /tmp/scrope-workable-smoke.rds \
  /tmp/scrope-rc-library/scRoPE
```

The script exits nonzero on a contract failure and otherwise prints
`Workable-version software smoke: PASS`.
