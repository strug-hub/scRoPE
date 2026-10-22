# Dependency and generated-interface notices

The MIT declaration applies to maintained scRoPE source contributions and the
inventoried NEBULA material covered by its separately retained MIT notices.
It does not relicense dependencies, generated-interface tooling, or their
compiled components. No dependency source or header subtree was added by
this change. Rcpp-generated wrapper markers remain intact.

The following are the installed dependency metadata observed on 22 October
2026, not grants from scRoPE or promises about future versions:

| Dependency | Observed version | Its declared license |
|:--|:--|:--|
| Rcpp | 1.1.2 | GPL (>= 2) |
| RcppEigen | 0.3.4.0.2 | GPL (>= 2) \| file LICENSE |
| Rfast | 2.1.5.2 | GPL (>= 2.0) |
| nloptr | 2.2.1 | LGPL (>= 3) |
| trust | 0.1-9 | MIT + file LICENSE |
| Matrix | 1.7-4 | GPL (>= 2) \| file LICENCE |
| parallelly | 1.48.0 | LGPL (>= 2.1) |
| doFuture | 1.2.2 | LGPL (>= 2.1) |
| future | 1.70.0 | LGPL (>= 2.1) |
| foreach | 1.5.2 | Apache License (== 2.0) |
| doRNG | 1.8.6.3 | GPL (>= 2) |
| Seurat | 5.3.0 | MIT + file LICENSE |
| SingleCellExperiment | 1.30.1 | GPL-3 |

RcppEigen's own LICENSE further distinguishes Eigen's MPL-2 material from
the GPL-2-or-later R integration; consult the exact installed headers and
notices. R and Rmath, system numerical libraries, optional plotting tools and
transitive dependencies retain their respective terms. This table is not a
complete transitive-dependency or binary-distribution analysis.

If distributing compiled artifacts or bundling dependencies, retain their
applicable notices and satisfy the terms of the actual components included.
This local source declaration is not an assertion that an entire compiled
installation or every combined binary is MIT-only. Dependency licenses are
not copied here and relabeled as scRoPE's work.
