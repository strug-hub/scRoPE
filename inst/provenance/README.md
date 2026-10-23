# Source and example-data provenance

This is a release candidate for a proposed 0.3.0 release. It does
not replace the historical source archives, assign results to figures, recover
missing runtime environments, or define the revised paper simulation/HLCA
workflows. No commit below is presented as a retrospective historical commit.

## Historical scRoPE source

The integration starts at existing commit
`30fb38c0753b3e8a876928d4f05fcb8e258e6008`. Its 121 selected package source
files match the preserved .9001 source. A present-day baseline import preserves
all 122 selected .9005 source files and leaves the parent's ancillary Git/build
files intact. The import is a 16-file delta; its per-file content identifiers
are in `source-0.2.5.9005.sha256` (relative to the package root at that import,
not the subsequently repaired tree).

Historical .9001, .9004 and .9005 source archives, together with the .9006
review candidate and their checksums, are maintained in separate research
records. They are distinct from this 0.3.0 release candidate.

The .9004 archive includes uncommitted work beyond the .9001 commit. Its exact
historical source is not represented by tagging that commit. The .9005 source
also required a snapshot import rather than a pre-existing exact Git commit.
Retained campaigns use different versions: earlier PGMM/HL sources include
.9001, selected replacement JLN sources use .9004, and newer campaigns use
.9005. Some aggregate views mix these versions. Their detailed source/result
bindings and available session records remain in the separate preserved
research records; this package contains no analysis data or private manifests.
No single package version or current local session reconstructs all campaigns.

## NEBULA and bundled examples

scRoPE incorporates code derived from NEBULA, by Liang He with contribution
from Raghav Sharma in the inspected archive:
[upstream source](https://github.com/lhe17/nebula),
[He et al. (2021)](https://doi.org/10.1038/s42003-021-02146-6).
The public repository starts with commit
`a0d38e45e5096187dbc6a7ef4a2bdfd9903a87fc` and does not contain the earlier
ancestry. Source provenance checks identify upstream NEBULA commit
`e17fb830b6fa0f692bb69e656b6a512c18212e7f`, "Version 1.5.3", dated
2024-02-15, as the development ancestor. Detailed development records are
maintained separately; this is not a Git parent relationship in the public
repository.

A separately preserved NEBULA 1.5.3 archive provides a source comparison.
Its ten R files match that upstream tree after line-ending normalization;
both example datasets and four shared native source files match exactly.
The source archive and the Git tree are distinct artifacts, not identical
byte streams. No old repository contents or history were merged into the
public package merely to establish this provenance.

Before the synthetic-example replacement in this local 0.3.0 preparation, the
two inherited example files were byte-identical to scRoPE's initial import and
the comparison archive. These are HISTORICAL hashes, not the maintained data:

| File | SHA-256 |
| --- | --- |
| data/sample_data.rda | `e7ff93ea0010078cbbcd374e8630b2600ad8cb2f0af6044c3221418e3208cb82` |
| data/sample_seurat.rda | `0ab12b1fdd041b775e2ee5aebfa92fed80335ed0ac0918ef42d3b0626c22d559` |

The historical `sample_data` contains 10 genes and 6,176 cells from 30 subject identifiers.
The inspected upstream documentation does not record its generation procedure
or study source; it must not be described as synthetic without further evidence.
The historical `sample_seurat` contains 1,000 genes and 1,000 cells. NEBULA attributes this
example to Satija Lab's *panc8.SeuratData: Eight Pancreas Datasets Across Five
Technologies* (2019; cited package version 3.0.2), distributed through
[SeuratData](https://github.com/satijalab/seurat-data). The exact subset-selection
procedure and dataset-specific redistribution record have not been recovered.


The maintained package now contains original, reproducibly generated teaching
examples: eight genes and 1,200 cells across 30 synthetic donors. Both objects
come from the same reference count matrix. The frozen parameters and pure
generator are in `inst/examples`; source regeneration and plotting scripts are
in `data-raw`. They use seed 20261020, with no seed search or gene selection.
Neither object contains the inherited NEBULA/pancreas data. The exact old files
remain available in historical source commits and preservation archives; no
historical result or source hash is relabeled. This replacement resolves the maintained example-data generation gap.
The separate 24-gene `toy_data` quick start is documented in
[its design notes](../examples/toy24/README.md); neither teaching design
contains historical biological observations.

## Historical licensing records

Upstream commit `e17fb830` declares GPL-3 in DESCRIPTION but contains GPL
Version 2 text in LICENSE. That historical LICENSE was inherited unchanged.
During earlier scRoPE development, DESCRIPTION changed from GPL-3 to GPL-2
without changing LICENSE.
The separate NEBULA 1.5.3 archive also declares GPL-3. These are historical
records; the maintained-source change below does not alter their declarations,
relabel old archives, or resolve which historical terms applied to old copies.

## Maintained-source MIT declaration

The maintained 0.3.0 candidate now declares `MIT + file LICENSE`. Lehang
Zhong is the existing author/maintainer and authorized this license for the
maintained scRoPE contributions. The new 2026 notice is for those contributions;
it does not invent historical ownership dates or an institutional rightsholder.
Authors@R remains unchanged. The complete source notice also retains exactly
`Copyright (c) 2026 nebula authors` for inherited material.

The technical basis is the inventoried inherited code's correspondence to:

- [NEBULA 1.5.8, 70593f9](https://github.com/lhe17/nebula/tree/70593f95342e623f721ec5919adc53868e25a2b0),
  whose DESCRIPTION and full notice declare MIT.
- [First MIT revision, e0a0da9](https://github.com/lhe17/nebula/commit/e0a0da9aaa8e1c4dbc8519de3efdf18c5d4f34a2),
  NEBULA 1.5.7, used to trace retained earlier variants.

The bounded reconciliation covered the 14 inherited R/native/text-build
paths and 53 declarations, including retained code behind retired dispatch.
The four differing bodies (`nebula`, `scToNeb`, `check_conv`, `opt_pml_nbm`)
were mapped at expression/block level. Retained preprocessing, convergence,
parallel and native extensions were traced separately from scRoPE changes.
No substantial unexplained old-upstream-only calculation was identified in
that inventory. Code was not replaced or numerically changed to obtain a
license-text match. A compact [source map](mit-source-map.csv) distinguishes
inherited expressions and local adaptations.

Upstream [LICENSE and LICENSE.md](../notices/nebula/README.md) are shipped
verbatim and remain in the installed package. The top-level MIT template
names the known scRoPE contributor, while the complete notice preserves the
upstream collective attribution without guessing individual rights holders.
The earlier GPL text and source declarations remain in historical commits and
preservation archives. MIT is applied to this maintained source; historical
simulation sources keep their original evidence.

This is a documented source-provenance and authorized maintained-source
change, not a comprehensive ownership opinion or clearance for every binary
distribution. No conflicting local holder or distinct per-file notice was
identified in the inspected R/src inventory. Git authorship and matching code
are evidence of lineage, not proof of every ownership interest.

### Dependencies and generated interfaces

[Dependency notices](../notices/DEPENDENCIES.md) are separate. scRoPE does
not relicense R, Rcpp, RcppEigen/Eigen or other dependencies. Their source,
headers, libraries and notices retain their own terms. The Rcpp-generated
wrapper markers and dependency references remain unchanged. No dependency
source/header subtree is vendored by this change. Binary redistribution and
bundling must account for the terms of the actual dependencies included.
