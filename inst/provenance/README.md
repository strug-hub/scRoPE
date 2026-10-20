# Source and example-data provenance

This is an unpublished local preparation for a proposed 0.3.0 release. It does
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

Historical source archives have these SHA-256 identifiers:

| Version | SHA-256 |
| --- | --- |
| 0.2.5.9001 | `5cf0a0d59b383d35d658d0a81bab51b8c8df88bdd9b35b632da9f819fe3c6377` |
| 0.2.5.9004 | `191ae565d59412846c22517882c7257f3fb780e35badc5e7b84d8310e8013f50` |
| 0.2.5.9005 | `763280750f8c36664f96acdbced4d371fcd6a6fbf6e54993cf6ce3954fa0f8a8` |
| 0.2.5.9006 (unpublished review candidate) | `8c60f73bbc285f58fde406e411b72d6718285b9ea9ec36ddf458ee46629a0571` |

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
ancestry. A recovered local development repository preserves that history:

- Upstream NEBULA commit `e17fb830b6fa0f692bb69e656b6a512c18212e7f`,
  "Version 1.5.3", dated 2024-02-15, is the direct parent of the first local
  customization, `7c6ede66b0ceca839269278a1b5637ddfe674fce` (2024-10-29).
- Package renaming is recorded at
  `48f8553348ef2a831e1d3451b6a777949ded2faf` and
  `0ecc0a47b12077097e117028a3ec93ad94d10614` (2024-11-07), followed by
  the scrope rename at `32c8ec11791917a28c42efd745c5735cc2e8eda6`
  (2025-09-17).
- The recovered development HEAD,
  `70dbf56a3337fbe05d5a7cc1f161ab10070d9906`, contains all 49 tracked
  paths selected for the public root import, from 121 historical tracked paths.
  Forty-eight of those are byte-identical; README
  differs by the omission of one local working-note line. This is a verified
  content bridge, not a Git parent relationship between the two repositories.

The recovered ancestry establishes NEBULA 1.5.3 as the development ancestor.
A separately preserved NEBULA 1.5.3 comparison archive has SHA-256
`764e6a1f775ccc2192990807a0ea8a75db3ebcf8ff5d807375523fa1096785aa`.
Its ten R files match that upstream tree after line-ending normalization;
both example datasets and four shared native source files match exactly.
The source archive and the Git tree are distinct artifacts, not identical
byte streams. No old repository contents or history were merged into the
public package merely to establish this provenance.

The two example files have remained byte-identical since scRoPE's initial
import and match the comparison archive:

| File | SHA-256 |
| --- | --- |
| data/sample_data.rda | `e7ff93ea0010078cbbcd374e8630b2600ad8cb2f0af6044c3221418e3208cb82` |
| data/sample_seurat.rda | `0ab12b1fdd041b775e2ee5aebfa92fed80335ed0ac0918ef42d3b0626c22d559` |

`sample_data` contains 10 genes and 6,176 cells from 30 subject identifiers.
The inspected upstream documentation does not record its generation procedure
or study source; it must not be described as synthetic without further evidence.
`sample_seurat` contains 1,000 genes and 1,000 cells. NEBULA attributes this
example to Satija Lab's *panc8.SeuratData: Eight Pancreas Datasets Across Five
Technologies* (2019; cited package version 3.0.2), distributed through
[SeuratData](https://github.com/satijalab/seurat-data). The exact subset-selection
procedure and dataset-specific redistribution record have not been recovered.

Licensing records remain inconsistent. Upstream commit `e17fb830` declares
GPL-3 in DESCRIPTION but contains GPL Version 2 text in LICENSE. That LICENSE
was inherited unchanged. The later scRoPE commit
`40d19b762c4d1a14f90cf730bb6d2a54bc2aef56` (2026-01-13) changed
DESCRIPTION from GPL-3 to GPL-2 and updated README; it did not change LICENSE.
The separate NEBULA 1.5.3 archive also declares GPL-3. This records where the
metadata changed without deciding which license terms apply or establishing a
separate data rights grant. Current license fields/files are left unchanged.
Clarify the conflicting upstream and later records before publication; local
test success does not resolve licensing provenance.
