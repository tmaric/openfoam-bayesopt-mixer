# Source geometry evidence and scientific framing

Status on 2026-07-28: **closed as an M10-inspired reconstruction**. Published
sources support the two-layer serpentine-crossing topology, nominal dimensions,
and the locations of inter-layer communication. They do not uniquely define the
fluid solid used in the original CFD. This repository therefore makes no exact
M10 reproduction claim.

## Evidence reviewed

The starting review is Raza, Hossain, and Kim, *Micromachines* 11 (2020) 455,
DOI [`10.3390/mi11050455`](https://doi.org/10.3390/mi11050455). Figure 10 and
Section 2.5 establish that M-10 has bottom N-shaped and top inverse-N-shaped
segments and that the streams communicate at vertical sections and crossing
intersections. The figure gives `H=1070 um`, `P=640 um`, `b=150 um`, `d=150
um`, and `w=300 um`.

The primary device paper is Hossain et al., *Chemical Engineering Journal* 327
(2017) 268-277, DOI
[`10.1016/j.cej.2017.06.106`](https://doi.org/10.1016/j.cej.2017.06.106).
Its model description confirms that the channels interconnect through the
middles of X shapes and the vertical sections. It reports two separately
fabricated PDMS microchannel layers and a 30 um nominal CFD element size, but
does not publish an exchangeable CAD model or fabrication mask.

The primary follow-up optimization by Hossain and Kim, *Chemical Engineering &
Technology* 40 (2017) 2212-2220, DOI
[`10.1002/ceat.201700437`](https://doi.org/10.1002/ceat.201700437), reuses the
device and confirms the same five nominal dimensions and inlet/outlet areas.
It uses only `w/P`, `H/P`, and `d/P` as design variables. Its flow description
reports low-Re layer transfer at crossing nodes and additional transfer through
vertical sections at high Re, which supports both connection families but does
not dimension their open masks.

Repository and web searches found no author-provided CAD, numerical-domain
mesh, fabrication mask, supplementary aperture drawing, or explicit aperture
area. The published schematics also do not uniquely specify inlet/outlet lead
transitions or whether internal unions were sharp, rounded, or fabrication-
smoothed.

## Inference and decision

Opening the complete projected overlap of face-adjacent channel layers is a
reasonable interpretation of a bonded two-layer PDMS device, and it is retained
as the baseline. It is still an inference rather than a source-defined
dimension. The present CAD additionally chooses short leads and sharp
rectangular unions.

Consequently:

- literature values are external mechanism and magnitude benchmarks, not an
  acceptance target that may be met by tuning unspecified geometry;
- fine-grid disagreement with the review is reported, not calibrated away;
- all generated manifests and campaign metadata identify the geometry as
  `m10_inspired_reconstruction` with `exact_reproduction_claim: false`;
- BO claims apply only to this parameterized M10-inspired family;
- comparison with the published M10 remains qualified by geometry and solver
  differences.

The existing folder name remains appropriate because it describes the actual
topology. Renaming it to imply exact M10 identity would make the scientific
claim less accurate.
