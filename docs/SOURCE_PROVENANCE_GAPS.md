# Source provenance and rights documentation gaps (postpartum-depression child, 2026-10-08 audit record)

This file records what this repository's own documents say about source hashes and terms. It is a documentation record, not a license verdict, not an integrity certificate, and not a clearance of any recorded rights hold. A matching claim below means two recorded hash strings are equal at a named field in pinned repository documents. It does not mean any assay bytes were re-downloaded or verified. No rights record was found for these sources, so reuse and redistribution are NOT cleared by this file.

Audited child commit: `4ddc85483ab56982d020a5f54faafe2316d26080`  
Compared shared-parent repo: `mega27-25-biomarkers-underserved-diseases` at `36ff87995f4fe0e1d08f9e2b2cd885d341b88937`

Rights: No source-specific rights check found in inspected disease-child documentation. Hash matches are not legal clearance or byte verification.

## Per-child correction

Recorded as written: no matching shared-parent record exists for this child's hashes. The child describes them as transferred from an earlier PPD branch that was not recovered, and they are not present in the inspected main shared-parent snapshot. They are child-specific pins, not matches.

## 1. Result-file hash records in this child (0)

Hash values are shown as 12-hex display prefixes only; the full value is in the cited file and field. Classes are kept separate on purpose: an expression/assay source hash, a GEO series-matrix hash (metadata that may include expression), a reference annotation, and a published artifact are different kinds of record.


## 3. Child-specific pins with no shared-parent match

- `projects/postpartum_depression/sources/GSE44132_analysis_summary.json` line(s) [4] hash `b2d715ec45c2...` (GSE44132); class: methylation_source_matrix; Not present in inspected main shared-parent snapshot; child describes transfer from earlier PPD branch, branch not recovered
- `projects/postpartum_depression/sources/GSE335141_analysis_summary.json` line(s) [4] hash `13002634be22...` (GSE335141); class: methylation_source_matrix; Not present in inspected main shared-parent snapshot; child describes transfer from earlier PPD branch, branch not recovered
- `projects/postpartum_depression/sensitivity/rebuild_gse45603.py` line(s) [22] hash `da120f052d21...` (); class: raw_series_matrix; No shared main-parent counterpart path in snapshot
- `projects/postpartum_depression/sensitivity/rebuild_gse45603.py` line(s) [23] hash `c914fdbe1130...` (); class: platform_reference_annotation; No shared main-parent counterpart path in snapshot

## 4. Metadata crosswalk tables (not assay bytes)

These per-sample tables carry hash columns describing GEO source-response metadata. They are not blanket expression-matrix integrity.

- `projects/postpartum_depression/sources/methylation_used_record_manifest.csv`: columns ['sha256'], 139 records; file identical to a shared-parent file: False; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/postpartum_depression/sources/GSE44132_source_crosswalk.csv`: columns ['source_sha256', 'series_sha256'], 55 records; file identical to a shared-parent file: False; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/postpartum_depression/sources/GSE45603_used_sample_crosswalk.csv`: columns ['source_sha256'], 48 records; file identical to a shared-parent file: True; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/postpartum_depression/sources/GSE335141_sample_crosswalk.csv`: columns ['source_sha256'], 82 records; file identical to a shared-parent file: True; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/postpartum_depression/sources/GSE290313_used_sample_crosswalk.csv`: columns ['source_sha256'], 119 records; file identical to a shared-parent file: True; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/postpartum_depression/sources/GSE335141_source_crosswalk.csv`: columns ['source_sha256', 'series_sha256'], 82 records; file identical to a shared-parent file: False; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity

## Coverage boundary

- Source accessions listed in the child manifest: 2.
- Of those, 2 have no mapped result-file payload hash in this audit: GSE290313, GSE45603.
- 0 hash fields inherited from other diseases' records are not counted toward this child.
- Records absent from child manifest can still have result-file hashes. Neither presence nor absence proves assay acquisition by this scout. Other-disease copied hashes never count toward this child.
