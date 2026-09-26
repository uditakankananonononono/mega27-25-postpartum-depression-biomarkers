# 2026 PPD PBMC paper raw-data access check

Checked September 27, 2026, 04:50 IST. The [September 16 Dove Medical Press paper](https://www.dovepress.com/peripheral-blood-transcriptomic-features-distinguishing-postpartum-dep-peer-reviewed-fulltext-article-IJWH) studies 11 PPD, 16 healthy postpartum controls and 20 MDD patients using **PBMC RNA-seq**. It reports exploratory unadjusted P <= 0.05 differential genes after few survive multiple-testing correction, including down FOLR3 and up TNNT1/MMP1; its data-sharing statement lists BioProject `PRJNA1456230` and a reviewer-token link. The paper's published biological claims are prior art, not our discovery; its PBMC assay and cohort definition are not interchangeable with the GSE45603 whole-blood array or GSE290313 whole-blood RNA-seq. Do not use the reviewer-token link as public access or publish/reuse its token.

At the checked time, NCBI E-utilities `esearch.fcgi` returned zero public IDs for this accession in `bioproject`, `sra` and `gds`; ENA Portal API returned `[]` for the public `study` and `read_run` queries. These are point-in-time **negative availability checks**, not proof that the data will never be released. No run/matrix was acquired or used, no sample or service count increment, and no independent validation/comparator score can be computed from the paper's unexposed raw data. The author-request route needs a separate communication decision; do not email or claim an accessible dataset based on a paper availability statement alone.

Public API checks:
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=PRJNA1456230&retmode=json`
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=PRJNA1456230&retmode=json`
- `https://www.ebi.ac.uk/ena/portal/api/search?result=study&query=study_accession%3D%22PRJNA1456230%22&fields=study_accession,secondary_study_accession,study_title&format=json`

Do not silently use this PBMC paper as a same-task benchmark for a whole-blood RNA/array result. A future comparison must align phenotype, sample timing, tissue, endpoints and training/test split, and distinguish authors' unadjusted exploratory DEG claims from a predictive metric.
