"""Audit sensitivity to generic-vs-explicit-euthymic GSE45603 controls.

Exploratory reanalysis of previously exposed cohort; no independent validation.
Input is the existing processed gene x GSM matrix, not a new discovery set.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'src'))
from ubiomark.stats import hedges_g, bh_fdr


def run(matrix_path, crosswalk, output_dir):
    matrix_path = Path(matrix_path)
    rows = list(csv.DictReader(open(crosswalk, newline='')))
    if len({r['gsm'] for r in rows}) != len(rows):
        raise ValueError('duplicate GSM in crosswalk')
    matrix = pd.read_pickle(matrix_path)
    if not matrix.index.is_unique or not matrix.columns.is_unique:
        raise ValueError('matrix gene or GSM identifiers are not unique')
    if set(matrix.columns) != {r['gsm'] for r in rows}:
        raise ValueError('matrix GSMs do not exactly match crosswalk')
    cases = [r['gsm'] for r in rows if r['source_condition'] == 'condition: PPD']
    euthymic = [r['gsm'] for r in rows if r['source_condition'] == 'condition: euthymic']
    generic = [r['gsm'] for r in rows if r['source_condition'] == 'condition: control']
    if len(cases) < 3 or len(euthymic) < 3 or len(generic) < 1:
        raise ValueError('insufficient or missing phenotype stratum')
    if len(cases) + len(euthymic) + len(generic) != len(rows):
        raise ValueError('unknown source condition')
    if len({r['person_token'] for r in rows}) != len(rows):
        raise ValueError('person token repeats across selected rows')
    if any(r['timepoint'] != 'time: postpartum' for r in rows):
        raise ValueError('non-postpartum sample')
    comparisons = {}
    for name, controls in [('broad', euthymic + generic), ('strict_euthymic', euthymic)]:
        g, v = hedges_g(matrix[cases].to_numpy(float), matrix[controls].to_numpy(float))
        comparisons[name] = (g, v)
    bg, bv = comparisons['broad']; sg, sv = comparisons['strict_euthymic']
    valid = np.isfinite(bg) & np.isfinite(sg) & np.isfinite(bv) & np.isfinite(sv)
    b_z, s_z = bg[valid]/np.sqrt(bv[valid]), sg[valid]/np.sqrt(sv[valid])
    from scipy.stats import norm
    bq = bh_fdr(2*norm.sf(np.abs(b_z)));sq = bh_fdr(2*norm.sf(np.abs(s_z)))
    both_nonzero = (bg[valid]!=0) & (sg[valid]!=0)
    flips = both_nonzero & (np.sign(bg[valid]) != np.sign(sg[valid]))
    out = pd.DataFrame({'gene':matrix.index[valid], 'broad_g':bg[valid], 'broad_v':bv[valid],
                        'strict_g':sg[valid], 'strict_v':sv[valid], 'broad_q':bq,
                        'strict_q':sq, 'sign_flip':flips})
    out = out.sort_values('gene', kind='mergesort')
    output_dir = Path(output_dir); output_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_dir / 'gse45603_control_scope_effects.csv', index=False)
    summary = {'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45603',
        'matrix_path_argument': str(matrix_path), 'matrix_sha256': hashlib.sha256(matrix_path.read_bytes()).hexdigest(),
        'case_people':len(cases), 'explicit_euthymic_people':len(euthymic),
        'generic_control_people':len(generic), 'finite_gene_comparisons':int(valid.sum()),
        'sign_flip_count':int(flips.sum()), 'sign_flip_rate_of_nonzero_genes':float(flips.sum()/both_nonzero.sum()),
        'broad_q_below_0_05':int((bq<.05).sum()), 'strict_q_below_0_05':int((sq<.05).sum()),
        'caveat':'Prior GSE45603 outcomes already exposed. Broad versus strict is exploratory source-label sensitivity, not an independent test or new PPD biomarker.'}
    (output_dir / 'gse45603_control_scope_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--matrix',required=True);p.add_argument('--crosswalk',default=str(HERE.parent/'sources'/'GSE45603_used_sample_crosswalk.csv'));p.add_argument('--output-dir',default=str(HERE))
    args=p.parse_args();print(json.dumps(run(args.matrix,args.crosswalk,args.output_dir),indent=2))
