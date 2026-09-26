import csv
import sys
from pathlib import Path
import pandas as pd
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gse45603_control_scope import run


def test_real_source_exploratory_result(tmp_path):
    root = Path(__file__).resolve().parents[3]
    x = run(root / 'data/processed/postpartum_depression__GSE45603.pkl.gz',
            root / 'projects/postpartum_depression/sources/GSE45603_used_sample_crosswalk.csv',tmp_path)
    assert (x['case_people'],x['explicit_euthymic_people'],x['generic_control_people']) == (16,27,5)
    assert x['finite_gene_comparisons'] == 9744
    assert x['sign_flip_count'] == 601
    assert x['sign_flips_both_abs_g_at_least_0_1'] == 2
    assert x['sign_flips_both_abs_g_at_least_0_2'] == 0
    assert x['broad_q_below_0_05'] == x['strict_q_below_0_05'] == 0


def test_unknown_or_duplicate_source_identity_blocks(tmp_path):
    root = Path(__file__).resolve().parents[3]
    matrix = root / 'data/processed/postpartum_depression__GSE45603.pkl.gz'
    rows=list(csv.DictReader(open(root/'projects/postpartum_depression/sources/GSE45603_used_sample_crosswalk.csv')))
    rows[1]['gsm']=rows[0]['gsm']
    p=tmp_path/'bad.csv'
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    with pytest.raises(ValueError,match='duplicate GSM'):
        run(matrix,p,tmp_path)


def test_cached_matrix_rebuilds_exactly_from_raw_geo():
    from rebuild_gse45603 import rebuild
    x=rebuild()
    assert x['raw_matrix_shape']==[15456,210]
    assert x['reconstructed_gene_shape']==[9744,48]
    assert x['exact_gene_and_value_match']
