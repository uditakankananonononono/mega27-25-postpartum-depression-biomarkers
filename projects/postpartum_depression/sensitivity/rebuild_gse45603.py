"""Reconstruct the exposed 48-column cached matrix from exact GEO raw sources.

The historic max-mean probe selection used all 210 longitudinal columns before
selecting 48 postpartum samples. This is a provenance reproduction, not a safe
train-fold preprocessing recipe for a prospective classifier.
"""
import csv
import gzip
import hashlib
import io
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from ubiomark.geo import parse_series_matrix, _clean_symbol, to_gene_level

GEO = {
    'data/raw/matrix/GSE45603_series_matrix.txt.gz': ('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45603/matrix/GSE45603_series_matrix.txt.gz', 'da120f052d21b45df26bc0c913cafc14839589ef7416f75579f854d1f1a9e931'),
    'data/raw/gpl/GPL10558.annot.gz': ('https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz', 'c914fdbe1130906ce3b9c97f5a75280591c89c617c168fb687c641f683e76b45')
}


def rebuild(root=ROOT):
    root=Path(root)
    for rel, (_, expected) in GEO.items():
        actual=hashlib.sha256((root/rel).read_bytes()).hexdigest()
        if actual!=expected:
            raise ValueError(f'GEO raw hash mismatch: {rel}')
    expr, ann, _ = parse_series_matrix(str(root/'data/raw/matrix/GSE45603_series_matrix.txt.gz'))
    if expr.shape!=(15456,210) or not all(ann['Sample_platform_id']=='GPL10558'):
        raise ValueError('unexpected GEO matrix dimensions or platform')
    with gzip.open(root/'data/raw/gpl/GPL10558.annot.gz','rt',errors='replace') as f:
        rows=[line for line in f if not line.startswith(('^','!','#'))]
    annot=pd.read_csv(io.StringIO(''.join(rows)),sep='\t',dtype=str,low_memory=False)
    if not {'ID','Gene symbol'}<=set(annot):
        raise ValueError('required GPL annotation columns missing')
    p2s=annot.set_index('ID')['Gene symbol'].map(lambda s:_clean_symbol(s,'Gene symbol')).dropna()
    p2s.index=p2s.index.astype(str)
    crosswalk=list(csv.DictReader(open(root/'projects/postpartum_depression/sources/GSE45603_used_sample_crosswalk.csv',newline='')))
    selected=[r['gsm'] for r in crosswalk]
    if len(selected)!=48 or len(set(selected))!=48 or not set(selected)<=set(expr.columns):
        raise ValueError('selected GSMs invalid')
    # Reproduce old pipeline: max-mean probe is chosen from *all* 210 samples.
    # This must not be mistaken for train-fold-only feature construction.
    gene_all=to_gene_level(expr,p2s)
    matrix=gene_all[selected]
    cached=pd.read_pickle(root/'data/processed/postpartum_depression__GSE45603.pkl.gz')
    if set(matrix.index)!=set(cached.index) or set(matrix.columns)!=set(cached.columns):
        raise ValueError('reconstructed gene/GSM set differs from cached matrix')
    a=matrix.loc[cached.index,cached.columns].to_numpy()
    b=cached.to_numpy()
    if not np.array_equal(a,b,equal_nan=True):
        raise ValueError('reconstructed values differ from cached matrix')
    return {'raw_matrix_shape':list(expr.shape),'reconstructed_gene_shape':list(matrix.shape),
            'cached_shape':list(cached.shape),'exact_gene_and_value_match':True,
            'source_urls_and_sha256':{k:{'url':v[0],'sha256':v[1]} for k,v in GEO.items()},
            'processing_caveat':'Original max-mean probe selection used all 210 longitudinal samples, then selected 48 postpartum GSMs. This is a faithful reproduction of the exposed descriptive result, not train-only preprocessing for a prospective classifier.'}


if __name__=='__main__':
    r=rebuild()
    p=ROOT/'projects/postpartum_depression/sensitivity/gse45603_rebuild_audit.json'
    p.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
