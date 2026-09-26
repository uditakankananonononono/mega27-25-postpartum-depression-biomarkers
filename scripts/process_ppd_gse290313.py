"""Registered P6: postpartum-only GSE290313 count files, no DE-prefix grouping."""
import gzip, tarfile, sys
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,'src')
from ubiomark import geo,stats
G='GSE290313';_,a,_=geo.parse_series_matrix(geo.download_matrices(G)[0])
t=a.Sample_characteristics_ch1_5.str.extract(r': (.*)')[0].str.replace('postpartumy','postpartum',regex=False)
g=a.Sample_characteristics_ch1_4.str.extract(r': (.*)')[0]
ss=a.Sample_characteristics_ch1_2.str.extract(r': (.*)')[0]
case=(t=='postpartum sample')&(g=='Depressive symptomes only postpartum')&(ss=='No')
ctrl=(t=='postpartum sample')&(g=='Control')&(ss=='No')
lab=pd.Series(np.where(case,'case',np.where(ctrl,'control','excluded')),index=a.index)
assert int((t=='postpartum sample').sum())==184
assert a.index.is_unique and a.Sample_description.is_unique
selected=a.index[lab!='excluded'];assert len(selected)>20 and lab.loc[selected].nunique()==2
h=pd.read_csv(geo.HGNC_PATH,sep='\t',dtype=str,usecols=['symbol','ensembl_gene_id']).dropna()
mapping=dict(zip(h.ensembl_gene_id,h.symbol))
files={}
with tarfile.open('/tmp/GSE290313_RAW.tar') as tar:
    for z in tar.getmembers():
        gsm=z.name.split('_')[0]
        if gsm not in selected:continue
        with gzip.GzipFile(fileobj=tar.extractfile(z)) as f:
            v=pd.read_csv(f,sep='\t',header=None,names=['gene','count'])
        v['gene']=v.gene.str.split('.').str[0].map(mapping)
        v=v.dropna(subset=['gene']).groupby('gene')['count'].sum()
        files[gsm]=v
assert set(files)==set(selected),(len(files),len(selected))
x=pd.DataFrame(files).fillna(0).astype(float)
assert np.isfinite(x.values).all() and (x.values>=0).all()
print('selected',lab.value_counts().to_dict(),'genes mapped',len(x),'full library/GSM unique',len(selected))
cpm=x.div(x.sum(0),axis=1)*1e6
keep=(cpm>1).mean(1)>=.2;z=np.log2(cpm.loc[keep]+1).astype('float32')
ca=lab.index[lab=='case'];co=lab.index[lab=='control']
effect,var=stats.hedges_g(z[ca].to_numpy(float),z[co].to_numpy(float))
Path('results/series').mkdir(exist_ok=True)
pd.DataFrame({'gene':z.index,'g':effect,'v':var}).dropna().to_csv('results/series/postpartum_depression__GSE290313.csv.gz',index=False)
pd.DataFrame({'gsm':a.index,'timepoint':t.values,'trajectory':g.values,'ssri':ss.values,'label':lab.values}).to_csv('results/series/postpartum_depression__GSE290313.labels.csv',index=False)
Path('data/processed').mkdir(exist_ok=True)
z.to_pickle('data/processed/postpartum_depression__GSE290313.pkl.gz')
print('effect genes',int(np.isfinite(effect).sum()),'case',len(ca),'control',len(co))
