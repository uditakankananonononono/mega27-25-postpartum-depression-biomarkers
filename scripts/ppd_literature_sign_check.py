"""Descriptive cross-study lookup of already-published PPD gene claims; no new validation."""
import pandas as pd
reported={'FOLR3':'down','TNNT1':'up','MMP1':'up'}
files=['results/series/postpartum_depression__GSE45603.csv.gz','results/series/postpartum_depression__GSE290313-PP.csv.gz']
rows=[]
for f in files:
 x=pd.read_csv(f).set_index('gene')
 for gene,direction in reported.items():
  if gene not in x.index:
   rows.append(dict(cohort=f.split('__')[-1].split('.')[0],gene=gene,published_direction=direction,measured=False))
  else:
   r=x.loc[gene];rows.append(dict(cohort=f.split('__')[-1].split('.')[0],gene=gene,published_direction=direction,measured=True,g=float(r.g),se=float(r.v**.5),same_sign=(r.g<0 if direction=='down' else r.g>0)))
pd.DataFrame(rows).to_csv('results/ppd_published_gene_context.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
