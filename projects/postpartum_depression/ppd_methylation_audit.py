#!/usr/bin/env python3
"""Source-level and matrix-level audit of two GEO PPD methylation studies; no clinical claims."""
import csv,gzip,hashlib,json,pathlib,os,urllib.request
import numpy as np
from scipy.stats import ttest_ind
R=pathlib.Path(__file__).resolve().parent/'sources'
CACHE=pathlib.Path(os.getenv('PPD_GEO_CACHE', '/tmp/ppd-geo-cache'))
CACHE.mkdir(parents=True,exist_ok=True)
EXPECTED_SHA256={'GSE44132':'b2d715ec45c250426c8b2a4e5317bf589f581f8425d1101d9f6d752ba6dc19b9','GSE335141':'13002634be22d3b58529402b07b89a57610479f18f3208ebd8ce94c6fc12b5b2'}
URLS={'GSE44132':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44132/matrix/GSE44132_series_matrix.txt.gz', 'GSE335141':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE335nnn/GSE335141/suppl/GSE335141_GEO_Upload_normalized_beta_matrix.tsv.gz'}
for acc,filename in [('GSE44132','GSE44132_series_matrix.txt.gz'),('GSE335141','GSE335141_beta_matrix.tsv.gz')]:
 p=CACHE/filename;cross=R/f'{acc}_source_crosswalk.csv'; records=list(csv.DictReader(cross.open()))
 assert len(records)==len({x['gsm'] for x in records}), 'GSM duplicate'
 for x in records:
  assert x['series']==acc and x['source_url'].startswith(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={x["gsm"]}')
 if not p.exists():
  with urllib.request.urlopen(URLS[acc],timeout=180) as response,p.open('wb') as out:
   import shutil;shutil.copyfileobj(response,out)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==EXPECTED_SHA256[acc], f'{acc}: source matrix digest differs; stop before analysis'
 if acc=='GSE44132':
  # GEO matrix columns are GSM; exclude all five PR01-084 cross-batch technical controls.
  with gzip.open(p,'rt') as f:
   for line in f:
    if line.startswith('"ID_REF"'):
     header=next(csv.reader([line],delimiter='\t'));break
  colmap={x['gsm']:x['gsm'] for x in records}
  keep=[x for x in records if x['patient_token']!='PR01-084']
  case=[x for x in keep if x['phenotype']=='yes'];control=[x for x in keep if x['phenotype']=='no']
  selected=case+control
 else:
  assert all(x.get('matrix_column') for x in records), 'Stored verified GSM-to-array map missing'
  with gzip.open(p,'rt') as f:header=f.readline().strip().split('\t')
  colmap={x['gsm']:x['matrix_column'] for x in records}
  case=[x for x in records if x['phenotype']=='PPD' and x['timepoint']=='T0']
  control=[x for x in records if x['phenotype']=='healthy' and x['timepoint']=='T0']
  selected=case+control
 assert len(header)==len(set(header)),(acc,'matrix duplicate columns')
 assert set(colmap.values())==set(header[1:]),(acc,'source-matrix mismatch',len(set(colmap.values())),len(header)-1)
 assert len({x['patient_token'] for x in selected})==len(selected)
 idx=[header.index(colmap[x['gsm']]) for x in selected]
 nvalid=0; nprobes=0; hits=[]; batch=[]; probeids=[]
 def flush_batch():
  global nvalid
  if not batch:return
  x=np.asarray(batch,dtype=float)
  good=np.isfinite(x).all(axis=1)&(x>=0).all(axis=1)&(x<=1).all(axis=1)
  nvalid+=int(good.sum()); x=x[good]; ids=[probeids[i] for i in np.flatnonzero(good)]
  if len(x):
   a=x[:,:len(case)];b=x[:,len(case):]; _,pvals=ttest_ind(a,b,equal_var=False,axis=1)
   diffs=a.mean(axis=1)-b.mean(axis=1)
   hits.extend((pid,float(d),float(pv)) for pid,d,pv in zip(ids,diffs,pvals) if np.isfinite(pv))
  batch.clear();probeids.clear()
 with gzip.open(p,'rt',newline='') as f:
  if acc=='GSE44132':
   for line in f:
    if line.startswith('"ID_REF"'):break
  rd=csv.reader(f,delimiter='\t')
  for line in rd:
   if not line or line[0].startswith('!series_matrix_table_end'):break
   nprobes+=1
   if len(line)!=len(header):continue
   try:vals=[float(line[j]) for j in idx]
   except ValueError:continue
   batch.append(vals);probeids.append(line[0]);
   if len(batch)>=5000:flush_batch()
   if nprobes%200000==0:print(acc,'parsed',nprobes,flush=True)
 flush_batch()
 from scipy.stats import false_discovery_control
 assert hits, f'{acc}: no finite probe tests'
 pvals=np.array([h[2] for h in hits]);qvals=false_discovery_control(pvals,method='bh')
 top=sorted(zip(hits,qvals),key=lambda x:(x[1],x[0][2]))[:30]
 with (CACHE/f'{acc}_exploratory_probes_reproduced.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['probe','case_minus_control_mean_beta','welch_p','bh_q']);w.writerows([(h[0],h[1],h[2],float(q)) for h,q in top])
 summary={'series':acc,'matrix_url':URLS[acc],'matrix_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'matrix_columns':len(header)-1,'GSM_records':len(records),'distinct_patient_tokens':len(set(x['patient_token'] for x in records)),'analysis_unit':'one person, one T0 draw' if acc=='GSE335141' else 'one person, antenatal draw; technical PR01-084 controls excluded','analyzed_case':len(case),'analyzed_control':len(control),'probes_seen':nprobes,'probes_complete_finite_beta':nvalid,'probes_tested':len(hits),'fdr_0.05':int(sum(qvals<.05)),'top_probe':top[0][0] if top else None,'caveat':'Exploratory in-cohort probe differences only; no held-out clinical prediction, no cross-study methylation validation, no new disease marker claim.'}
 (CACHE/f'{acc}_analysis_summary_reproduced.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
