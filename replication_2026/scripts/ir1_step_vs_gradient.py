"""Is IR1 unit count a STEP at the panel boundary or a CONTINUOUS function of ancestry?

Host selection predicts a gradient: IR1 should track genotype PCs continuously, including WITHIN
a continental group. A per-panel virus stock predicts a step: the continental/population label
should capture the variance and ancestry PCs should add nothing within a group.

Uses MAGE's published genotype PCs (same PCs used for its eQTL mapping).
"""
import csv, glob, math, gzip, statistics as st
from collections import defaultdict
BASE='/work/nfp8/EBV'; RU=22737/3072.0

ir1={}
for f in glob.glob(f'{BASE}/ir1/*.tsv'):
    r=list(csv.DictReader(open(f), delimiter='\t'))
    if r:
        try:
            v=float(r[0]['ir1_copies_vs_unique'])
            if v>0: ir1[r[0]['sample']]=v*RU
        except (ValueError,KeyError): pass
wt={}
for r in csv.DictReader(open(f'{BASE}/ebv_cn_1kGP_2504.tsv'), delimiter='\t'):
    try: wt[r['sample']]=float(r['wt_fraction'])
    except ValueError: pass
coinf={s for s,v in wt.items() if v>0.01}
meta={m['sample_kgpID']:(m['population'],m['continentalGroup']) for m in csv.DictReader(
    open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sample.metadata.MAGE.v1.0.txt'),
    delimiter='\t')}
Z='/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data'
with gzip.open(f'{Z}/QTL_results/eQTL_results/expression_quants/eQTL_covariates.tab.gz','rt') as fh:
    rows=[l.rstrip('\n').split('\t') for l in fh]
ids=rows[0][1:]; cov={r[0]:r[1:] for r in rows[1:]}
PCS=['PC1','PC2','PC3','PC4','PC5']

ss=[s for s in ids if s in ir1 and s not in coinf and s in meta]
idx={s:ids.index(s) for s in ss}
y=[ir1[s] for s in ss]
X=[[float(cov[p][idx[s]]) for p in PCS] for s in ss]
print(f'n = {len(ss)} lines with IR1 + genotype PCs\n')

def ols_r2(X, y):
    """Multiple R^2 via normal equations with an intercept."""
    k=len(X[0]); n=len(y)
    A=[[1.0]+row for row in X]; p=k+1
    XtX=[[sum(A[i][a]*A[i][b] for i in range(n)) for b in range(p)] for a in range(p)]
    Xty=[sum(A[i][a]*y[i] for i in range(n)) for a in range(p)]
    M=[XtX[i][:]+[Xty[i]] for i in range(p)]
    for c in range(p):
        piv=max(range(c,p), key=lambda r: abs(M[r][c]))
        if abs(M[piv][c])<1e-12: continue
        M[c],M[piv]=M[piv],M[c]
        d=M[c][c]
        M[c]=[v/d for v in M[c]]
        for r in range(p):
            if r!=c and abs(M[r][c])>1e-15:
                f=M[r][c]; M[r]=[v-f*w for v,w in zip(M[r],M[c])]
    beta=[M[i][p] for i in range(p)]
    pred=[sum(b*v for b,v in zip(beta,A[i])) for i in range(n)]
    my=st.mean(y)
    ssr=sum((a-b)**2 for a,b in zip(y,pred)); sst=sum((a-my)**2 for a in y)
    return 1-ssr/sst if sst else 0.0

# continental indicator as dummies
cgs=sorted({meta[s][1] for s in ss})
D=[[1.0 if meta[s][1]==c else 0.0 for c in cgs[1:]] for s in ss]
print(f'variance in IR1 explained:')
print(f'  genotype PC1-5 (continuous ancestry) : R2 = {ols_r2(X,y):.4f}')
print(f'  continental-group dummies (step)     : R2 = {ols_r2(D,y):.4f}')
print(f'  both                                 : R2 = {ols_r2([a+b for a,b in zip(X,D)],y):.4f}')

# the decisive part: WITHIN Europe, do PCs predict IR1?
for grp in cgs:
    sub=[s for s in ss if meta[s][1]==grp]
    if len(sub)<40: continue
    ysub=[ir1[s] for s in sub]
    Xsub=[[float(cov[p][idx[s]]) for p in PCS] for s in sub]
    popd=sorted({meta[s][0] for s in sub})
    Dsub=[[1.0 if meta[s][0]==p else 0.0 for p in popd[1:]] for s in sub]
    print(f'\n  within {grp} (n={len(sub)}, {len(popd)} populations): '
          f'PCs R2 = {ols_r2(Xsub,ysub):.4f}   population dummies R2 = {ols_r2(Dsub,ysub):.4f}')
