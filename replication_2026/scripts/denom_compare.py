"""Does the robust 15 Mb / 150-window denominator change CN versus the pooled 2.4 Mb one?

Also builds a per-sample GC-bias metric from the per-window depths (highest-GC vs
lowest-GC windows) and asks whether CN depends on it. EBV is 59.4% GC vs 40.8% for the
autosomal panel, so GC bias would act differentially on numerator and denominator.
"""
import csv, glob, math, os, statistics as st

gc = {}
with open('/work/nfp8/EBV/metadata/window_gc.tsv') as fh:
    for line in fh:
        reg, v = line.split('\t')
        c, rng = reg.split(':'); s = int(rng.split('-')[0]) - 1
        gc[f'{c}_{s//100000}'] = float(v)

old = {}
for f in glob.glob('/work/nfp8/EBV/cn/*.tsv'):
    with open(f) as fh:
        r = list(csv.DictReader(fh, delimiter='\t'))
    if r:
        try:
            old[r[0]['sample']] = (float(r[0]['auto_depth']), float(r[0]['ebv_depth']),
                                   float(r[0]['ebv_cn_depth_ratio']))
        except ValueError: pass

lo_k = [k for k, _ in sorted(gc.items(), key=lambda kv: kv[1])[:30]]
hi_k = [k for k, _ in sorted(gc.items(), key=lambda kv: kv[1])[-30:]]
print(f'GC-bias metric: mean depth of 30 highest-GC windows '
      f'({100*st.mean([gc[k] for k in hi_k]):.1f}% GC) / 30 lowest '
      f'({100*st.mean([gc[k] for k in lo_k]):.1f}% GC)\n')

rows = []
for f in glob.glob('/work/nfp8/EBV/denom/*.windows.tsv'):
    s = os.path.basename(f).replace('.windows.tsv', '')
    if s not in old: continue
    d = {}
    with open(f) as fh:
        for line in fh:
            k, n, dep = line.split('\t')
            if int(n) > 50000: d[k] = float(dep)
    if len(d) < 140: continue
    med = st.median(d.values())
    hi = [d[k] for k in hi_k if k in d]; lo = [d[k] for k in lo_k if k in d]
    rows.append((s, old[s][0], med, old[s][1], old[s][2],
                 st.mean(hi)/st.mean(lo) if lo and hi else float('nan'),
                 st.pstdev(d.values())/med))

if not rows:
    print('no overlapping samples yet'); raise SystemExit

print(f'n = {len(rows)} samples with both denominators\n')
print(f'{"sample":<10} {"old_24w":>8} {"new_150w":>9} {"ratio":>7} {"CN_old":>8} {"CN_new":>8} {"GCbias":>7} {"winCV":>6}')
diffs, cno, cnn = [], [], []
for s, o, m, ebv, c_old, gcb, cv in sorted(rows):
    c_new = 2 * ebv / m
    diffs.append(m / o); cno.append(c_old); cnn.append(c_new)
    print(f'{s:<10} {o:>8.2f} {m:>9.2f} {m/o:>7.3f} {c_old:>8.2f} {c_new:>8.2f} {gcb:>7.3f} {cv:>6.3f}')

print(f'\ndenominator ratio new/old: mean {st.mean(diffs):.4f}  sd {st.pstdev(diffs):.4f} '
      f' range {min(diffs):.3f}-{max(diffs):.3f}')
mo, mn = st.mean(cno), st.mean(cnn)
do, dn = [v-mo for v in cno], [v-mn for v in cnn]
r = sum(a*b for a,b in zip(do,dn))/math.sqrt(sum(a*a for a in do)*sum(b*b for b in dn))
print(f'CN_old vs CN_new Pearson r = {r:.5f}')
print(f'mean |relative CN change| = {100*st.mean([abs(a-b)/b for a,b in zip(cno,cnn)]):.2f}%')
gcs = [r_[5] for r_ in rows]
print(f'\nGC-bias metric across samples: mean {st.mean(gcs):.3f}  sd {st.pstdev(gcs):.3f}')
print(f'within-sample per-window CV: mean {st.mean([r_[6] for r_ in rows]):.3f}')
