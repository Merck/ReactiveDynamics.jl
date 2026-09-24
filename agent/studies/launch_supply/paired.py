import sys; sys.argv=[sys.argv[0]]; sys.path.insert(0,'.')
from analyze import runs, W
def npv(a,w,s): return runs[a][w][s]['npv']
print("both-world per-seed difference vs its wait reference, and the odds-weighted delta if that path were the both-world path (other worlds: seed 1)")
for a,ref in [("arm_cmo","model"),("arm_cmo_tern_only","model"),("arm_build_now","model"),("arm_build_after","model"),("arm_cmo_no2","model_no2"),("arm_cmo_tern_only_no2","model_no2")]:
    rb=ref.replace("_no2","") ; ab=a.replace("_no2","")
    rest=sum(W[w]*(npv(ab,w,"s1")-npv(rb,w,"s1")) for w in ("konly","tonly","neither"))
    ds=[npv(a,"both",s)-npv(ref,"both",s) for s in ("s1","s2","s3")]
    print(f"  {a:24s} vs {ref:10s} both-world d s1..s3 {ds[0]:+7.1f} {ds[1]:+7.1f} {ds[2]:+7.1f} | weighted per path {W['both']*ds[0]+rest:+6.1f} {W['both']*ds[1]+rest:+6.1f} {W['both']*ds[2]+rest:+6.1f}")
pv5=1.09**(-14/12)*5
print(f"walk5: PV of $5M at t=14 {pv5:.3f}; arm_cmo +{pv5*W['neither']:.2f} (P doomed = neither {W['neither']:.2f}); arm_cmo_tern_only +{pv5*(W['konly']+W['neither']):.2f} (P doomed = Tern fails {W['konly']+W['neither']:.2f})")
