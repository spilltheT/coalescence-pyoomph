import numpy as np, os, re, sys, glob
if len(sys.argv)!=4:
    sys.exit(f"usage: {os.path.basename(sys.argv[0])} <base_dir> <theta_deg> <out.npz>")
base=sys.argv[1]; theta_deg=float(sys.argv[2]); out=sys.argv[3]
theta=theta_deg*np.pi/180.0; HP=1e-4
dom=os.path.join(base,"domain"); files=sorted(glob.glob(os.path.join(dom,"domain_*.txt")))
if not files:
    sys.exit(f"no domain_*.txt snapshots found in {dom!r}; check the run path / that the simulation has produced output")
def sub(a,n=600):
    a=np.asarray(a)
    return a if len(a)<=n else a[np.linspace(0,len(a)-1,n).astype(int)]
def hdr_time(fp):
    with open(fp) as f: h=f.readline()
    return float(re.search(r"@time=([0-9.eE+-]+)",h).group(1))
def load(fp):
    t=hdr_time(fp); d=np.loadtxt(fp,skiprows=1); x=d[:,0]; h=d[:,1]
    s=np.argsort(x); return t,x[s],h[s]
times=np.array([hdr_time(fp) for fp in files]); o=np.argsort(times)
files=[files[i] for i in o]; times=times[o]
nearest=lambda tt:int(np.argmin(np.abs(times-tt)))
prof_times=[0.3,1,3,10,30,60,100]
P=[load(files[nearest(tt)]) for tt in prof_times if times.max()>=tt*0.99]
H=[]
for fp in files:
    t,x,h=load(fp); m=np.abs(x)<1.0; H.append((t,float(h[m].min())))
H=np.array(H); t_all=H[:,0]; h0_all=H[:,1]
ratio_all=np.where(t_all>0,h0_all/(theta**4*t_all),np.inf)
sel=np.where((ratio_all>=0.27)&(ratio_all<=0.30)&(h0_all>3*HP)&(t_all>0))[0]
if len(sel)>14:
    lt=np.log10(t_all[sel]); g=np.linspace(lt.min(),lt.max(),14)
    sel=np.unique([sel[int(np.argmin(np.abs(lt-q)))] for q in g])
C=[]
for i in sel:
    t,x,h=load(files[i]); m1=np.abs(x)<1.0
    j=np.argmin(h[m1]); x0=x[m1][j]; h0=h[m1][j]
    xi=(x-x0)*theta/h0; Hn=h/h0
    m=np.abs(xi)<=5.5; s=np.argsort(xi[m])
    C.append((t,float(ratio_all[i]),xi[m][s],Hn[m][s]))
np.savez(out, theta_deg=theta_deg,
  prof_t=np.array([p[0] for p in P]),
  prof_x=np.array([sub(p[1]) for p in P],dtype=object),
  prof_h=np.array([sub(p[2]) for p in P],dtype=object),
  coll_t=np.array([c[0] for c in C]),
  coll_ratio=np.array([c[1] for c in C]),
  coll_xi=np.array([sub(c[2]) for c in C],dtype=object),
  coll_H=np.array([sub(c[3]) for c in C],dtype=object),
  h0t=H)
print(f"wrote {out} | nfiles {len(files)} | tmax {times.max():.1f} | coll snaps {len(C)} | "
      f"t[{(C[0][0] if C else 0):.2f}..{(C[-1][0] if C else 0):.1f}]")
