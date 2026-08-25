"""Robustness of the PC to denominator + background choices (Ryan's questions, 2026-08-23):
 A. cytoplasm measure: mean (current) vs guard-band mean (exclude <0.5um around granules) vs median
 B. background: global p1 / p3 (current) / p5 / outside-worm (>10um from any nucleus) median
Whole-gonad, lamin envelopes, clean-13 gonads. Reports per-variant male HS-vs-noHS delta + p."""
import json, os, sys, re
import numpy as np, pandas as pd
from scipy import ndimage as ndi
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis"); sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import chload, crescent_axis as ca, tifffile

SP=chload.SP; SPACING=tuple(SP); VVOL=float(SP.prod())
SMOOTH_UM,BG_UM,K,VMIN,VMAX,CYTO=0.15,0.513,4.48,0.003,15.0,2.5
DB=np.arange(0.0,CYTO+1e-6,0.25)
OUT=r"C:/Users/ryane/coloc_analysis/denominator_check.csv"
CLEAN=[i for i in pd.read_csv(r"C:/Users/ryane/coloc_analysis/acquisition_metadata.csv").iid if not chload.is_excluded(i)]  # exclusions.json is the single source

def pc(syp,gran,outside,dt,bg,agg):
    num=w=0.0
    for lo,hi in zip(DB[:-1],DB[1:]):
        b=(dt>lo)&(dt<=hi); g=gran&b; o=outside&b
        if g.sum()<20 or o.sum()<50: continue
        den=(float(np.median(syp[o])) if agg=="median" else float(syp[o].mean()))-bg
        if den<=0: continue
        num+=g.sum()*((float(syp[g].mean())-bg)/den); w+=g.sum()
    return round(num/w,4) if w else None

rows=[]
if os.path.exists(OUT):
    rows=pd.read_csv(OUT).to_dict("records"); print(f"resuming, {len(rows)} rows already done", flush=True)
done={r["image_id"] for r in rows}
for iid in CLEAN:
    if iid in done: continue
    dapi,syp_c,lam,lab_c,fb=ca.load_crops(iid)
    rd=ca.find_run(iid); csv=os.path.join(rd,f"{iid}__nuclei.csv"); nc=pd.read_csv(csv)
    germ=[int(x) for x in nc[nc.in_germline.astype(bool)].nucleus_id]
    full=tifffile.imread(os.path.join(rd,f"{iid}__nuclei_labels.tif")); gn=np.isin(full,germ)
    obj=ndi.find_objects(gn.astype(np.uint8))[0]
    sl=tuple(slice(max(0,o.start-ca.PAD),min(dim,o.stop+ca.PAD)) for o,dim in zip(obj,gn.shape))
    ch=chload.load(ca.find_nd2(iid,csv)); pgl=ch["pgl"][sl].astype(np.float32); syp=ch["syp"][sl].astype(np.float32); del ch
    env=lam>0; dt=ndi.distance_transform_edt(~env,sampling=SPACING)
    cyto=ndi.binary_fill_holes(dt<=CYTO)&~env
    sm=ndi.gaussian_filter(pgl,sigma=SMOOTH_UM/SP); th=sm-ndi.gaussian_filter(sm,sigma=BG_UM/SP)
    rr=th[cyto]; neg=rr[rr<0]; noise=float(np.sqrt(np.mean(neg**2))) if neg.size>50 else float(np.std(rr))
    lab2,_=ndi.label((th>K*noise)&cyto); vol=np.bincount(lab2.ravel())*VVOL
    keep=np.where((vol>=VMIN)&(vol<=VMAX))[0]; keep=keep[keep!=0]
    gran=np.isin(lab2,keep)
    out_std=cyto&~gran
    guard=ndi.binary_dilation(gran,structure=np.ones((3,9,9),bool))   # ~0.3um z, 0.43um xy
    out_guard=cyto&~guard
    bgs={"p1":float(np.percentile(syp,1)),"p3":float(np.percentile(syp,3)),
         "p5":float(np.percentile(syp,5)),"outworm":float(np.median(syp[dt>10.0])) if (dt>10).sum()>1e4 else None}
    row={"image_id":iid,"sex":"male" if "_male" in iid.lower() else "herm",
         "treat":"HS" if re.search(r"_HS_",iid) else "noHS","n_gran":len(keep)}
    row["PC_current"]=pc(syp,gran,out_std,dt,bgs["p3"],"mean")
    row["PC_guard"]=pc(syp,gran,out_guard,dt,bgs["p3"],"mean")
    row["PC_median"]=pc(syp,gran,out_std,dt,bgs["p3"],"median")
    for name,b in bgs.items():
        if b is not None: row[f"PC_bg_{name}"]=pc(syp,gran,out_std,dt,b,"mean")
    rows.append(row); pd.DataFrame(rows).to_csv(OUT,index=False)
    print(f"{iid[-20:]}: cur={row['PC_current']} guard={row['PC_guard']} med={row['PC_median']} "
          f"bg1={row.get('PC_bg_p1')} bg5={row.get('PC_bg_p5')} bgout={row.get('PC_bg_outworm')}", flush=True)
print("DONE", flush=True)
