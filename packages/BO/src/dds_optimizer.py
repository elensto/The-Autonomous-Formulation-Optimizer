"""DGEBA/4,4'-DDS Bayesian target-window optimizer.

Adapts Nannan's Gaussian-process surrogate and target-interval acquisition;
separates (a) existing cured-network selection from (b) stoichiometry-and-cure
requests that still require molecular network generation.

Never trains on literature/synthetic labels in real mode. Synthetic mode is
allowed ONLY for labeled EXAMPLE fixtures with --allow-synthetic-demo.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE/'reference'))
from nannan_original import GPSurrogate, log_prob_in_window

APPROVED = 'real_classical_MD'
FEATURES = {'conversion':['conversion_pct'], 'composition':['ratio','conversion_pct']}
REQUIRED_OBS = ['observation_id','candidate_id','n_DGEBA','n_DDS','ratio','conversion_pct','tg_C','tg_se_C','physcheck_decision','bo_trainable','source_type','experiment_status']
REQUIRED_CAND = ['candidate_id','structure_id','n_DGEBA','n_DDS','ratio','conversion_pct','structure_status','file_path','physcheck_decision','job_status']
ACTIVE = {'queued','running','submitted','pending','reserved'}
COMPLETE = {'completed','done','accepted'}


def read_csv(path, columns):
    df=pd.read_csv(path,keep_default_na=False)
    missing=set(columns)-set(df.columns)
    if missing: raise ValueError(f'{path}: missing columns {sorted(missing)}')
    return df


def ratio_from_counts(n_dgeba, n_dds):
    """DDS has 4 reactive amine H; DGEBA has 2 epoxy functions."""
    if int(n_dgeba)<=0 or int(n_dds)<=0: raise ValueError('Positive integer counts required')
    return 4*int(n_dds)/(2*int(n_dgeba))


def verify_counts(df, name):
    for _,r in df.iterrows():
        n1=float(r.n_DGEBA);n2=float(r.n_DDS)
        if not(n1.is_integer() and n2.is_integer() and n1>0 and n2>0):
            raise ValueError(f'{name}: invalid molecule counts in {r.candidate_id}')
        expected=ratio_from_counts(n1,n2)
        if not np.isfinite(float(r.ratio)) or abs(float(r.ratio)-expected)>1e-5:
            raise ValueError(f'{name}: incorrect stoichiometry {r.candidate_id} (expected {expected:.6g})')
        if not 0<=float(r.conversion_pct)<=100:
            raise ValueError(f'{name}: impossible conversion for {r.candidate_id}')


def accepted_observations(df, allow_synthetic_demo=False):
    permitted={APPROVED}
    if allow_synthetic_demo: permitted.add('SYNTHETIC_DEMO_NOT_MD')
    df=df.copy()
    ok=(df.source_type.isin(permitted)
        & df.physcheck_decision.astype(str).str.upper().eq('ACCEPT')
        & df.bo_trainable.astype(str).str.lower().eq('true')
        & df.experiment_status.astype(str).str.lower().isin(COMPLETE))
    a=df.loc[ok].copy()
    if len(a):
        verify_counts(a,'observations')
        for key in ['ratio','conversion_pct','tg_C','tg_se_C']:
            a[key]=pd.to_numeric(a[key],errors='coerce')
        a=a.replace([np.inf,-np.inf],np.nan).dropna(subset=['ratio','conversion_pct','tg_C','tg_se_C'])
        a=a[a.tg_se_C>=0]
    return a


def build_composition_pool(n_dgeba=100,dds_counts=(40,45,50,55,60),conversions=(60,75,85,90),ratio_min=.8,ratio_max=1.2):
    """Planning grid only: each point requires actual packing, reaction and equilibration."""
    rows=[]
    for nd in dds_counts:
        ratio=ratio_from_counts(n_dgeba,nd)
        if ratio<ratio_min-1e-9 or ratio>ratio_max+1e-9:continue
        # For any starting stoichiometry, maximum attainable conversion of
        # INITIAL EPOXIDES cannot exceed min(1, amine-H/epoxide ratio).
        max_epoxy_conv=min(1.,ratio)*100
        for a in conversions:
            if a>max_epoxy_conv+1e-9:continue
            id=f'DDS_N{n_dgeba}_H{nd}_A{a:03}'
            rows.append({'candidate_id':id,'structure_id':'','n_DGEBA':n_dgeba,'n_DDS':nd,
                         'ratio':round(ratio,6),'conversion_pct':a,'structure_status':'NEEDS_NETWORK_GENERATION',
                         'file_path':'','physcheck_decision':'NOT_EVALUATED','job_status':'not_started',
                         'request_type':'generate_and_cure_network', 'source_type':'proposed_design_not_observation'})
    return pd.DataFrame(rows)


def eligible_candidates(df,mode):
    d=df.copy()
    if mode=='conversion':
        ok=(d.structure_status.astype(str).str.upper().eq('CURED_NETWORK_VALIDATED')
            & d.physcheck_decision.astype(str).str.upper().eq('ACCEPT')
            & d.file_path.astype(str).str.len().gt(0))
        d=d.loc[ok].copy()
        # A path must be real (no placeholder file names)
        d=d[d.file_path.apply(lambda p: Path(p).is_file() or (HERE / str(p)).is_file())]
        if len(d):
            # avoid mixing stoichiometries when optimizing only conversion
            pairs=d[['n_DGEBA','n_DDS']].drop_duplicates()
            if len(pairs)!=1:raise ValueError('Conversion mode requires single fixed stoichiometry; filter manifest first')
    else:
        # Composition requests may be for structures not yet generated.
        d=d[d.structure_status.isin(['NEEDS_NETWORK_GENERATION','CURED_NETWORK_VALIDATED'])].copy()
    return d


def remove_duplicates(d,observed,all_observations, mode):
    if d.empty:return d.reset_index(drop=True)
    used=set(all_observations.loc[all_observations.experiment_status.astype(str).str.lower().isin(ACTIVE|COMPLETE),'candidate_id'].astype(str))
    # Also check validated observations, even if status inconsistent
    used.update(observed.candidate_id.astype(str))
    d=d[~d.candidate_id.isin(used)]
    d=d[~d.job_status.astype(str).str.lower().isin(ACTIVE|COMPLETE)]
    if mode=='conversion':
        # Conversions are distinct design locations; replicas can be studied
        # separately after an explicit revisit policy is added.
        used_conversion=set(observed.conversion_pct.astype(float).round(5))
        d=d[~d.conversion_pct.astype(float).round(5).isin(used_conversion)]
        d=d.sort_values(['conversion_pct','candidate_id']).drop_duplicates('conversion_pct',keep='first')
    else:
        measured_keys=set(zip(observed.ratio.astype(float).round(5),observed.conversion_pct.astype(float).round(5)))
        d=d[[tuple(k) not in measured_keys for k in zip(d.ratio.astype(float).round(5),d.conversion_pct.astype(float).round(5))]]
        d=d.drop_duplicates(subset=['ratio','conversion_pct'])
    return d.reset_index(drop=True)


def farthest_first(candidates,obs,features,n_points):
    """Initial deterministic design selection for zero/few valid MD labels."""
    allx=candidates[features].to_numpy(float)
    scales=np.ptp(allx,axis=0)
    scales=np.where(scales>0,scales,1.)
    X=allx/scales
    anchor=obs[features].to_numpy(float)/scales if len(obs) else np.empty((0,len(features)))
    picked=[]
    for i in range(min(n_points,len(candidates))):
        refs=np.vstack([anchor,X[picked]]) if picked else anchor
        if len(refs):
            distance=np.min(np.linalg.norm(X[:,None,:]-refs[None,:,:],axis=2),axis=1)
        else:
            # begin near reference stoichiometry or center conversion for a sensible seed
            mid=np.mean(X,axis=0)
            distance=-np.linalg.norm(X-mid,axis=1)
        distance[picked]=-np.inf
        j=int(np.argmax(distance));picked.append(j)
    return picked


def propose(candidates,all_obs,mode,tg_min=150.,tg_max=220.,batch_size=3,seed=42,
            allow_synthetic_demo=False,min_gp=4):
    if mode not in FEATURES:raise ValueError('mode must be conversion/composition')
    if tg_max<=tg_min:raise ValueError('Tg interval reversed')
    verify_counts(candidates,'candidates')
    obs=accepted_observations(all_obs,allow_synthetic_demo)
    if mode=='conversion':
        if len(candidates):
            first=candidates.iloc[0]
            obs=obs[(obs.n_DGEBA.astype(int)==int(first.n_DGEBA))&(obs.n_DDS.astype(int)==int(first.n_DDS))]
    usable=eligible_candidates(candidates,mode)
    usable=remove_duplicates(usable,obs,all_obs,mode)
    if not len(usable):
        return pd.DataFrame(),{'status':'NO_ELIGIBLE_CANDIDATES','mode':mode,'n_approved_observations':len(obs),'message':('Import and validate cured network files before conversion BO' if mode=='conversion' else 'No untested formulation designs remain')}
    features=FEATURES[mode]
    n_unique=len(obs[features].drop_duplicates()) if len(obs) else 0
    method='diversity_exploration';meta={'mode':mode,'n_approved_observations':len(obs),'n_distinct_training_points':n_unique,'n_candidates_eligible':len(usable),'target_C':[tg_min,tg_max],'source_mode':'SYNTHETIC_DEMO_NOT_MD' if allow_synthetic_demo else 'real_MD_only'}
    if n_unique>=min_gp:
        training=obs.groupby(features,as_index=False).agg({'tg_C':'mean','tg_se_C':lambda s:float(max(np.sqrt((s.astype(float)**2).sum())/len(s),1.0))})
        scaler=StandardScaler().fit(training[features].to_numpy(float))
        gp=GPSurrogate(seed,2).fit(scaler.transform(training[features].to_numpy(float)),training.tg_C.to_numpy(float),training.tg_se_C.to_numpy(float))
        posterior=gp.posterior(scaler.transform(usable[features].to_numpy(float)))
        prob=np.exp(log_prob_in_window(posterior,tg_min,tg_max))
        # Posterior interval satisfaction as in original notebook. Batch diversification
        # here uses normalized distance, not Nannan's fantasy conditioning.
        method='GP_target_interval_probability'
        rank=np.argsort(-np.nan_to_num(prob,nan=-1.))
        used=[]
        scale=np.ptp(usable[features].to_numpy(float),axis=0);scale=np.where(scale>0,scale,1.)
        for j in rank:
            if not used or np.min(np.linalg.norm((usable.iloc[used][features].to_numpy(float)-usable.iloc[j][features].to_numpy(float))/scale,axis=1))>=.12:
                used.append(int(j))
            if len(used)>=batch_size:break
        selected=usable.iloc[used].copy()
        selected['predicted_tg_C']=posterior.mean[used]
        selected['tg_predictive_sd_C']=posterior.predictive_sd[used]
        selected['p_tg_in_window']=prob[used]
        meta['gp_diagnostics']=gp.diagnostics()
    else:
        indexes=farthest_first(usable,obs,features,batch_size)
        selected=usable.iloc[indexes].copy()
        selected['predicted_tg_C']=np.nan
        selected['tg_predictive_sd_C']=np.nan
        selected['p_tg_in_window']=np.nan
        meta['reason']=f'Only {n_unique} distinct approved training points; require {min_gp} for exploratory GP fit'
    selected['selection_method']=method
    selected['proposed_job_type']=np.where(selected.structure_status.eq('CURED_NETWORK_VALIDATED'),'existing_network_MD','generate_and_cure_network')
    selected['tg_target_min_C']=tg_min
    selected['tg_target_max_C']=tg_max
    selected['proposal_status']='proposed_not_executed'
    meta['status']='OK';meta['method']=method;meta['n_selected']=len(selected)
    return selected.reset_index(drop=True),meta


def run(args):
    if args.mode=='composition' and args.generate_pool:
        candidates=build_composition_pool(args.n_dgeba,args.dds_counts,args.conversions,args.ratio_min,args.ratio_max)
        if args.write_pool:
            Path(args.write_pool).parent.mkdir(parents=True,exist_ok=True)
            candidates.to_csv(args.write_pool,index=False)
    else:
        if not args.candidates:raise ValueError('--candidates is required unless --mode composition --generate-pool')
        candidates=read_csv(args.candidates,REQUIRED_CAND)
    obs=read_csv(args.observations,REQUIRED_OBS)
    result,metadata=propose(candidates,obs,args.mode,args.tg_min,args.tg_max,args.batch_size,args.seed,args.allow_synthetic_demo,args.min_gp)
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True);result.to_csv(out,index=False)
    out.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2,default=str))
    print(json.dumps(metadata,indent=2,default=str))
    if len(result):print(result[['candidate_id','ratio','conversion_pct','proposed_job_type','selection_method']].to_string(index=False))
    return metadata


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['conversion','composition'],required=True)
    p.add_argument('--candidates')
    p.add_argument('--observations',default=str(HERE/'data/md_tg_observations.csv'))
    p.add_argument('--output',default=str(HERE/'outputs/proposals.csv'))
    p.add_argument('--tg-min',type=float,default=150.)
    p.add_argument('--tg-max',type=float,default=220.)
    p.add_argument('--batch-size',type=int,default=3)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--min-gp',type=int,default=4)
    p.add_argument('--allow-synthetic-demo',action='store_true')
    p.add_argument('--generate-pool',action='store_true')
    p.add_argument('--write-pool')
    p.add_argument('--n-dgeba',type=int,default=100)
    p.add_argument('--dds-counts',nargs='+',type=int,default=[40,45,50,55,60])
    p.add_argument('--conversions',nargs='+',type=int,default=[60,75,85,90])
    p.add_argument('--ratio-min',type=float,default=.8)
    p.add_argument('--ratio-max',type=float,default=1.2)
    return p

if __name__=='__main__':run(parser().parse_args())
