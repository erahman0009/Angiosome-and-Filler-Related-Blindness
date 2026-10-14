"""Reproducible reduced-order hydraulic transport calculation, seed 42.
Run: python cfd_reproducibility_script.py --output ../
Dependencies: numpy, pandas, scipy, matplotlib, scikit-learn.
This uncalibrated model is conditional on intravascular exposure.
"""
import argparse,json,platform,sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import brentq
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split
import sklearn
from decision_audits import threshold_audit, analytical_audit
from publication_figures import plots

SEED=42; N=12000; MMHG=133.322; AREA=np.pi*(.0047/2)**2
TOL=1e-10; MAX_ITER=500
ROUTES=[('Glabella supratrochlear',.92,17.5,24,.64),('Nasal dorsum dorsal nasal',.86,22,31,.58),('Nasolabial angular',.74,31,44,.46),('Forehead supraorbital',.80,27,36,.52)]

def sample():
    rng=np.random.RandomState(SEED); rows=[]
    for i in range(N):
        route,d,l,a,p=ROUTES[i%4]; sbp=rng.normal(124,13); dbp=rng.normal(76,8)
        rows.append(dict(profile_id=i,route=route,diameter_mm=np.clip(rng.normal(d,.19),.38,1.55),length_mm=np.clip(rng.normal(l,4.8),9,43),angle_deg=np.clip(rng.normal(a,11),5,82),force_N=rng.uniform(6,30),volume_ml=rng.uniform(.015,.18),duration_s=rng.uniform(1.5,8),eta0_Pa_s=rng.uniform(75,190),systolic_mmHg=sbp,diastolic_mmHg=dbp,tone=rng.beta(2.2,2.5),true_anastomosis=rng.binomial(1,p),choke_index=rng.beta(1.8,3.8),pressure_coupling=rng.uniform(.025,.14),intravascular_fraction=rng.uniform(.22,.88)))
    return pd.DataFrame(rows)

def viscosity(eta0,gamma,alpha):
    filler=.018+(eta0-.018)*(1+(.22*np.maximum(gamma,1e-3))**2)**(-.38)
    return alpha*filler+(1-alpha)*.0035

def solve(data,overrides=None,c=.35,site='upstream',tol=TOL):
    x=data.copy()
    for k,v in (overrides or {}).items():x[k]=v
    x=x.drop(columns=['route','profile_id']); x={k:x[k].to_numpy(float) for k in x}
    map_pa=(x['diastolic_mmHg']+(x['systolic_mmHg']-x['diastolic_mmHg'])/3)*MMHG
    inlet=map_pa+x['pressure_coupling']*x['force_N']/AREA
    back=map_pa+(8+10*x['tone'])*MMHG; dp=inlet-back
    r=np.maximum(x['diameter_mm']*.001/2,5e-5); l=x['length_mm']*.001
    absent=x['true_anastomosis']<.5; r=np.where(absent,np.minimum(r,.105e-3),r); l=np.where(absent,l*1.35,l)
    if site=='upstream':r=r*np.maximum(.05,1-c*x['choke_index'])
    k=1+.85*(1-np.cos(np.deg2rad(x['angle_deg']))); b=8*l*k/(np.pi*r**4)
    alpha=np.clip(x['intravascular_fraction'],.05,.95); mu=x['eta0_Pa_s'].copy(); q=np.zeros(len(r))
    done=np.zeros(len(r),bool); iterations=np.zeros(len(r),int)
    for it in range(1,MAX_ITER+1):
        qnew=np.maximum(dp,0)/(b*mu); munew=viscosity(x['eta0_Pa_s'],4*qnew/(np.pi*r**3),alpha)
        dq=np.abs(qnew-q)/np.maximum(np.abs(qnew),1e-15); dm=np.abs(munew-mu)/np.maximum(np.abs(munew),1e-15)
        converged=(dq<tol)&(dm<tol); new=converged&(~done); iterations[new]=it; done=converged; q=qnew; mu=munew
        if done.all():break
    if not done.all():raise RuntimeError(f'Nonconverged profiles: {(~done).sum()}')
    # Recompute Q with returned viscosity and report closure residual explicitly.
    q=np.maximum(dp,0)/(b*mu); mucheck=viscosity(x['eta0_Pa_s'],4*q/(np.pi*r**3),alpha)
    residual=np.abs(mucheck-mu)/np.maximum(np.abs(mu),1e-15)
    if np.max(residual)>tol:raise RuntimeError('Closure residual exceeds tolerance')
    dead=np.pi*r**2*l; available=x['volume_ml']*1e-6*x['intravascular_fraction']; transport=np.minimum(q*x['duration_s'],available)
    terminal=(.0014+.0028*x['tone'])*1e-6
    if site=='downstream':terminal*=1+3*c*x['choke_index']
    oa=(dp>0)&(transport>=.85*dead); cra=oa&(transport>=dead+terminal)
    transit=dead/np.maximum(q,1e-15); reynolds=2*1060*q/(np.pi*r*mu)
    return pd.DataFrame(dict(inlet_pressure_mmHg=inlet/MMHG,ophthalmic_pressure_mmHg=back/MMHG,pressure_surplus_mmHg=dp/MMHG,mu_eff_Pa_s=mu,hydraulic_resistance_Pa_s_m3=b*mu,flow_ml_s=q/1e-6,path_volume_ml=dead/1e-6,transported_ml=transport/1e-6,fill_fraction=transport/dead,oa_entry=oa.astype(int),cra_exposure=cra.astype(int),cra_fraction=np.clip((transport-dead)/terminal,0,1),transit_uncapped_s=transit,transit_s=np.minimum(transit,120),transit_capped=(transit>=120).astype(int),reynolds=reynolds,wall_shear_pa=4*mu*q/(np.pi*r**3),effective_radius_m=r,effective_length_m=l,angle_multiplier=k,solver_iterations=iterations,solver_residual=residual))

def summarize(df):
    return dict(n=len(df),OA_entry=df.oa_entry.mean(),CRA_exposure=df.cra_exposure.mean(),pressure_surplus_mmHg=df.pressure_surplus_mmHg.mean(),flow_ml_s=df.flow_ml_s.mean(),transit_s=df.transit_s.median(),transported_ml=df.transported_ml.mean())

def run(out):
    data=out/'data'; fig=out/'figures'; data.mkdir(parents=True,exist_ok=True);fig.mkdir(exist_ok=True)
    inputs=sample(); base=pd.concat([inputs,solve(inputs)],axis=1); base.to_csv(data/'baseline_simulation_outputs.csv',index=False)
    q75=inputs.choke_index.quantile(.75)
    specs=[('S0 Baseline',{}),('S1 Low-pressure boundary',{'force_N':inputs.force_N.quantile(.25)}),('S2 High-pressure boundary',{'force_N':inputs.force_N.quantile(.75)}),('S3 High-viscosity filler',{'eta0_Pa_s':inputs.eta0_Pa_s.quantile(.75)}),('S4 No true anastomosis',{'true_anastomosis':0}),('S5 Forced arterial continuity',{'true_anastomosis':1}),('S6 High upstream choke resistance',{'choke_index':q75}),('S7 Route-specific high-risk extremes',{'force_N':inputs.force_N.quantile(.9),'eta0_Pa_s':inputs.eta0_Pa_s.quantile(.1),'true_anastomosis':1,'choke_index':0})]
    frames=[]; summaries=[]
    for name,overrides in specs:
        s=solve(inputs,overrides); s.insert(0,'Scenario',name); s.insert(1,'profile_id',inputs.profile_id);s.insert(2,'route',inputs.route);frames.append(s)
        summaries.append({'Scenario':name,**summarize(s)})
    scenarios=pd.concat(frames,ignore_index=True); scenarios.to_csv(data/'scenario_evaluations.csv',index=False)
    stress=pd.DataFrame(summaries); stress['Change OA']=stress.OA_entry-stress.OA_entry.iloc[0];stress['Change CRA']=stress.CRA_exposure-stress.CRA_exposure.iloc[0];stress.to_csv(data/'stress_test_outputs.csv',index=False)
    routes=pd.DataFrame([{'Route':r,**summarize(g)} for r,g in base.groupby('route')]);routes.to_csv(data/'route_summary.csv',index=False)
    bands=pd.cut(base.cra_fraction,[-.001,.05,.35,.75,np.inf],labels=['Very low','Low','Moderate','High'])
    band=pd.DataFrame([{'CRA band':str(b),'n':len(g),'Pressure surplus (mmHg)':g.pressure_surplus_mmHg.mean(),'Diameter (mm)':g.diameter_mm.mean(),'Effective viscosity (Pa s)':g.mu_eff_Pa_s.mean(),'Flow (mL/s)':g.flow_ml_s.mean(),'Transported volume (mL)':g.transported_ml.mean()} for b,g in base.assign(band=bands).groupby('band',observed=False)])
    band.to_csv(data/'baseline_CRA_bands.csv',index=False)
    grid=[]; gridframes=[]
    for site in ['upstream','downstream']:
        ref=solve(inputs,{'choke_index':q75},0,site)
        for c in [0,.15,.35,.6,.9]:
            s=solve(inputs,{'choke_index':q75},c,site);grid.append({'Action site':site,'Coefficient c':c,'Index x':q75,'Modifier c x':c*q75,'OA entry':s.oa_entry.mean(),'CRA exposure':s.cra_exposure.mean(),'Change OA':s.oa_entry.mean()-ref.oa_entry.mean(),'Change CRA':s.cra_exposure.mean()-ref.cra_exposure.mean()})
            s.insert(0,'profile_id',inputs.profile_id);s.insert(1,'Action site',site);s.insert(2,'Coefficient c',c);gridframes.append(s)
    grid=pd.DataFrame(grid);grid.to_csv(data/'choke_sensitivity_summary.csv',index=False);pd.concat(gridframes,ignore_index=True).to_csv(data/'choke_sensitivity_full_grid.csv',index=False)
    decision_report=threshold_audit(inputs,scenarios,data)
    analytical_report=analytical_audit(inputs,solve,AREA,MMHG,data)
    s7=scenarios[scenarios.Scenario.str.startswith('S7')];pd.DataFrame([{'Route':r,**summarize(g)} for r,g in s7.groupby('route')]).to_csv(data/'S7_route_outputs.csv',index=False)
    # Held-out surrogate screening; this is not Sobol analysis or clinical calibration.
    feats=['diameter_mm','length_mm','angle_deg','force_N','volume_ml','duration_s','eta0_Pa_s','tone','true_anastomosis','choke_index','pressure_coupling','intravascular_fraction']; sens=[]; scores={}
    for outcome in ['oa_entry','cra_exposure']:
        train,test=train_test_split(np.arange(N),test_size=.25,random_state=42,stratify=base[outcome])
        model=RandomForestClassifier(n_estimators=200,min_samples_leaf=3,class_weight='balanced',random_state=42,n_jobs=1)
        model.fit(base.iloc[train][feats],base.iloc[train][outcome]);scores[outcome]=model.score(base.iloc[test][feats],base.iloc[test][outcome])
        importance=permutation_importance(model,base.iloc[test][feats],base.iloc[test][outcome],n_repeats=10,random_state=42,n_jobs=1,scoring='roc_auc')
        for j,f in enumerate(feats):sens.append({'Outcome':outcome,'Parameter':f,'Permutation AUROC decrease':importance.importances_mean[j],'SD':importance.importances_std[j],'Spearman rho':spearmanr(base[f],base[outcome]).statistic})
    pd.DataFrame(sens).to_csv(data/'sensitivity_screening.csv',index=False)
    # Numerical verification is separate from external validation.
    alternative=solve(inputs,tol=1e-8); rel=abs(alternative.flow_ml_s-base.flow_ml_s)/np.maximum(base.flow_ml_s,1e-15)
    ids=np.random.RandomState(43).choice(N,200,replace=False); errors=[]
    for i in ids:
        row=base.iloc[i];dp=row.pressure_surplus_mmHg*MMHG;r=row.effective_radius_m;l=row.effective_length_m;k=row.angle_multiplier;a=np.clip(row.intravascular_fraction,.05,.95);b=8*l*k/(np.pi*r**4)
        low=max(dp,0)/(b*row.eta0_Pa_s);high=max(dp,0)/(b*(a*.018+(1-a)*.0035))
        root=brentq(lambda q:q*b*viscosity(row.eta0_Pa_s,4*q/(np.pi*r**3),a)-max(dp,0),low,high,xtol=1e-30,rtol=1e-13)
        errors.append(abs(row.flow_ml_s*1e-6-root)/max(root,1e-30))
    zero=solve(inputs,{'pressure_coupling':0})
    audit={'seed':SEED,'baseline_profiles':N,'scenario_evaluations':len(scenarios),'choke_profile_evaluations':120000,'threshold_decision_evaluations':2400000,'tolerance':TOL,'iteration_limit':MAX_ITER,'max_iterations_all_hydraulic_evaluations':int(max(scenarios.solver_iterations.max(),max(f.solver_iterations.max() for f in gridframes))),'max_residual_all_hydraulic_evaluations':float(max(scenarios.solver_residual.max(),max(f.solver_residual.max() for f in gridframes))),'baseline_transport_max_error_ml':float(abs(base.transported_ml-np.minimum(base.flow_ml_s*base.duration_s,base.volume_ml*base.intravascular_fraction)).max()),'tolerance_comparison_max_relative_flow_error':float(rel.max()),'tolerance_comparison_OA_changes':int((alternative.oa_entry!=base.oa_entry).sum()),'tolerance_comparison_CRA_changes':int((alternative.cra_exposure!=base.cra_exposure).sum()),'independent_brent_profiles':200,'max_brent_relative_flow_error':float(max(errors)),'zero_pressure_surplus_transport_nonzero':int((zero.transported_ml!=0).sum()),'baseline_re_median':float(base.reynolds.median()),'baseline_re_p95':float(base.reynolds.quantile(.95)),'baseline_re_p99':float(base.reynolds.quantile(.99)),'baseline_re_above2000_n':int((base.reynolds>2000).sum()),'baseline_transit_capped_n':int(base.transit_capped.sum()),'surrogate_test_accuracy':scores,'external_validation_performed':False,'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'scikit-learn':sklearn.__version__}}
    (data/'numerical_audit.json').write_text(json.dumps(audit,indent=2));plots(base,stress,routes,grid,fig)
    exclusion=[]
    for state,g in scenarios.groupby('Scenario',sort=False):
        clean=g[g.reynolds<=2000]
        exclusion.append({'Scenario':state,'n':len(g),'excluded_Re_above_2000':len(g)-len(clean),'OA_all':g.oa_entry.mean(),'OA_laminar_subset':clean.oa_entry.mean(),'CRA_all':g.cra_exposure.mean(),'CRA_laminar_subset':clean.cra_exposure.mean()})
    pd.DataFrame(exclusion).to_csv(data/'laminar_exclusion_audit.csv',index=False)
    export_checks=[]
    for tag,frame in [('scenario',scenarios),('choke',pd.concat(gridframes,ignore_index=True))]:
        grouped=frame.groupby('Scenario',sort=False) if tag=='scenario' else frame.groupby(['Action site','Coefficient c'],sort=False)
        for condition,g in grouped:
            bb=base.iloc[g.profile_id.to_numpy()]
            vt=np.minimum(g.flow_ml_s.to_numpy()*bb.duration_s.to_numpy(),bb.volume_ml.to_numpy()*bb.intravascular_fraction.to_numpy())
            oa=((g.pressure_surplus_mmHg>0)&(g.transported_ml>=.85*g.path_volume_ml)).astype(int)
            terminal=.0014+.0028*bb.tone.to_numpy()
            if tag=='choke' and condition[0]=='downstream':terminal*=1+3*condition[1]*q75
            cra=((oa>0)&(g.transported_ml>=g.path_volume_ml+terminal)).astype(int)
            export_checks.append({'population':tag,'condition':str(condition),'volume_error_ml':float(abs(vt-g.transported_ml).max()),'OA_mismatches':int((oa!=g.oa_entry).sum()),'CRA_mismatches':int((cra!=g.cra_exposure).sum())})
    if any(x['volume_error_ml']>1e-12 or x['OA_mismatches'] or x['CRA_mismatches'] for x in export_checks):
        raise RuntimeError('Exported transport outcomes are inconsistent')
    pd.DataFrame(export_checks).to_csv(data/'export_consistency_audit.csv',index=False)
    print(stress[['Scenario','OA_entry','CRA_exposure','Change OA']].to_string(index=False));print(json.dumps(audit,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]);args=p.parse_args();run(args.output.resolve())
