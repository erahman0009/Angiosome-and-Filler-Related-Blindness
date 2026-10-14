"""Decision-rule sensitivity and analytical verification for the CFD model."""
import json
import numpy as np
import pandas as pd
from scipy.stats import beta

OA_THRESHOLDS=[.50,.70,.85,.95,1.00]
CRA_SCALES=[.50,.75,1.,1.25,1.50]

def threshold_audit(inputs,scenarios,folder):
    rows=[]
    for name,g in scenarios.groupby('Scenario',sort=False):
        tone=inputs.set_index('profile_id').loc[g.profile_id,'tone'].to_numpy()
        for threshold in OA_THRESHOLDS:
            for scale in CRA_SCALES:
                oa=(g.pressure_surplus_mmHg.to_numpy()>0)&(g.transported_ml.to_numpy()>=threshold*g.path_volume_ml.to_numpy())
                cra=oa&(g.transported_ml.to_numpy()>=g.path_volume_ml.to_numpy()+scale*(.0014+.0028*tone))
                rows.append({'Scenario':name,'OA threshold':threshold,'CRA volume scale':scale,'n':len(g),'OA entry':float(oa.mean()),'CRA exposure':float(cra.mean())})
    full=pd.DataFrame(rows);full.to_csv(folder/'threshold_sensitivity_all_scenarios.csv',index=False)
    baseline=full[full.Scenario.str.startswith('S0')].copy()
    baseline.to_csv(folder/'threshold_sensitivity_baseline_full_grid.csv',index=False)
    baseline[baseline['CRA volume scale']==1].to_csv(folder/'Table_S10_Panel_A.csv',index=False)
    baseline[baseline['OA threshold']==.85].to_csv(folder/'Table_S10_Panel_B.csv',index=False)
    baseline.to_csv(folder/'threshold_sensitivity.csv',index=False)
    reference=full[(full['OA threshold']==.85)&(full['CRA volume scale']==1)]
    order=[]
    for (th,scale),g in full.groupby(['OA threshold','CRA volume scale'],sort=False):
        row={'OA threshold':th,'CRA volume scale':scale}
        for outcome in ['OA entry','CRA exposure']:
            ranking=g.sort_values(outcome,ascending=False).Scenario.str[:2].tolist()
            expected=reference.sort_values(outcome,ascending=False).Scenario.str[:2].tolist()
            row[outcome+' ordering']=' > '.join(ranking)
            row[outcome+' preserved']=ranking==expected
        order.append(row)
    pd.DataFrame(order).to_csv(folder/'threshold_scenario_ordering.csv',index=False)
    quantiles=[]
    for p in [.25,.5,.75,.90,.95]:
        empirical=float(inputs.choke_index.quantile(p))
        quantiles.append({'Percentile':p,'Theoretical index':float(beta.ppf(p,1.8,3.8)),'Sample index':empirical,'Modifier c=0.35':.35*empirical,'Modifier c=0.90':.9*empirical})
    pd.DataFrame(quantiles).to_csv(folder/'Table_S11_Choke_Quantiles.csv',index=False)
    audit={'threshold_combinations':25,'scenario_states':8,'profile_setting_evaluations':2400000,
           'OA_order_preserved':bool(all(x['OA entry preserved'] for x in order)),
           'CRA_order_preserved':bool(all(x['CRA exposure preserved'] for x in order)),
           'baseline_OA_range':[float(baseline['OA entry'].min()),float(baseline['OA entry'].max())],
           'baseline_CRA_range':[float(baseline['CRA exposure'].min()),float(baseline['CRA exposure'].max())],
           'reference_OA':float(reference.iloc[0]['OA entry']),'reference_CRA':float(reference.iloc[0]['CRA exposure']),
           'hydraulics_recomputed_for_thresholds':False,'classification_only':True}
    (folder/'threshold_audit.json').write_text(json.dumps(audit,indent=2))
    return audit

def analytical_audit(inputs,solve,area,mmhg,folder):
    # eta0=eta_inf makes the implemented Carreau filler exactly Newtonian.
    x=inputs.iloc[:1000].copy();x['eta0_Pa_s']=.018
    result=solve(x,c=0,site='upstream')
    r=x.diameter_mm.to_numpy()*.001/2
    absent=x.true_anastomosis.to_numpy()<.5
    r=np.where(absent,np.minimum(r,.105e-3),r)
    length=x.length_mm.to_numpy()*.001*np.where(absent,1.35,1.)
    angle=1+.85*(1-np.cos(np.deg2rad(x.angle_deg.to_numpy())))
    alpha=np.clip(x.intravascular_fraction.to_numpy(),.05,.95)
    mu=alpha*.018+(1-alpha)*.0035
    dp=np.maximum(x.pressure_coupling.to_numpy()*x.force_N.to_numpy()/area-(8+10*x.tone.to_numpy())*mmhg,0)
    exact=dp*np.pi*r**4/(8*mu*length*angle)
    error=float(np.max(np.abs(result.flow_ml_s.to_numpy()*1e-6-exact)/np.maximum(exact,1e-30)))
    if error>1e-12:raise RuntimeError('Newtonian Poiseuille limit failed')
    report={'Newtonian_limit_profiles':len(x),'maximum_relative_flow_error':error,'analytical_limit_passed':True,'external_validation':False}
    (folder/'analytical_verification.json').write_text(json.dumps(report,indent=2))
    return report
