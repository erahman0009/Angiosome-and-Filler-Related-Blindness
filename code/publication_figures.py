"""Publication figures drawn exclusively from converged simulation outputs."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BLUE='#176585';ORANGE='#CE7135';GREEN='#438B73';PURPLE='#806AA3'
ROUTE_COLORS=[BLUE,ORANGE,GREEN,PURPLE]

def plots(base,stress,routes,grid,folder):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
        'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,
        'savefig.dpi':350,'pdf.fonttype':42})
    def save(f,name):
        for suffix in ['png','pdf','tiff']:
            tmp=folder/('.'+name+'.'+suffix)
            kwargs={'bbox_inches':'tight'}
            if suffix!='pdf':kwargs['dpi']=350
            if suffix=='tiff':kwargs['pil_kwargs']={'compression':'tiff_lzw'}
            f.savefig(tmp,**kwargs)
            tmp.replace(folder/(name+'.'+suffix))
        plt.close(f)
    def mark(ax,label,title):ax.set_title(label+'  '+title,loc='left',fontweight='bold',pad=10)
    def ecdf(ax,values,label,color,log=False):
        v=np.sort(np.asarray(values));v=v[v>0] if log else v
        ax.plot(v,np.arange(1,len(v)+1)/len(v),label=label,color=color,lw=1.8)
        if log:ax.set_xscale('log')
    from matplotlib.patches import FancyBboxPatch
    fig,axes=plt.subplots(3,1,figsize=(11,6.5),layout='constrained')
    for a in axes:a.set_xlim(0,11);a.set_ylim(0,1);a.axis('off')
    def box(a,x,w,text,color=BLUE):
        a.add_patch(FancyBboxPatch((x,.18),w,.57,boxstyle='round,pad=0.02',facecolor='white',edgecolor=color,lw=1.5))
        a.text(x+w/2,.465,text,ha='center',va='center',fontsize=9)
    box(axes[0],.1,2.1,'Injection pressure\nMAP + coupling × F / area')
    box(axes[0],3.0,4.5,'Effective facial-to-ophthalmic conduit\nR = 8 μ L K / (π r⁴); Q = max(ΔP, 0) / R')
    box(axes[0],8.2,2.6,'Ophthalmic boundary\nMAP + tone offset')
    for x,y in [(2.25,3),(7.55,8.2)]:axes[0].annotate('',(y,.46),(x,.46),arrowprops={'arrowstyle':'->','color':BLUE,'lw':1.6})
    axes[0].text(.1,.97,'A  Reduced-order hydraulic calculation',fontweight='bold',fontsize=11)
    box(axes[1],.1,3.4,'Transported volume\nmin(Q × duration, available bolus)')
    box(axes[1],4.2,2.6,'OA entry\nPositive surplus; ≥ h × route volume')
    box(axes[1],7.5,3.3,'CRA exposure\nOA entry; ≥ route + s × terminal volume')
    for x,y in [(3.55,4.2),(6.85,7.5)]:axes[1].annotate('',(y,.46),(x,.46),arrowprops={'arrowstyle':'->','color':BLUE,'lw':1.6})
    axes[1].text(.1,.97,'B  Volume-based conditional outcomes',fontweight='bold',fontsize=11)
    box(axes[2],.1,5.1,'Upstream implementation\nr → r(1 − c × index); modifies resistance and route volume')
    box(axes[2],5.7,5.1,'Downstream implementation\nTerminal volume × (1 + 3c × index); upstream flow unchanged',ORANGE)
    axes[2].text(.1,.97,'C  Distinct choke hypotheses',fontweight='bold',fontsize=11)
    save(fig,'Figure_1_model_architecture')
    grouped=list(base.groupby('route',sort=True))
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    maprows=[]
    xp=np.linspace(0,2000,17);yp=np.linspace(.38,1.55,14)
    for j,(ax,(route,g)) in enumerate(zip(axes.flat,grouped)):
        n,_,_=np.histogram2d(g.pressure_surplus_mmHg,g.diameter_mm,bins=[xp,yp])
        hits,_,_=np.histogram2d(g.pressure_surplus_mmHg,g.diameter_mm,bins=[xp,yp],weights=g.oa_entry)
        z=np.divide(hits,n,out=np.full_like(hits,np.nan),where=n>=5)
        im=ax.pcolormesh(xp,yp,z.T,vmin=0,vmax=1,cmap='viridis',rasterized=True)
        mark(ax,'ABCD'[j],route);ax.set(xlabel='Pressure surplus (mmHg)',ylabel='Sampled route diameter (mm)')
        ax.text(.02,.97,'n = 3,000; bins with <5 profiles masked',transform=ax.transAxes,va='top',fontsize=8,
            bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
        for a in range(n.shape[0]):
            for b in range(n.shape[1]):maprows.append([route,xp[a],xp[a+1],yp[b],yp[b+1],int(n[a,b]),int(hits[a,b]),z[a,b]])
    fig.colorbar(im,ax=list(axes.flat),label='Conditional OA transport proportion',shrink=.9)
    pd.DataFrame(maprows,columns=['Route','Pressure lower mmHg','Pressure upper mmHg','Diameter lower mm','Diameter upper mm','n','OA events','Proportion']).to_csv(folder.parent/'data'/'Figure_2_bin_data.csv',index=False)
    save(fig,'Figure_2_pressure_geometry')

    fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    names=[r.replace(' supratrochlear','').replace(' dorsal nasal','').replace(' angular','').replace(' supraorbital','') for r,g in grouped]
    y=np.arange(4)
    for delta,field,color,label in [(-.16,'oa_entry',BLUE,'OA entry'),(.16,'cra_exposure',ORANGE,'CRA exposure')]:
        props=np.array([g[field].mean() for r,g in grouped]);n=3000;z=1.96
        mid=(props+z*z/(2*n))/(1+z*z/n);half=z*np.sqrt(props*(1-props)/n+z*z/(4*n*n))/(1+z*z/n)
        ax[0,0].errorbar(props,y+delta,xerr=[props-(mid-half),(mid+half)-props],fmt='o',color=color,label=label,capsize=3)
    ax[0,0].set(yticks=y,yticklabels=names,xlabel='Conditional proportion',xlim=(0,.55));ax[0,0].legend(frameon=False)
    mark(ax[0,0],'A','Route-level conditional transport')
    for (r,g),color,name in zip(grouped,ROUTE_COLORS,names):ecdf(ax[0,1],g.transported_ml/(g.volume_ml*g.intravascular_fraction),name,color)
    ax[0,1].set(xlabel='Fraction of available intravascular bolus transported',ylabel='Cumulative proportion',xlim=(0,1));mark(ax[0,1],'B','Delivered bolus fraction')
    for (r,g),color,name in zip(grouped,ROUTE_COLORS,names):ecdf(ax[1,0],g.flow_ml_s,name,color,True)
    ax[1,0].set(xlabel='Flow (mL/s)',ylabel='Cumulative proportion');ax[1,0].legend(frameon=False,fontsize=8);mark(ax[1,0],'C','Route flow distributions')
    h=ax[1,1].hexbin(np.log10(base.flow_ml_s),base.cra_fraction,gridsize=40,mincnt=1,bins='log',cmap='viridis')
    ax[1,1].set(xlabel='log₁₀ flow (mL/s)',ylabel='Terminal CRA filling fraction');fig.colorbar(h,ax=ax[1,1],label='Profiles (log colour scale)');mark(ax[1,1],'D','Flow and terminal filling')
    save(fig,'Figure_3_route_outputs')

    fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    shear=4*base.flow_ml_s*1e-6/(np.pi*base.effective_radius_m**3)
    h=ax[0,0].hexbin(np.log10(shear),np.log10(base.mu_eff_Pa_s),gridsize=40,mincnt=1,bins='log',cmap='viridis')
    ax[0,0].set(xlabel='log₁₀ apparent shear rate (s⁻¹)',ylabel='log₁₀ effective viscosity (Pa·s)');fig.colorbar(h,ax=ax[0,0],label='Profiles');mark(ax[0,0],'A','Converged rheological states')
    x=np.log10(base.path_volume_ml);z=np.log10(base.transported_ml)
    h=ax[0,1].hexbin(x,z,gridsize=45,mincnt=1,bins='log',cmap='viridis');lim=[float(x.min()),float(x.max())]
    ax[0,1].plot(lim,np.array(lim)+np.log10(.85),'--',color=ORANGE,lw=1.6,label='85% OA filling criterion')
    ax[0,1].plot(lim,lim,':',color='black',lw=1.2,label='Full route dead volume');ax[0,1].legend(frameon=False,fontsize=8)
    ax[0,1].set(xlabel='log₁₀ route dead volume (mL)',ylabel='log₁₀ transported volume (mL)');fig.colorbar(h,ax=ax[0,1],label='Profiles');mark(ax[0,1],'B','Transport and route volume')
    for field,color,label in [('transit_uncapped_s',BLUE,'Uncapped hydraulic transit'),('transit_s',ORANGE,'Reported transit, capped at 120 s')]:ecdf(ax[1,0],base[field],label,color,True)
    ax[1,0].axvline(120,color='black',ls=':',lw=1);ax[1,0].set(xlabel='Transit time (s)',ylabel='Cumulative proportion');ax[1,0].legend(frameon=False,fontsize=8);mark(ax[1,0],'C','Transit reporting and censoring')
    ecdf(ax[1,1],base.reynolds,'Baseline',BLUE,True);ax[1,1].axvline(2000,color=ORANGE,ls='--',lw=1.3,label='Re = 2,000 audit cutoff')
    ax[1,1].set(xlabel='Reynolds number',ylabel='Cumulative proportion');ax[1,1].legend(frameon=False,fontsize=8);mark(ax[1,1],'D','Laminar-flow audit')
    save(fig,'Figure_4_transport_physics')

    scenario=pd.read_csv(folder.parent/'data'/'scenario_evaluations.csv');reference=scenario[scenario.Scenario.str.startswith('S0')].sort_values('profile_id')
    fig,ax=plt.subplots(2,2,figsize=(11,9),layout='constrained');labels=stress.Scenario.str[:2];y=np.arange(8)
    vals=stress[['OA_entry','CRA_exposure']].to_numpy();im=ax[0,0].imshow(vals,vmin=0,vmax=1,cmap='viridis',aspect='auto')
    ax[0,0].set(xticks=[0,1],xticklabels=['OA entry','CRA exposure'],yticks=y,yticklabels=labels)
    for i in range(8):
        for j in range(2):ax[0,0].text(j,i,f'{vals[i,j]:.4f}',ha='center',va='center',color='white' if vals[i,j]<.6 else 'black')
    fig.colorbar(im,ax=ax[0,0],label='Conditional proportion');mark(ax[0,0],'A','Scenario transport proportions')
    for offset,outcome,col,color,label in [(-.15,'oa_entry','Change OA',BLUE,'OA entry'),(.15,'cra_exposure','Change CRA',ORANGE,'CRA exposure')]:
        means=[];ses=[]
        for name in stress.Scenario:
            g=scenario[scenario.Scenario==name].sort_values('profile_id');d=g[outcome].to_numpy()-reference[outcome].to_numpy();means.append(d.mean());ses.append(1.96*d.std(ddof=1)/np.sqrt(len(d)))
        ax[0,1].errorbar(means,y+offset,xerr=ses,fmt='o',color=color,label=label,capsize=2)
    ax[0,1].axvline(0,color='black',lw=.7);ax[0,1].set(yticks=y,yticklabels=labels,xlabel='Paired change from S0');ax[0,1].legend(frameon=False);mark(ax[0,1],'B','Counterfactual changes')
    for name,g in scenario.groupby('Scenario',sort=False):ecdf(ax[1,0],g.transported_ml,name[:2],plt.get_cmap('tab10')(int(name[1])))
    ax[1,0].set(xlabel='Transported volume (mL)',ylabel='Cumulative proportion');ax[1,0].legend(frameon=False,ncol=4,fontsize=8);mark(ax[1,0],'C','Scenario volume distributions')
    t=pd.read_csv(folder.parent/'data'/'threshold_sensitivity_all_scenarios.csv')
    for name,g in t.groupby('Scenario',sort=False):
        sub=g[g['CRA volume scale']==1].sort_values('OA threshold');ax[1,1].plot(sub['OA threshold'],sub['OA entry'],marker='o',ms=3,label=name[:2],color=plt.get_cmap('tab10')(int(name[1])))
    ax[1,1].set(xlabel='OA filling threshold',ylabel='Conditional OA proportion');ax[1,1].legend(frameon=False,ncol=4,fontsize=8);mark(ax[1,1],'D','Decision-threshold robustness')
    save(fig,'Figure_5_stress_tests')

    full=pd.read_csv(folder.parent/'data'/'choke_sensitivity_full_grid.csv')
    fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for site,color,style in [('upstream',BLUE,'-'),('downstream',ORANGE,'--')]:
        g=grid[grid['Action site']==site]
        for a,outcome in [(ax[0,0],'OA entry'),(ax[0,1],'CRA exposure')]:a.plot(g['Modifier c x'],g[outcome],style,color=color,marker='o',label=site)
    for a,label,title in [(ax[0,0],'A','OA entry'),(ax[0,1],'B','CRA exposure')]:
        a.axhline(base['oa_entry' if label=='A' else 'cra_exposure'].mean(),color=GREEN,ls=':',label='S0 sampled choke state')
        a.set(xlabel='Coefficient × fixed choke index',ylabel='Conditional proportion');a.legend(frameon=False,fontsize=8);mark(a,label,title)
    upstream=full[full['Action site']=='upstream'];ref=upstream[upstream['Coefficient c']==0].sort_values('profile_id')
    mods=[];ratios=[];volume=[];lost=[];gained=[]
    for c,g in upstream.groupby('Coefficient c',sort=True):
        g=g.sort_values('profile_id');mods.append(c*grid['Index x'].iloc[0]);ratios.append(np.median(g.flow_ml_s.to_numpy()/ref.flow_ml_s.to_numpy()));volume.append(np.median(g.path_volume_ml.to_numpy()/ref.path_volume_ml.to_numpy()))
        lost.append(np.mean((ref.oa_entry.to_numpy()==1)&(g.oa_entry.to_numpy()==0)));gained.append(np.mean((ref.oa_entry.to_numpy()==0)&(g.oa_entry.to_numpy()==1)))
    ax[1,0].plot(mods,ratios,'o-',color=BLUE,label='Median flow ratio');ax[1,0].plot(mods,volume,'s-',color=ORANGE,label='Route dead-volume ratio')
    ax[1,0].set(xlabel='Realised upstream radius reduction',ylabel='Ratio to zero-constriction reference');ax[1,0].legend(frameon=False,fontsize=8);mark(ax[1,0],'C','Flow-volume competition')
    ax[1,1].plot(mods,lost,'o-',color=BLUE,label='OA entry lost');ax[1,1].plot(mods,gained,'s-',color=ORANGE,label='OA entry gained')
    ax[1,1].set(xlabel='Realised upstream radius reduction',ylabel='Fraction of paired profiles');ax[1,1].legend(frameon=False,fontsize=8);mark(ax[1,1],'D','Paired OA-state changes')
    pd.DataFrame({'Radius modifier':mods,'Median flow ratio':ratios,'Dead-volume ratio':volume,'OA entry lost':lost,'OA entry gained':gained}).to_csv(folder.parent/'data'/'Figure_6_mechanism_data.csv',index=False)
    save(fig,'Figure_6_choke_sensitivity')
