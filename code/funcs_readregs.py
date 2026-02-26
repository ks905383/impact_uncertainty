# Being developed in read_regresults.ipynb

import xarray as xr
import numpy as np
import pandas as pd
import os
import re
import glob
import warnings
from tqdm.notebook import tqdm
from dask.diagnostics import ProgressBar

import distributed

from funcs_support import get_filepaths,get_params,utility_save
dir_list = get_params()


def coeffs_gdppc(obstreat='sep'):

    fns = pd.DataFrame({'path':np.sort(glob.glob(dir_list['regs']+'gdppc/*.csv'))})
    fns['proj_base'] = [re.split('post1980',re.split(r'\.csv',re.split('allyears',re.split(r'\_',re.split(r'\/',fn)[-1])[-1])[0])[0])[0] for fn in fns.path.values]
    fns['timeframe'] = [re.split(r'\.',re.split(row[1]['proj_base'],row[1]['path'])[1])[0] for row in fns.iterrows()]
    fns['datatype'] = [re.split(r'\_',re.split(r'\/estimated',fn)[-1])[0] for fn in fns.path.values]
    fns = fns.set_index(['proj_base','datatype','timeframe'])
    
    def get_coeffs_gdppc(output,obstreat,fns=fns):
        ''' Reading Burke, Hsiang, and Miguel regression output
    
        Parameters
        --------------
        output : 'coeffs' or 'vcv'
    
        obstreat : str, by default 'sep'
            'sep': regressions run on each reanalysis separately
            'joint': regressions run jointly interacted with reanalysis product
            
    
        Returns
        -------------
        coeffs or vcv
    
        '''
        if output not in ['coeffs','vcv']:
            raise KeyError(f'`output` must be "coeffs" or "vcv", but is "{str(output)}"')
        
        if obstreat not in ['sep','joint']:
            raise KeyError(f'`obstreat` must be "sep" or "joint", but is "{str(obstreat)}"')
        
        dss = []
        if obstreat == 'sep':
            
            for proj_base in [mod for mod in fns.index.unique(level='proj_base') if mod != 'joint']:
                for typ in fns.index.unique(level='timeframe'):
                    col_idxs = ['tassumpoly_1','tassumpoly_2']
    
                    if output == 'coeffs':
                        df = pd.read_csv(fns.loc[(proj_base,'Coefficients',typ)].path)
                        df.columns = col_idxs
                        ds_tmp = df.to_xarray().to_dataarray(dim='pred',name='coeffs').isel(index=0,drop=True)
                    elif output == 'vcv':
                        df = pd.read_csv(fns.loc[(proj_base,'VCV',typ)].path)
                        df.columns = col_idxs 
                        df.index = pd.Index(col_idxs)
                        ds_tmp = df.to_xarray().to_dataarray(dim='pred_1',name='vcv').rename({'index':'pred_2'})
                    
                    ds_tmp = ds_tmp.expand_dims({'proj_base':[proj_base],
                                                     'timeframe':[typ]})
                    dss.append(ds_tmp)
            dss = xr.combine_by_coords(dss)
        elif obstreat == 'joint':
            # Get order it's encoded in by STATA
            modlist = np.sort([k for k in fns.index.unique(level='proj_base') if k != 'joint'])
        
            for typ in fns.index.unique(level='timeframe'):
                if output == 'coeffs':
                    df = pd.read_csv(fns.loc[('joint','Coefficients',typ)].path)
                    df = df.iloc[[0],0:(len(modlist)*2)]
                    df.columns = pd.MultiIndex.from_product([['tassumpoly_1','tassumpoly_2'],modlist],
                                               names = ['pred','proj_base'])
                    ds_tmp = df.melt(value_name='coeffs').set_index(['pred','proj_base']).to_xarray()
                    
                elif output == 'vcv':
                    df = pd.read_csv(fns.loc[('joint','VCV',typ)].path)
                    df = df.iloc[0:(len(modlist)*2),0:(len(modlist)*2)]
                    
                    df.columns = pd.MultiIndex.from_product([['tassumpoly_1','tassumpoly_2'],modlist],
                                                           names = ['pred_1','proj_base_1'])
                    df.index = pd.MultiIndex.from_product([['tassumpoly_1','tassumpoly_2'],modlist],
                                                           names = ['pred_2','proj_base_2'])
                    
                    ds_tmp = (df.melt(ignore_index=False,value_name='vcv').reset_index().
                              set_index(['pred_1','proj_base_1','pred_2','proj_base_2']).to_xarray())
                    
                ds_tmp = ds_tmp.expand_dims({'timeframe':[typ]})
                dss.append(ds_tmp)
            dss = xr.combine_by_coords(dss)#.transpose('pred','proj_base','timeframe')
    
        return dss

    coeffs,vcv = get_coeffs_gdppc('coeffs',obstreat),get_coeffs_gdppc('vcv',obstreat)

    # Turn ERA5-025 into ERA5
    coeffs['proj_base'] = [re.sub(r'^ERA5-025$','ERA5',mod) for mod in coeffs.proj_base.values]
    if 'proj_base_1' in vcv.sizes:
        for k in [1,2]:
            vcv[f'proj_base_{str(k)}'] = [re.sub(r'^ERA5-025$','ERA5',mod) for mod in vcv[f'proj_base_{str(k)}'].values]
    else:
        vcv['proj_base'] = [re.sub(r'^ERA5-025$','ERA5',mod) for mod in vcv.proj_base.values]
    

    return coeffs,vcv


def coeffs_mortality(obstreat='sep',regtyp='pooled'):
    fns = pd.DataFrame({'path':np.sort(glob.glob(dir_list['regs']+'mortality/*.csv'))})
    fns['proj_base'] = [re.split(r'public\_',re.split(r'\.csv',re.split(r'\/',fn)[-1])[0])[-1] for fn in fns.path.values]
    fns['spectype'] = [re.split(r'\_response\_',re.split(r'\/',fn)[-1])[0] for fn in fns.path.values]
    fns['datatype'] = ['VCV' if '_VCV_' in fn else 'Coefficients' for fn in fns.path.values]
    fns = fns.set_index(['proj_base','datatype','spectype'])
    
    fns = fns.iloc[['nototal' not in fn for fn in fns.path.values],:]
    
    fns = fns.sort_index()

    def get_coeffs_mortality(output,obstreat,regtyp='pooled',fns=fns,drop_intersect = True):
        ''' Reading Carleton et al. mortality regression output
    
        Parameters
        --------------
        output : 'coeffs' or 'vcv'
    
        obstreat : str, by default 'sep'
            'sep': regressions run on each reanalysis separately
            'joint': regressions run jointly interacted with reanalysis product
    
        regtyp : str, by default 'pooled'
            'pooled': single regression for the globe
            'Agespec_interaction': regression interacted with GDPpc / long-run climate
    
        drop_intersect : bool, by default True
            If True, then intersect is removed from output
    
        Returns
        -------------
        coeffs or vcv
    
        '''
        if output not in ['coeffs','vcv']:
            raise KeyError(f'`output` must be "coeffs" or "vcv", but is "{str(output)}"')
        
        if obstreat not in ['sep','joint']:
            raise KeyError(f'`obstreat` must be "sep" or "joint", but is "{str(obstreat)}"')
    
        if regtyp not in ['pooled','Agespec_interaction']:
            raise KeyError(f'`regtyp` must be "pooled" or "Agespec_interaction", but is "{str(regtyp)}"')
            
        dss = []
        if obstreat == 'sep':
            for proj_base in [mod for mod in fns.index.unique(level='proj_base') if mod != 'joint']:
                if regtyp == 'pooled':
                    col_idxs = [*[f'tassumpoly_{k}' for k in np.arange(1,5)],'intersect']
    
                    if output == 'coeffs':
                        df = pd.read_csv(fns.loc[(proj_base,'Coefficients',regtyp)].path)
                        df.columns = col_idxs
                        ds_tmp = df.to_xarray().to_dataarray(dim='pred',name='coeffs').isel(index=0,drop=True)
                        ds_tmp = ds_tmp.expand_dims({'covar':['1'],'ageset':['total']})
    
                    elif output == 'vcv':
                        df = pd.read_csv(fns.loc[(proj_base,'VCV',regtyp)].path)
                        df.columns = col_idxs
                        df.index = col_idxs
                        ds_tmp = df.to_xarray().to_dataarray(dim='pred_1',name='vcv').rename({'index':'pred_2'})
    
                    ds_tmp = ds_tmp.expand_dims({'proj_base':[proj_base],
                                                  'spectype':[regtyp]}) 
                    ds_tmp = ds_tmp.to_dataset()
                
                    dss.append(ds_tmp)
                elif regtyp == 'Agespec_interaction':
                    # Set columns from fully interacted model reg call
                    #* Fully interacted model
                    # reghdfe deathrate_w99 c.tavg_poly_1#i.agegroup_ c.tavg_poly_2#i.agegroup_ ///
                    #				c.tavg_poly_3#i.agegroup_ c.tavg_poly_4#i.agegroup_ ///
                    #				c.tavg_poly_1#c.loggdppc_adm1_avg#i.agegroup_ c.tavg_poly_2#c.loggdppc_adm1_avg#i.agegroup_ ///
                    #				c.tavg_poly_3#c.loggdppc_adm1_avg#i.agegroup_ c.tavg_poly_4#c.loggdppc_adm1_avg#i.agegroup_ ///
                    #				c.tavg_poly_1#c.lr_tavg_adm1_avg#i.agegroup_ c.tavg_poly_2#c.lr_tavg_adm1_avg#i.agegroup_ ///
                    #				c.tavg_poly_3#c.lr_tavg_adm1_avg#i.agegroup_ c.tavg_poly_4#c.lr_tavg_adm1_avg#i.agegroup_ 
                    # In this order: "0-4"' `"5-64"' `"65+"' `"total"'
                    levels = ['young','middle','old','total']
                    tpreds = [f'tassumpoly_{k}' for k in np.arange(1,5)]
                    covars = ['1','tasmean','gdppc']
                    
                    tpreds_full = list(np.repeat(tpreds,len(levels)))*len(covars)
                    levels_full = (levels*len(tpreds))*len(covars)
                    covars_full = np.repeat(covars,len(tpreds)*len(levels))
                    
                    col_idxs = pd.DataFrame({'pred':tpreds_full,
                                              'covar':covars_full,
                                              'ageset':levels_full})
                    # Really should just instead repeat it x times across agesets
                    col_idxs = pd.concat([col_idxs,
                                          pd.DataFrame({'pred':['intersect'],
                                                       'covar':['1'],
                                                       'ageset':['intersect']})])
                    col_idxs = pd.MultiIndex.from_frame(col_idxs)
    
                    # Load coefficients
                    if output == 'coeffs':
                        df = pd.read_csv(fns.loc[(proj_base,'Coefficients',regtyp)].path)
                        df.columns = pd.MultiIndex.from_frame(col_idxs)
                        
                        # Prepare for xarray
                        df['proj_base'] = proj_base
                        df['spectype'] = regtyp
                        df = df.set_index(['proj_base','spectype'])
                        df = df.melt(ignore_index=False,value_name='coeffs').reset_index().set_index(['proj_base','spectype','pred','covar','ageset'])
                        
                        # To xarray
                        ds_tmp = df.to_xarray()
                    elif output == 'vcv':
                        df = pd.read_csv(fns.loc[('ERA5-025','VCV','Agespec_interaction')].path)
                        df.columns = col_idxs.set_names(['pred_1','covar_1','ageset_1'])
                        df.index = col_idxs.set_names(['pred_2','covar_2','ageset_2'])
    
                        ds_tmp = (df.melt(ignore_index=False,value_name='vcv').reset_index().
                                        set_index(['pred_1','covar_1','ageset_1','pred_2','covar_2','ageset_2']).to_xarray())
                        
                        ds_tmp = ds_tmp.expand_dims({'proj_base':[proj_base],
                                                                 'spectype':[regtyp]}) 
    
                    dss.append(ds_tmp)
        elif obstreat == 'joint':
            if regtyp == 'pooled':
                # Get order it's encoded in by STATA
                modlist = np.sort([k for k in fns.index.unique(level='proj_base') if k != 'joint'])
                
                col_idxs = pd.MultiIndex.from_product([[f'tassumpoly_{k}' for k in np.arange(1,5)],modlist],
                                                       names = ['pred','proj_base'])
                col_idxs = col_idxs.append(pd.MultiIndex.from_tuples([('intersect','intersect')],names=['pred','proj_base']))
                
                if output == 'coeffs':
                    df = pd.read_csv(fns.loc[('joint','Coefficients',regtyp)].path)
                    df.columns = col_idxs
                    
                    ds_tmp = df.melt(value_name='coeffs').set_index(['pred','proj_base']).to_xarray()
                    ds_tmp = ds_tmp.expand_dims({'covar':['1'],'ageset':['total']})
                    
    
                elif output == 'vcv':
                    df = pd.read_csv(fns.loc[('joint','VCV',regtyp)].path)
                    
                    df.columns = col_idxs.set_names(['pred_1','proj_base_1'])
                    df.index = col_idxs.set_names(['pred_2','proj_base_2'])
            
                    ds_tmp = (df.melt(ignore_index=False,value_name='vcv').reset_index().
                                  set_index(['pred_1','proj_base_1','pred_2','proj_base_2']).to_xarray())
                ds_tmp = ds_tmp.expand_dims({'spectype':[regtyp]}) 
                dss.append(ds_tmp)
            
            elif regtyp == 'Agespec_interaction':
                raise KeyError('joint-reanalysis fully interacted model not run yet')
                        
        dss = xr.merge(dss,join='outer',compat='no_conflicts')
    
        if drop_intersect:
            if output == 'coeffs':
                for dim in ['pred','ageset','proj_base']:
                    dss = dss.sel({dim:[d for d in dss[dim].values if d != 'intersect']})
            elif output == 'vcv':
                dss = dss.sel({f'pred_{k}':[p for p in dss[f'pred_{k}'].values if p != 'intersect'] for k in [1,2]})
                if 'proj_base_1' in dss.sizes:
                    dss = dss.sel({f'proj_base_{k}':[p for p in dss[f'proj_base_{k}'].values if p != 'intersect'] for k in [1,2]})
    
        # Make singleton dimensions empty coordinates instead
        for dim in [d for d in dss.sizes if d not in ['proj_base','pred']]:
            if dss.sizes[dim] == 1:
                dss = dss.isel({dim:0})
        
        return dss

    coeffs,vcv = get_coeffs_mortality('coeffs',obstreat,regtyp),get_coeffs_mortality('vcv',obstreat,regtyp)

    # Turn ERA5-025 into ERA5
    coeffs['proj_base'] = [re.sub(r'^ERA5-025$','ERA5',mod) for mod in coeffs.proj_base.values]
    if 'proj_base_1' in vcv.sizes:
        for k in [1,2]:
            vcv[f'proj_base_{str(k)}'] = [re.sub(r'^ERA5-025$','ERA5',mod) for mod in vcv[f'proj_base_{str(k)}'].values]
    else:
        vcv['proj_base'] = [re.sub(r'^ERA5-025$','ERA5',mod) for mod in vcv.proj_base.values]
    

    return coeffs,vcv


def coeffs_payroll(obstreat='sep',popweight=True):
    fns = pd.DataFrame({'path':np.sort(glob.glob(dir_list['regs']+r'payroll/regresults*tasmax*.csv'))})
    fns = fns.iloc[[not re.search(r'\_os\_',fn) for fn in fns['path'].values]]
    fns['proj_base'] = [re.split(r'\_',re.split(r'\/',fn)[-1])[3] for fn in fns['path'].values]
    fns['wghttyp'] = [re.split(r'\.csv',re.split(r'\_',re.split(r'\/',fn)[-1])[4])[0] for fn in fns['path'].values]
    fns = fns.set_index(['proj_base','wghttyp'])

    def get_coeffs_payroll(output,obstreat,popweight=popweight,fns=fns):
        '''

        '''

        if output not in ['coeffs','vcv']:
            raise KeyError(f'`output` must be "coeffs" or "vcv", but is "{str(output)}"')
        
        if obstreat not in ['sep','joint']:
            raise KeyError(f'`obstreat` must be "sep" or "joint", but is "{str(obstreat)}"')

        if popweight:
            wghttyp = 'pop'
        else:
            wghttyp = 'raw'
        dss = []
        
        if obstreat == 'sep':
            for proj_base in [mod for mod in fns.index.unique(level='proj_base') if mod != 'joint']:
        
                df = pd.read_csv(fns.loc[(proj_base,wghttyp),'path'],sep='\t',skiprows=1)
                df.columns = ['param','vals']
                df['varname'] = ['b','se','p']*int(len(df)/3)
                df['param'] = np.repeat(list(df.loc[np.arange(0,len(df),3),'param'].values),3)
                df = df.pivot(index='param',values = 'vals',columns ='varname')
        
                df = df.loc[['days90plus']]
        
                # Put in format that can be propagated easily, 
                # across all temperature bins above 90 degrees
                bin_names = [f'tasmaxbinF_{str(bin)}' for bin in np.arange(92.5,127.6,5)]
                df = pd.concat([df.rename({'days90plus':bn}) for bn in bin_names])
        
                # To xarray
                ds_tmp = df[['b','se']].to_xarray()
                ds_tmp = ds_tmp.rename({'b':'coeffs','se':'vcv','param':'pred'})
        
                if output == 'coeffs':
                    ds_tmp = ds_tmp[['coeffs']]
                elif output == 'vcv':
                    # Since just a single predictor (days > 90) that we're 
                    # awkwardly expanding out across bins to make the syntax
                    # clearer in `3_project_impacts`, make the VCV a 
                    # corresponding diagonal matrix
                    ds_tmp = xr.DataArray(np.diag(ds_tmp['vcv'].values),
                     dims = ['pred_1','pred_2'],
                     coords = {'pred_1':(('pred_1'),ds_tmp.pred.values),
                               'pred_2':(('pred_2'),ds_tmp.pred.values)}).to_dataset(name='vcv')
        
                 # Add proj_base dimension
                ds_tmp = ds_tmp.expand_dims({'proj_base':[proj_base]})
        
                # Add to full dataset
                dss.append(ds_tmp)
        
            # Make into single ds
            dss = xr.merge(dss,join='outer',compat='no_conflicts')
        
            return dss

    coeffs,vcv = get_coeffs_payroll('coeffs',obstreat,popweight),get_coeffs_payroll('vcv',obstreat,popweight)

    # Turn ERA5-025 into ERA5
    coeffs['proj_base'] = [re.sub(r'^ERA5-025$','ERA5',mod) for mod in coeffs.proj_base.values]
    if 'proj_base_1' in vcv.sizes:
        for k in [1,2]:
            vcv[f'proj_base_{str(k)}'] = [re.sub(r'^ERA5-025$','ERA5',mod) for mod in vcv[f'proj_base_{str(k)}'].values]
    else:
        vcv['proj_base'] = [re.sub(r'^ERA5-025$','ERA5',mod) for mod in vcv.proj_base.values]

    return coeffs,vcv

rng = np.random.default_rng()
def load_dmgfs(proc,load_kwargs = {},
              rng = rng,
              nsamples = 1000):
    ''' Load MC'ed damage function, sampled at temperatures Tref

    Parameters
    ---------------
    proc : str
        Name of impact variable

    load_kwargs : dict
        Piped into `coeffs_*()`

    rng : numpy random generator

    nsamples : int
        Number of MC samples of VCV to use

    Returns
    ---------------
    dmgf : xr.DataArray
        `Tref` x `proj_base` x `sidx` datarray giving 
        coefficients at temperatures `Tref`


    '''
    
    # Load regression coefficients
    if proc == 'gdppc':
        coeffs,vcv = coeffs_gdppc(**load_kwargs)
        coeffs = coeffs.sel(timeframe='allyears')
        vcv = vcv.sel(timeframe='allyears')
    elif proc == 'mortality':
        coeffs,vcv = coeffs_mortality(**load_kwargs)
    elif proc == 'payroll':
        coeffs,vcv = coeffs_payroll(**load_kwargs)
    else:
        raise KeyError()

    # Get MC samples of coefficients
    if 'proj_base_1' in vcv:
        # If jointly estimated across parameters
        coeffs_tmp = coeffs.coeffs.stack(predproj = ['pred','proj_base'])
        vcv_tmp = vcv.vcv.stack(predproj_1 = ['pred_1','proj_base_1'],
                                 predproj_2 = ['pred_2','proj_base_2'])
        
        coeffs_mvn = rng.multivariate_normal(coeffs_tmp,
                                             vcv_tmp,
                                             size=nsamples)
        
        coeffs_mvn = xr.DataArray(coeffs_mvn,
                                  dims = ['sidx','predproj'],
                                  coords = {'predproj':coeffs_tmp.predproj})
        coeffs_mvn = coeffs_mvn.unstack('predproj')
    else:
        # If separate fit for each proj_base
        coeffs_mvn = xr.apply_ufunc(rng.multivariate_normal,
                                   coeffs.coeffs,
                                   vcv.vcv,
                                   input_core_dims = [['pred'],['pred_1','pred_2']],
                                   output_core_dims = [['sidx','pred']],
                                   kwargs = {'size':nsamples},
                                           dask = 'parallelized',
                                           dask_gufunc_kwargs = {'output_sizes':{'sidx':nsamples}},
                                   vectorize=True)


    if 'sumpoly' in coeffs_mvn.pred.values[0]:
        ts = np.arange(-20,51)
        var_stem = 'tassumpoly_'
        poly_degs = [int(re.split(var_stem,v)[-1]) for v in coeffs_mvn.pred.values if var_stem in v]
        t_vals = xr.DataArray(np.vstack([ts**k for k in poly_degs]),
                     dims = ['pred','bin_c'],
                     coords = {'pred':[v for v in coeffs_mvn.pred.values if var_stem in v],
                               'bin_c':ts})
    
        # Get damage function
        dmgf = (coeffs_mvn*t_vals).sum('pred')
        
        if proc == 'mortality':
            # Normalize over lowest value
            dmgf = dmgf - dmgf.sel(bin_c = slice(10,30)).min('bin_c')
        elif proc == 'gdppc':
            # Normalize over highest value
            dmgf = dmgf - dmgf.sel(bin_c = slice(10,30)).max('bin_c')
            #dmgf = dmgf - dmgf.sel(bin_c = 13)
        
        dmgf = dmgf.rename({'bin_c':'Tref'})
    elif 'bin' in coeffs_mvn.pred.values[0]:
        # Constructing the single-bin version
        dmgf = xr.apply_ufunc(rng.normal,
                   coeffs.isel(pred=0),
                   vcv.vcv.isel(pred_1 = 0,pred_2 = 0),
                   output_core_dims = [['sidx']],
                   vectorize=True,
                   kwargs = {'size':nsamples})
        
        dmgf = dmgf.expand_dims({'bin_f':np.arange(float(re.split(r'tasmaxbinF\_',str(coeffs.pred.values[0]))[-1]),
                                             125,5)})
        
        dmgf = xr.concat([xr.zeros_like(dmgf.isel(bin_f=0,drop=True)).expand_dims({'bin_f':np.arange(float(re.split(r'tasmaxbinF\_',str(coeffs.pred.values[0]))[-1])-20,
                                                                                   float(re.split(r'tasmaxbinF\_',str(coeffs.pred.values[0]))[-1]),
                                                                                   5)}),
                          dmgf],
                         dim='bin_f')
    
        dmgf = dmgf.rename({'bin_f':'Tref'}).coeffs

    return dmgf
