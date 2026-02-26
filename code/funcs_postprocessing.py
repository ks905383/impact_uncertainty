# Functions to put regression results into the right form


import xarray as xr
import numpy as np
import pandas as pd
import re

def read_csvv_new(path,store_meta = False,
              remove_constant = True):
    ''' Turns CSVV files of regression results into xarray datasets

    Parameters
    --------------
    path : :py:meth:`str`
        Path to csvv file.

    store_meta : :py:meth:`bool`, by default False
        Whether to store parameter metadata as variables (off by default 
        because since the variables are a dimension ('pred') in the final 
        dataset, this involves in the creation of extra variables just 
        for the metadata. A little clunky, basically).

    remove_constant : :py:meth:`bool`, by default True
        If True, removes the constant term from the gamma and vcv matrix. 
        Useful if you're not using it for anything, since it's a bit of 
        a clunky add-on in the xarray context.

    Returns
    ---------------
    ds : :py:meth:`xr.Dataset`

    ds_vcv : :py:meth:`xr.Dataset`
    '''

    with open(path, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('------')]
    
    attrs = {}
    variables = {}
    section = None
    
    #-------------------------------------------
    #--------- Scrape data  --------------------
    #-------------------------------------------
    
    # Expected section titles
    section_list = ['observations','prednames','covarnames','gamma',
                    'gammavcv','residvcv','FE','coefnames']
    
    # Get structure for section contents
    content = {v:[] for v in section_list if v != 'observations'}
    
    for line in lines:
        if (line == '---') or (line == '...'):
            continue
        # Header metadata
        if ':' in line and section is None:
            key, val = line.split(':', 1)
            key = key.strip().lower()
            val = val.strip()
            if key == "variables":
                section = "variables"
            else:
                attrs[key] = val
            continue
    
        # If the line is a section title, then set the section
        # and move on to the next line
        if line.lower() in section_list:
            section = line.lower()
            continue
        
        # Section parsing
        if section == "variables":
            if ':' in line:
                var, rest = line.split(':', 1)
                match = re.match(r"(.*)\s+\[(.*)\]$", rest.strip())
                longname = rest.strip()
                units = ''
                if match:
                    longname, units = match.groups()
                variables[var.strip()] = {'longname': longname.strip(), 'units': units.strip()}
        elif section == 'observations':
            attrs['nobs'] = int(line)
        elif section in ['prednames','covarnames','coefnames']:
            content[section] += [k.strip(',') for k in line.split()]
        elif section in ['gamma','gammavcv','residvcv']:
            vals = []
            for val in line.strip(', ').split(','):
                try: 
                    vals.append(float(val))
                except:
                    vals.append(np.nan)
            content[section] += vals
    
    # Determine shapes
    var_lengths = {v:len(content[v]) for v in ['prednames','coefnames','gamma']}
    
    #-------------------------------------------
    #--------- Manage indices  -----------------
    #-------------------------------------------
    # Turn "coefnames" (i.e., identifiers of gammas) into dataframe
    df_gamma = pd.DataFrame([re.split(r'\#',cn) for cn in content['coefnames']])
    df_gamma.columns = ['idx'+str(n) for n in df_gamma.columns]
    # Add gammas
    df_gamma = df_gamma.join(pd.DataFrame(content['gamma'],columns=['gamma']))
    
    # Group indices by type (to allow simplifying names, removing 
    # STATA-specific indicators, etc.)
    idx_defs = {}
    idx_names = [c for c in df_gamma.columns if 'idx' in c]
    idx_renames = {c:c for c in df_gamma.columns}
    for idxn in idx_names:
        # Figure out if a given column is categorical variables, 
        # by flagging all where the index is [##](b).[name] (b is base
        # in STATA's estout parlance)
        if np.all([bool(re.search(r'(^[0-9].*[b]{0,1}\.)|(\_cons)',cidx)) 
                   for cidx in df_gamma.loc[:,idxn].values if cidx is not None]):
            idx_defs[idxn] = 'factor'
    
        # Figure out if a given column is continuous (and thus the "c." 
        # can be removed) 
        elif np.all([bool(re.search(r'^c\.',cidx)) 
                             for cidx in df_gamma.loc[:,idxn].values if cidx is not None]):
            idx_defs[idxn] = 'cont'
        else:
            idx_defs[idxn] = 'other'
    
    # Rename factor index 
    for idxn in idx_names:
        if idx_defs[idxn] == 'factor':
            factorvar_name = [re.split(r'\.',k)[1] for k in df_gamma.loc[:,idxn].values if k != '_cons']
            factorvar_name = np.unique(factorvar_name)
            # If truly unique, single (allowing for _cons as an option) factor
            # name identifier, use that to simplify column
            if len(factorvar_name) == 1:
                factorvar_name = factorvar_name[0]
    
                # Remove common variable name from individual index values
                df_gamma.loc[:,idxn] = [re.sub(r'\.'+factorvar_name,'',k) for k in df_gamma.loc[:,idxn].values]
    
                # Replace column name
                idx_renames[idxn] = factorvar_name
            else:
                continue
        elif idx_defs[idxn] == 'cont':
            # Remove the "c." at the start, for convenience
            df_gamma.loc[:,idxn] = [re.sub(r'^c\.','',k) for k in df_gamma.loc[:,idxn].astype(str).values]
    
        # If all the predictors are in this column, this is the pred var column
        if np.all([pn in df_gamma.reset_index().loc[:,idxn].values 
                   for pn in content['prednames']]):
            # So, replace column name
            idx_renames[idxn] = 'pred'
        
            # Make 'None' and None into '_cons'
            df_gamma.loc[:,idxn] = [re.sub(r'^None$','_cons',k) for k in df_gamma.loc[:,idxn].astype(str).values]
    
        # If all the covars are in this column, this is the covvar column
        if np.all([re.sub('^1$','None',pn) in df_gamma.reset_index().loc[:,idxn].values 
                   for pn in content['covarnames']]):
            # So, replace column name
            idx_renames[idxn] = 'covar'
    
            # Make 'None' and None into '1'
            df_gamma.loc[:,idxn] = [re.sub(r'^None$','1',k) for k in df_gamma.loc[:,idxn].astype(str).values]
    
    # Rename columns, if necessary
    df_gamma.columns = [idx_renames[c] for c in df_gamma.columns]
    
    # Set updated index columns as index
    df_gamma = df_gamma.set_index([idx_renames[idxn] for idxn in idx_names])
    
    # Change idx_names to the renamed ones, to use later if needed
    idx_names = [idx_renames[k] for k in idx_names]
    
    # Create dataset for gammas
    ds_gamma = df_gamma.to_xarray()
    
    # Use gamma index columns to label vcv matrix columns/rows
    df_gammavcv = pd.DataFrame(np.reshape(content['gammavcv'],(var_lengths['gamma'],var_lengths['gamma'])),
                                index = df_gamma.index,
                                 columns = df_gamma.index)
    
    # Create dataset for vcv matrix
    ds_gammavcv = xr.DataArray(
        df_gammavcv.values,
        dims=['gamma_flat_i', 'gamma_flat_j'],
        coords={
            'gamma_flat_i': pd.MultiIndex.from_tuples(df_gammavcv.index, names=[k+'_i' for k in df_gammavcv.index.names]),
            'gamma_flat_j': pd.MultiIndex.from_tuples(df_gammavcv.columns, names=[k+'_j' for k in df_gammavcv.columns.names]),
        },
        name='gammavcv')
    
    #-------------------------------------------
    #--------- Output  -------------------------
    #-------------------------------------------
    # Get unified ds
    ds_out = xr.merge([ds_gamma,ds_gammavcv])
    
    # Housekeeping
    ds_out.attrs = attrs
    
    # Optionally store variable metadata in attributes
    if store_meta:
        for varname, meta in variables.items():
            ds_out[f'varmeta_{varname}'] = xr.DataArray([], attrs=meta)
            
    # Remove constant, if desired 
    if remove_constant:
        # remove constant from gammas
        ds_out = ds_out.sel({k:[kv for kv in ds_out[k].values if kv != '_cons']
                             for k in idx_names})
    
        # remove constant from vcv
        ds_out = ds_out.sel({k:ds_out[k].values[['_cons' not in kv for kv in ds_out[k].values]]
                             for k in ds_out.dims if 'gamma_flat' in k})
    
    # Now, split out the vcv 
    ds_out_vcv = ds_out[['gammavcv']].copy()
    ds_out = ds_out.drop_vars(['gammavcv','gamma_flat_j','gamma_flat_i',
                               *[k+'_i' for k in idx_names],
                               *[k+'_j' for k in idx_names]])
    
    return ds_out,ds_out_vcv

def read_csvv(path,store_meta = False,
              dim_order = ['pred','covar','ageset'],
              covarnames = None,
              remove_constant = True):
    ''' Turns CSVV files of regression results into xarray datasets

    Parameters
    --------------
    path : :py:meth:`str`
        Path to csvv file.

    store_meta : :py:meth:`bool`, by default False
        Whether to store parameter metadata as variables (off by default 
        because since the variables are a dimension ('pred') in the final 
        dataset, this involves in the creation of extra variables just 
        for the metadata. A little clunky, basically).

    dim_order = ['pred','covar','ageset']
        Manually set the stacking order for the betas, vcv matrix. Needed
        because the csvv does not separate out 'ageset' as a dimension, 
        and thus you can't get it from metadata alone

    covarnames : by default None
        Temporary hack while the csvv generation code doesn't print the 
        covarnames. Should be ['1','climtas','loggdppc'], I believe. 

    Returns
    ---------------
    ds : :py:meth:`xr.Dataset`

    ds_vcv : :py:meth:`xr.Dataset`

    '''
    
    with open(path, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('------')]
    
    attrs = {}
    variables = {}
    section = None

    # Expected section titles
    section_list = ['observations','prednames','covarnames','gamma',
                    'gammavcv','residvcv','FE','coefnames']

    # Get structure for section contents
    content = {v:[] for v in section_list if v != 'observations'}
    
    for line in lines:
        if (line == '---') or (line == '...'):
            continue
        # Header metadata
        if ':' in line and section is None:
            key, val = line.split(':', 1)
            key = key.strip().lower()
            val = val.strip()
            if key == "variables":
                section = "variables"
            else:
                attrs[key] = val
            continue
    
        # If the line is a section title, then set the section
        # and move on to the next line
        if line.lower() in section_list:
            section = line.lower()
            continue
        
        # Section parsing
        if section == "variables":
            if ':' in line:
                var, rest = line.split(':', 1)
                match = re.match(r"(.*)\s+\[(.*)\]$", rest.strip())
                longname = rest.strip()
                units = ''
                if match:
                    longname, units = match.groups()
                variables[var.strip()] = {'longname': longname.strip(), 'units': units.strip()}
        elif section == 'observations':
            attrs['nobs'] = int(line)
        elif section in ['prednames','covarnames','coefnames']:
            content[section] += [k.strip(',') for k in line.split()]
        elif section in ['gamma','gammavcv','residvcv']:
            vals = []
            for val in line.strip(', ').split(','):
                try: 
                    vals.append(float(val))
                except:
                    vals.append(np.nan)
            content[section] += vals

    # Add covarnames in hack, if not readable from csvv
    if covarnames is not None:
        content['covarnames'] = covarnames
    else:
        # If the covarnames slot is empty, fill it with '1' (i.e., level, no interaction)
        if len(content['covarnames']) == 0:
            content['covarnames'] = '1'
    
    # Determine shapes
    var_lengths = {v:len(content[v]) for v in ['prednames','covarnames','gamma']}
    if 'Age‑combined' not in attrs['oneline']:
        extra_dim = {'ageset':['young','older','oldest']}
    else:
        extra_dim = {'ageset':['combined']}
    extra_dim_name= [k for k in extra_dim][0]
    var_lengths['extra'] = len(extra_dim[extra_dim_name])

    # If the gammas are exactly one longer than the number of named
    # variables, then there's a constant in here 
    if var_lengths['gamma'] % var_lengths['prednames'] == 1:
        content['prednames'] += ['constant']
    
        # Pad the 'constant' slot with extra nans, if there is
        # more than one covariate
        content['gamma'] += [np.nan]*(var_lengths['covarnames'] - 1)
    
        # Pad the gammavcv slot with extra nans, if there is 
        # more than one covariate. This is a bit more complicated,
        # since it means adding nans every len(preds)+2:len(preds)+len(covars)
        gammavcv_withnans = []
        for slot in range(0,len(content['gammavcv']),var_lengths['prednames']*var_lengths['covarnames']+1):
            # First, copy the slot of gammas into the new list
            gammavcv_withnans.extend(content['gammavcv'][slot:(slot+var_lengths['prednames']*var_lengths['covarnames']+1)])
            # Now, add padded nans
            gammavcv_withnans.extend([np.nan]*(var_lengths['covarnames'] - 1))
        # Add remaining rows of nans at the end
        gammavcv_withnans.extend([np.nan]*(((var_lengths['prednames']+1)*var_lengths['covarnames'])*(var_lengths['covarnames'] - 1)))
        
        # Now, replace the original vcv with the new padded one
        content['gammavcv'] = gammavcv_withnans
    
        # Fix varlengths
        for v in ['prednames','gamma']:
            var_lengths[v] = len(content[v])

    assert var_lengths['prednames'] * np.max([1,var_lengths['covarnames']]) * np.max([1,var_lengths['extra']]) == var_lengths['gamma']

    # Create xarray.Dataset
    ds = xr.Dataset(
        data_vars={
            'gamma': ([extra_dim_name,'pred', 'covar'], 
                      np.array(content['gamma']).reshape(*[var_lengths[v]                                   
                                                           for v in ['extra','prednames','covarnames']])),
            'gammavcv': (['gamma_flat_i', 'gamma_flat_j'], 
                         np.array(content['gammavcv']).reshape(var_lengths['gamma'],var_lengths['gamma'])),
            'residvcv': ([], np.array(content['residvcv'])[0] 
                         if len(content['residvcv']) == 1 
                         else np.array(content['residvcv']))
        },
        coords={
            'pred': content['prednames'],
            'covar': np.atleast_1d(content['covarnames']),
            extra_dim_name : extra_dim[extra_dim_name],
            'gamma_flat_i': np.arange(var_lengths['gamma']),
            'gamma_flat_j': np.arange(var_lengths['gamma']),
        },
        attrs=attrs
    )

    # Optionally store variable metadata in attributes
    if store_meta:
        for varname, meta in variables.items():
            ds[f'varmeta_{varname}'] = xr.DataArray([], attrs=meta)

    # Now, split out the vcv 
    ds_vcv = ds[['gammavcv']].copy()
    ds = ds.drop_vars(['gammavcv','gamma_flat_j','gamma_flat_i'])

    # Make multiindex identifying elements of the VCV matrix
    for k in ['i','j']:
        ds_vcv['gamma_flat_'+k] = (ds[dim_order].stack({'gamma_flat_'+k :dim_order})['gamma_flat_'+k].
                                   rename({d:d+'_'+k for d in dim_order}))

    return ds,ds_vcv


def wrapper_mvn(mean,cov,ndraws,rng=np.random.default_rng()):
    ''' Draw coefficients from a multivariate normal distribution

    (an apply_ufunc-compatible wrapper for :py:meth:`rng.multivariate_normal`)

    Sample usage: 

        xr.apply_ufunc(wrapper_mvn,
                       ds.gamma,
                       ds_vcv.gammavcv,
                       ndraws,
                       input_core_dims = [['pred','covar','ageset'],['gamma_flat_i','gamma_flat_j'],[]],
                       output_core_dims = [['ndraws','pred','covar','ageset']],
                       vectorize=True)
        
    Parameters
    --------------
    mean : coefficient estimates (N x 1) (or split up into other dimensions, 
           that the final coefficients are reshaped to)

    cov : variance-covariance matrix (N x N)

    ndraws : how many draws to get

    rng : random seed

    Returns
    --------------
    

    '''
    
    # Get dims of coefficient matrix
    shapes = mean.shape

    # Remove nans, assuming that the nans are padded at the end
    mean = mean.flatten()
    if np.any(np.isnan(mean)):
        nonan_len = np.argmax(np.isnan(mean))
        nnans = len(mean) - nonan_len
        mean = mean[0:nonan_len]
        cov = cov[0:nonan_len,0:nonan_len]
    else:
        nnans = 0

    # Draw betas
    out = rng.multivariate_normal(mean.flatten(),cov,ndraws)

    # Pad back nans if needed
    out = np.pad(out,((0,0),(0,nnans)),constant_values=np.nan)
    
    # Reshape back into original dimensions
    out = out.reshape(ndraws,*shapes)
    
    return out