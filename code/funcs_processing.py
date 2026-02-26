import xarray as xr
import numpy as np
import pandas as pd


def calc_uncerts(da,dim_order = ['sidx','proj_base','exprun','model'],
                 names = {'sidx':'uncreg','proj_base':'uncobs','exprun':'uncint','model':'uncmod'}):
    ''' Calculate variances

    Parameters
    -------------
    da : :py:class:`xr.DataArray`

    dim_order : :py:class:`list`

    names : :py:class:`dict`

    '''
    # Figure out which dimensions are in the dataset
    dim_order = [dim for dim in dim_order if dim in da.sizes]

    # Get functions to apply
    funcs = [np.roll([np.var,*([np.mean]*(len(dim_order)-1))],shift=shift) for shift in range(len(dim_order))]

    # Get names of output variables
    varnames = [names[dim] for dim in dim_order]

    # Apply 
    def apply_functions_sequentially(da, funcs, dims):
        """
        Apply a sequence of functions to specified dimensions in order.
    
        Parameters:
        - da (xarray.DataArray or xarray.Dataset): Input data.
        - func_dim_mapping (dict): Ordered mapping {function: dimension} specifying 
          the operations to apply in sequence.
    
        Returns:
        - xarray.DataArray or xarray.Dataset: Transformed data after sequential application.
        """
        for func, dim in zip(funcs,dims):
            da = func(da, axis=da.get_axis_num(dim))
        return da
    
    uncerts = xr.merge([apply_functions_sequentially(da,funclist,dim_order).to_dataset(name=varname)
                        for funclist,varname in zip(funcs,varnames)])

    # Get each as a fraction of their sum (`.to_array()` turns the dataset into a dataarray
    # with dimension 'variable')
    uncerts_norm = uncerts / uncerts.to_array().sum('variable')


    uncerts = xr.concat([uncerts,uncerts_norm],
          dim = pd.Index(['raw','norm'],name='typ'))
    

    return uncerts
