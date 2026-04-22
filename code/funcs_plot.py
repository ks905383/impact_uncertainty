import xarray as xr
import pandas as pd
import geopandas as gpd
import numpy as np
from matplotlib import patches as mpatches
from matplotlib import pyplot as plt
import matplotlib as mpl
import cmocean
from cartopy import crs as ccrs

colors = {'ERA5':'tab:red','JRA-3Q':'tab:orange','MERRA2':'tab:green','GMFD':'tab:blue','PRISM':'yellow','UDel':'yellow'}

proc_names = {'mortality':'mortality rate [/100,000]',
              'gdppc':'log GDP per capita growth',
              'payroll':'payroll'}
proc_titles = {'mortality':'annual mortality rate',
               'gdppc':'GDP per capita growth',
               'payroll':'annual payroll'}
proc_xlims = {'mortality':[-15,45],
              'gdppc':[-2.5,32.5],
              'payroll':[70,115]}
proc_ylims = {'mortality':[-1,8],
              'gdppc':[-0.375,0.025],
              'payroll':[-0.00075,0.0001]}
def plot_dmgfs(dmgf,proc,ax=None,
               colors = colors,
               add_grid = True,
               proc_names = proc_names,
               proc_titles = proc_titles,
               proc_xlims = proc_xlims,
               proc_ylims = proc_ylims,
               simple = False):
    '''
    Damage function plotting code 

    '''

    if ax is None:
        ax = plt.subplot()
    
    for proj_base in dmgf.proj_base.values:
        plot_data = dmgf.sel(proj_base = proj_base)
        ax.fill_between(plot_data.Tref,
                        plot_data.quantile(0.025,dim='sidx'),
                        plot_data.quantile(0.975,dim='sidx'),
                        color=colors[proj_base],alpha=0.1)
        plot_data.mean('sidx').plot(color=colors[proj_base],linewidth=1,label=proj_base,ax=ax)
    
    
    ax.set_xlim(proc_xlims[proc])
    ax.set_ylim(proc_ylims[proc])

    if not simple:
        if add_grid:
            plt.grid()

    if simple:
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title('')
    else:    
        ax.set_xlabel(r'Temperature')# [$^\circ$C]')
        ax.set_ylabel('Change in '+proc_names[proc])
        ax.set_title('Impact of day at T on\n'+proc_titles[proc])

    ax.axhline(0,color='k',linewidth=1)

    return ax 

def add_colorbar(cbar_params, 
                 orientation='vertical', label=None, 
                 axpos=None,
                 fig=None,
                 **kwargs):
    """
    Create a colorbar on its own axis using parameters in cbar_params.

    Parameters
    ----------
    cbar_params : dict
        Dict containing at least {'vmin', 'vmax', 'cmap'}.
        Optionally may include 'levels' (either explicit or as nlevels) for 
        discrete shading.
    orientation : str
        'vertical' (default) or 'horizontal'.
    label : str, optional
        Label for the colorbar.
    axpos : 
    **kwargs :
        Additional kwargs passed to plt.colorbar().

    Returns
    -------
    cbar : matplotlib.colorbar.Colorbar
    """
    if axpos is None:
        if orientation == 'horizontal':
            axpos = [0.15,0.065,0.7,0.025]
        elif orientation == 'vertical':
            axpos = [0.875, 0.15, 0.025, 0.7]
    cax = fig.add_axes(axpos)

    if fig is None:
        fig = plt.gcf()

    
    vmin = cbar_params['vmin']
    vmax = cbar_params['vmax']
    cmap = cbar_params['cmap']

    if 'levels' in cbar_params:
        # Discrete colormap
        levels = cbar_params['levels']

        if type(levels) == int:
            levels = mpl.ticker.MaxNLocator(nbins=cbar_params['levels']).tick_values(cbar_params['vmin'],cbar_params['vmax'])
            
        norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    else:
        # Continuous colormap
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    # Create a ScalarMappable and colorbar
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(sm, cax=cax, orientation=orientation, **kwargs)

    if label is not None:
        cbar.set_label(label)

    return cbar

# Function to convert integer to Roman values
def printRoman(number):
    # from https://www.geeksforgeeks.org/python-program-to-convert-integer-to-roman/
    num = [1, 4, 5, 9, 10, 40, 50, 90,
        100, 400, 500, 900, 1000]
    sym = ["I", "IV", "V", "IX", "X", "XL",
        "L", "XC", "C", "CD", "D", "CM", "M"]
    i = 12

    output = ''
    while number:
        div = number // num[i]
        number %= num[i]
 
        while div:
            output = output + sym[i]
            div -= 1
        i -= 1

    return output


def plot_map(da,gdf,proj=ccrs.PlateCarree(),fig=None,ax=None,cat_cmap=None,join_by_color = False,**plot_kwargs,
             ):
    '''

    Parameters
    -------------

    cat_cmap : dict 
        
    '''
    if fig is None:
        fig = plt.figure()
    
    gdf_tmp = gdf.copy()

    gdf_tmp['plotvar'] = da.values.flatten()

    if ax is None:
        ax = plt.subplot(projection=proj)

    if cat_cmap is None:    
        ax = gdf_tmp.plot(column = 'plotvar',ax=ax,transform=ccrs.PlateCarree(),**plot_kwargs)
    else:
        for label, group in gdf_tmp.groupby('plotvar'):
            ax = group.plot(ax=ax,transform=ccrs.PlateCarree(),
                              color = cat_cmap[label],
                              **plot_kwargs)
        
    return fig,ax

def create_segcmap(plot_data,cbar_params,return_all = False,centered = False):
    '''

    Parameters
    ------------------
    plot_data : 1-D xr.DataArray

    cbar_params : dict
        Containing:
            'vmin'
            'vmax'
            'levels'
            'cmap'

    return_all : bool
        If True, returns bin bounds and discrete 
        colors used

    centered : bool
        If True, then `vmin` and `vmax` are 
        assumed to be the _center_ of the first
        and last bins

    Returns
    ------------------
    colors : xr.DataArray
        With the same dims, coords as plot_data

    bounds : np.Array
        If `return_bounds` only

    colors_discrete : list
        If `return_bounds` only
    '''
    if centered:
        range = cbar_params['vmax'] - cbar_params['vmin']
        binwidth = (range + 1)/cbar_params['levels']
        
        bounds = np.linspace(cbar_params['vmin']-binwidth/2,cbar_params['vmax']+binwidth/2,
                             cbar_params['levels']+1)
    else:
        bounds = np.linspace(cbar_params['vmin'],cbar_params['vmax'],cbar_params['levels']+1)
    
    # Create segmented colormap
    norm = mpl.colors.BoundaryNorm(bounds,cbar_params['cmap'].N)

    # Get color of each value of plot_data
    colors = [mpl.colors.to_hex(cbar_params['cmap'](norm(v))) for v in plot_data.values]
    # Turn into xarray to match up dimensions
    colors = xr.DataArray(colors,
             dims = plot_data.dims,
             coords = plot_data.coords)

    if return_all:
        # Get reference colors used (as the color of the 
        # center of each discrete bin)
        colors_discrete = [cbar_params['cmap'](norm(v)) for v in bounds[0:-1]+np.diff(bounds)[0]/2]

        return colors,bounds,colors_discrete
    else:
        return colors
    


def plot_bauble_histogram(data,x,hue=None,
                          bins='auto',binwidth=None,
                          palette=None,
                          cmap=None,
                          xlims = None,
                          adjust_xlims = True,goal_aspect = 1.2,
                          yticks=True,
                          legend=True,
                          legend_kwargs = {},
                          ax=None):
    """ Plot bauble histogram
    
    Plot a histogram that explicitly shows the number of data points in each bin,
    as a 'bauble', stacked on top of the other baubles in each bin. 
    
    NB: The inputs are vaguely modeled after `seaborn`'s inputs, so you'll find
    many of the commands to, e.g., `sns.histplot()` etc. copied here as well.
    
    Parameters
    -------------------
    data : :class:`pandas.DataFrame`
        Input data structure
        
    x : key in ``data``
        Variable that specifies positions on the x axis
        
    hue : None or key in ``data`
        Which column to use to assign colors 
        
    palette : None or list of strs, floats, :class:`matplotlib.colors.Colormap`
        Method for choosing colors when using `hue`. If None,
        then the palette is automatically created using `cmap`. 
        
    cmap : None or :class:`matplotlib.colors.Colormap`, by default cmocean's "phase"
        Which colormap to use to assign colors if no explicit palette
        is applied. If None, then `cmocean.cm.phase` is used as the colormap; 

    bins : str, number, vector, or a pair of such values; by default 'auto'
        Passed to :func:`numpy.histogram_bin_edges`, can be the name of a 
        reference rule, the number of bins, or the breaks of the bins.
        
    binwidth : number; by default `None`
        Alternate to `bins` input; used to explicitly set binwidths
        
    xlims : None or [float,float]
        Used to set custom xlimits if desired. These are used to set 
        binwidths as well, if binwidth is not None. 
        
    adjust_xlims : bool, by default True
        If True, then the aspect ratio of the plot is set to `goal_aspect`
        
    goal_aspect : float, by default 1.2
        If `adjust_xlims==True`, then the width of the plot is set to 
        `goal_aspect*[(max obs in a given bins) + 2]`. 
        
    yticks : bool, by default True
        If False, then yticks (showing the number of obs in a given 
        bin) are removed. 
        
    legend : bool, by default True
        If True, then a legend is added
        
    legend_kwargs : dict, by default {}
        Piped into `ax.legend()`
        
    
    
    Returns
    -------------------
    ax : the map axis
    
    
    Author: Kevin Schwarzwald
    Last update: 07/11/2023
    License: CC BY-SA
    """
    #----------- Input checks -----------
    if type(data) != pd.DataFrame:
        raise NotImplementedError('`data` must be a `pd.DataFrame`, alternate methods not yet supported.')
        
    if type(x) != str:
        raise NotImplementedError('`x` must be a key to `data`, alternate methods not yet supported.')
        
    if hue is not None:
        if type(hue) != str:
            raise NotImplementedError('`hue` must be a key to `data`, alternate methods not yet supported.')
        

    #----------- Setup -----------
    # Get rid of index
    data = data.reset_index()

    # Set placeholder column if no hue is set
    if hue is None:
        data['hue'] = 'all'
        hue = 'hue'
        legend = False
    # Find unique values to color by
    hue_levels = np.unique(data[hue])


    # Set palette if None
    if palette is None:
        if cmap is None:
            palette = cmocean.cm.phase(np.arange(0,len(hue_levels))/len(hue_levels))
        else:
            palette = cmap(np.arange(0,len(hue_levels))/len(hue_levels))

    # Get histogram bin edges
    if xlims is None:
        if binwidth is None:
            bin_edges = np.histogram_bin_edges(data[x],bins=bins)
        else:
            bin_edges = np.arange(data[x].min(),data[x].max(),binwidth)
    else:
        bin_edges = np.arange(*xlims,binwidth)

    # Set binwidth, if not specified
    if binwidth is None:
        binwidth = bin_edges[1]-bin_edges[0]

    # Get empty array for bin population
    bincounts = np.zeros((len(hue_levels),len(bin_edges)-1))
    # Get empty array for marker population
    binmarkers = np.zeros((len(hue_levels),len(bin_edges)-1))

    #----------- Populate bins -----------
    for hue_level,level_idx in zip(hue_levels,np.arange(0,len(hue_levels))):
        bincounts[level_idx,:] = np.histogram(data.loc[data[hue]==hue_level,x],bins=bin_edges)[0]

    #----------- Plot -----------
    if ax is None:
        ax = plt.subplot()

    ## Plot bauble histogram
    # Get one handle per hue level, to use as legend handles
    hdls = [None]*len(hue_levels)
    # Iterate through histogram bins
    for bin_idx in np.arange(0,len(bin_edges)-1):
        # Iterate through hue categories
        for level_idx in np.arange(0,len(hue_levels)):
            # Iterate through number of members of hue category in 
            # a given bin 
            for marker_idx in np.arange(0,bincounts[level_idx,bin_idx]):
                # Set offset from previous hue levels in the same
                # histogram bin
                if level_idx > 0:
                    bauble_offset = np.sum(bincounts[0:level_idx,bin_idx])
                else:
                    bauble_offset = 0

                # Add patch for each point to plot
                hdls[level_idx] = ax.add_patch(mpatches.Ellipse((bin_edges[bin_idx]+binwidth/2,
                                                                    bauble_offset+marker_idx+0.5),
                                                width=0.8*binwidth,height=0.8,facecolor=palette[level_idx]))
    ## Annotate
    # Set axis limits
    if xlims is None:
        ax.set_xlim((bin_edges[0]-binwidth,bin_edges[-1]+binwidth))
        xlims = ax.get_xlim()
    else:
        ax.set_xlim(xlims)
    # Set y limit to 2 bauble-widths above the top 
    ax.set_ylim((0,np.max(np.sum(bincounts,0))+2))
    # Set aspect so the run markers look circular
    ax.set_aspect(binwidth)

    # Adjust xlims to a desired aspect ratio for readabiilty,
    # if desired
    if adjust_xlims:
        ax.set_xlim(np.mean(bin_edges)-goal_aspect*binwidth*(np.max(np.sum(bincounts,0))+2)/2,
                    np.mean(bin_edges)+goal_aspect*binwidth*(np.max(np.sum(bincounts,0))+2)/2)

    # Remove y axis ticks
    if not yticks:
        ax.set_yticks([])

    # Vertical line at 0
    if (xlims[0]<0) and (xlims[-1]>0):
        ax.axvline(0,color='k',zorder=0)

    # Grid and ticks on both sides, for ease of reading
    plt.grid(True)
    ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False,labelsize=10)

    # Set xlabel
    ax.set_xlabel(x)

    ## Legend
    # Add legend
    if legend:
        ax.legend(hdls,hue_levels,ncol=4,**legend_kwargs)
        
    #----------- Return -----------
    return ax

def vals_to_color_triangle(xs,dim='source',col_names = {'r':'internal','g':'scenario','b':'model'}):
    ''' Convert bound values along three relative axes into color triangle colors

    '''
    if type(xs)==xr.core.dataarray.DataArray:
        out_colors = np.array([mpl.colors.to_hex([1,1,1])]*xs.values.shape[1])

        # > 2/3 internal: tab:red
        out_colors[xs.sel({dim:col_names['r']})>(2/3)] = mpl.colors.to_hex([214/255,39/255,40/255])
        # > 2/3 scenario: tab:green
        out_colors[xs.sel({dim:col_names['g']})>(2/3)] = mpl.colors.to_hex([40/255,160/255,44/255])
        # > 2/3 model: tab:blue
        out_colors[xs.sel({dim:col_names['b']})>(2/3)] = mpl.colors.to_hex([31/255,119/255,180/255])

        # 1/3-2/3 internal & 1/3-2/3 scenario: tab:orange
        out_colors[(xs.sel({dim:col_names['r']})>(1/3)) & (xs.sel({dim:col_names['r']})<=(2/3)) & 
                   (xs.sel({dim:col_names['g']})>(1/3)) & (xs.sel({dim:col_names['g']})<=(2/3))] = mpl.colors.to_hex([255/255,187/255,120/255])
        # 1/3-2/3 internal & 1/3-2/3 model: tab:purple
        out_colors[(xs.sel({dim:col_names['r']})>(1/3)) & (xs.sel({dim:col_names['r']})<=(2/3)) & 
                   (xs.sel({dim:col_names['b']})>(1/3)) & (xs.sel({dim:col_names['b']})<=(2/3))] = mpl.colors.to_hex([247/255,182/255,210/255])
        # 1/3-2/3 scenario & 1/3-2/3 model: tab:cyan
        out_colors[(xs.sel({dim:col_names['g']})>(1/3)) & (xs.sel({dim:col_names['g']})<=(2/3)) & 
                   (xs.sel({dim:col_names['b']})>(1/3)) & (xs.sel({dim:col_names['b']})<=(2/3))] = mpl.colors.to_hex([158/255,218/255,229/255])

        # Light blue
        out_colors[(xs.sel({dim:col_names['r']})<(1/3)) & (xs.sel({dim:col_names['g']})<(1/3)) 
                   & (xs.sel({dim:col_names['b']})<=(2/3))] = mpl.colors.to_hex([174/255,199/255,232/255])
        # Light red
        out_colors[(xs.sel({dim:col_names['r']})<=(2/3)) & (xs.sel({dim:col_names['g']})<(1/3)) 
                   & (xs.sel({dim:col_names['b']})<(1/3))] = mpl.colors.to_hex([255/255,152/255,150/255])
        # Light green
        out_colors[(xs.sel({dim:col_names['r']})<(1/3)) & (xs.sel({dim:col_names['g']})<=(2/3)) 
                   & (xs.sel({dim:col_names['b']})<(1/3))] = mpl.colors.to_hex([152/255,223/255,138/255])

        #(this overlaps a bit with the above)
        # where all three are within 0.1, grey
        out_colors[(xs.sel({dim:col_names['r']})>3/12) & (xs.sel({dim:col_names['r']})<5/12) & 
                   (xs.sel({dim:col_names['g']})>3/12) & (xs.sel({dim:col_names['g']})<5/12) &
                   (xs.sel({dim:col_names['b']})>3/12) & (xs.sel({dim:col_names['b']})<5/12)] = mpl.colors.to_hex([0.8,0.8,0.8])

        # White for nans
        out_colors[np.isnan(xs.sel({dim:col_names['r']}))] = mpl.colors.to_hex([1,1,1])
    elif type(xs) in [gpd.geodataframe.GeoDataFrame,pd.core.frame.DataFrame]:
        out_colors = np.array(['white']*len(xs))

        # > 2/3 internal: tab:red
        out_colors[xs[col_names['r']]>(2/3)] = mpl.colors.to_hex([214/255,39/255,40/255])
        # > 2/3 scenario: tab:green
        out_colors[xs[col_names['g']]>(2/3)] = mpl.colors.to_hex([40/255,160/255,44/255])
        # > 2/3 model: tab:blue
        out_colors[xs[col_names['b']]>(2/3)] = mpl.colors.to_hex([31/255,119/255,180/255])

        # 1/3-2/3 internal & 1/3-2/3 scenario: tab:orange
        out_colors[(xs[col_names['r']]>(1/3)) & (xs[col_names['r']]<=(2/3)) & 
                   (xs[col_names['g']]>(1/3)) & (xs[col_names['g']]<=(2/3))] = mpl.colors.to_hex([255/255,187/255,120/255])
        # 1/3-2/3 internal & 1/3-2/3 model: tab:purple
        out_colors[(xs[col_names['r']]>(1/3)) & (xs[col_names['r']]<=(2/3)) & 
                   (xs[col_names['b']]>(1/3)) & (xs[col_names['b']]<=(2/3))] = mpl.colors.to_hex([247/255,182/255,210/255])
        # 1/3-2/3 scenario & 1/3-2/3 model: tab:cyan
        out_colors[(xs[col_names['g']]>(1/3)) & (xs[col_names['g']]<=(2/3)) & 
                   (xs[col_names['b']]>(1/3)) & (xs[col_names['b']]<=(2/3))] = mpl.colors.to_hex([158/255,218/255,229/255])

        # Light blue
        out_colors[(xs[col_names['r']]<(1/3)) & (xs[col_names['g']]<(1/3)) 
                   & (xs[col_names['b']]<=(2/3))] = mpl.colors.to_hex([174/255,199/255,232/255])
        # Light red
        out_colors[(xs[col_names['r']]<=(2/3)) & (xs[col_names['g']]<(1/3)) 
                   & (xs[col_names['b']]<(1/3))] = mpl.colors.to_hex([255/255,152/255,150/255])
        # Light green
        out_colors[(xs[col_names['r']]<(1/3)) & (xs[col_names['g']]<=(2/3)) 
                   & (xs[col_names['b']]<(1/3))] = mpl.colors.to_hex([152/255,223/255,138/255])

        #(this overlaps a bit with the above)
        # where all three are within 0.1, white
        out_colors[(xs[col_names['r']]>3/12) & (xs[col_names['r']]<5/12) & 
                   (xs[col_names['g']]>3/12) & (xs[col_names['g']]<5/12) &
                   (xs[col_names['b']]>3/12) & (xs[col_names['b']]<5/12)] = mpl.colors.to_hex([0.8,0.8,0.8])

        # White for nans
        out_colors[np.isnan(xs[col_names['r']])] = mpl.colors.to_hex([1,1,1])
    
    return out_colors