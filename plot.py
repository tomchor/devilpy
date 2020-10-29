def animate_hv(darr, simulation=None, saveas='', xy_dims=['x', 'y'], 
        cmap=None, clabel='', figsize=(7, 7),
        clim=None, logscale=False, grid=True, interpolation=None, fps=15,
        aspect='equal', cbar=True, dpi=120, **kwargs):
    import xarray as xr
    import holoviews as hv
    if saveas!='':
        import matplotlib
        matplotlib.use('Agg')
#    hv.extension('matplotlib')

    renderer=hv.renderer('matplotlib')
    renderer.dpi=dpi
    renderer.fps=fps

    if type(darr)==xr.core.dataarray.DataArray:
        dsarr = xr.Dataset({clabel:darr})
    elif type(darr)==xr.core.dataarray.Dataset:
        dsarr = darr
    else:
        print('submit dataarray or dataset')

    options = dict(fig_inches=figsize, colorbar=cbar,
            cmap=cmap, interpolation=interpolation,
            clims=clim, show_grid=grid, 
            aspect=aspect, **kwargs)
    if logscale==True:
        from matplotlib.colors import LogNorm
        options['norm']=LogNorm

    print('Generating images ...')
    ds = hv.Dataset(dsarr)
    images = ds.to(hv.Image, xy_dims).options('Image', **options) 

    #-----
    # save output
    print('Saving ...')
    if saveas!='':
        if saveas.endswith('html'):
            renderer.save(images, saveas[:-5])
#            with open(saveas, 'w') as fou:
#                fou.write(renderer.static_html(images))
        elif saveas.endswith('mp4'):
            renderer.save(images, saveas[:-4], 'mp4')
        else:
            renderer.save(images, saveas, saveas.split('.')[-1])

    #-----
    return images



