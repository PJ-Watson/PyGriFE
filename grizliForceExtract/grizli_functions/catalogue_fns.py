import numpy as np
import astropy.units as u
import astropy.io.fits as pf
import astropy.wcs as pywcs
from pathlib import Path
import os
import glob
import sep
import grizli
# from grizli import utils, prep, jwst_utils, multifit, fitting
from grizli import prep, utils
from grizli.pipeline import auto_script
# import eazy

def make_SEP_catalog(
    root="",
    sci=None,
    wht=None,
    threshold=1.8,
    get_background=True,
    bkg_only=False,
    bkg_params={'bw': 32, 'bh': 32, 'fw': 3, 'fh': 3},
    verbose=True,
    phot_apertures=None,
    aper_segmask=False,
    rescale_weight=True,
    err_scale=-np.inf,
    use_bkg_err=False,
    column_case=str.upper,
    save_to_fits=True,
    include_wcs_extension=True,
    source_xy=None, 
    compute_auto_quantities=True,
    autoparams=[2.5, 0.35*u.arcsec, 2.4, 3.8],
    flux_radii=[0.2, 0.5, 0.9],
    subpix=0,
    mask_kron=False,
    max_total_corr=2,
    detection_params=None,
    bkg_mask=None,
    pixel_scale=0.06,
    log=False,
    gain=2000., 
    extract_pixstack=int(3e7),
    sub_object_limit=4096,
    exposure_footprints=None,
    seg_image=None,
    in_dir=None,
    out_dir=None,
    **kwargs,
):
    """Make a catalog from drizzle products using the SEP implementation of SourceExtractor
    
    Parameters
    ----------
    root : str
        Rootname of the FITS images to use for source extraction.  This
        function is designed to work with the single-image products from 
        `drizzlepac`, so the default data/science image is searched by 
        
        >>> drz_file = glob.glob(f'{root}_dr[zc]_sci.fits*')[0]
        
        Note that this will find and use gzipped versions of the images, 
        if necessary.
        
        The associated weight image filename is then assumed to be
        
        >>> weight_file = drz_file.replace('_sci.fits', '_wht.fits')
        >>> weight_file = weight_file.replace('_drz.fits', '_wht.fits')
        
    sci, wht : str
        Filenames to override `drz_file` and `weight_file` derived from the
        ``root`` parameter.

    threshold : float
        Detection threshold for `sep.extract`
    
    get_background : bool
        Compute the background with `sep.Background`
    
    bkg_only : bool
        If `True`, then just return the background data array and don't run
        the source detection
    
    bkg_params : dict
        Keyword arguments for `sep.Background`.  Note that this can include
        a separate optional keyword ``pixel_scale`` that indicates that the 
        background sizes `bw`, `bh` are set for a paraticular pixel size.  
        They will be scaled to the pixel dimensions of the target images using
        the pixel scale derived from the image WCS.
    
    verbose : bool
        Print status messages
            
    phot_apertures : str or array-like
        Photometric aperture *diameters*. If given as a string then assume
        units of pixels. If an array or list, can have units, e.g.,
        `astropy.units.arcsec`.

    aper_segmask : bool
        If true, then run SEP photometry with segmentation masking.  This
        requires the sep fork at https://github.com/gbrammer/sep.git, 
        or `sep >= 1.10.0`.
    
    rescale_weight : bool
        If true, then a scale factor is calculated from the ratio of the 
        weight image to the variance estimated by `sep.Background`.  
    
    err_scale : float
        Explicit value to use for the weight scaling, rather than calculating
        with `rescale_weight`.  Only used if ``err_scale > 0``
    
    use_bkg_err : bool
        If true, then use the full error array derived by `sep.Background`.
        This is turned off by default in order to preserve the pixel-to-pixel 
        variation in the drizzled weight maps.
        
    column_case : func
        Function to apply to the catalog column names.  E.g., the default 
        `str.upper` results in uppercase column names
    
    save_to_fits : bool
        Save catalog FITS file ``{root}.cat.fits``
    
    include_wcs_extension : bool
        An extension will be added to the FITS catalog with the detection 
        image WCS
    
    source_xy : (x, y) or (ra, dec) arrays
        Force extraction positions.  If the arrays have units, then pass them
        through the header WCS.  If no units, positions are *zero indexed*
        array coordinates.

        To run with segmentation masking (`1sep > 1.10``), also provide 
        `aseg` and `aseg_id` arrays with `source_xy`, like

            >>> source_xy = ra, dec, aseg, aseg_id
    
    compute_auto_quantities : bool
        Compute Kron/auto-like quantities with 
        `~grizli.prep.compute_SEP_auto_params`

    autoparams : list
        Parameters of Kron/AUTO calculations with 
        `~grizli.prep.compute_SEP_auto_params`.
        
    flux_radii : list
        Light fraction radii to compute with 
        `~grizli.prep.compute_SEP_auto_params`, e.g., ``[0.5]`` will calculate
        the half-light radius (``FLUX_RADIUS``)
    
    subpix : int
        Pixel oversampling
        
    mask_kron : bool
        Not used
        
    max_total_corr : float
        Not used
    
    detection_params : dict
        Parameters passed to `sep.extract`
    
    bkg_mask : array
        Additional mask to apply to `sep.Background` calculation
        
    pixel_scale : float
        Not used
    
    log : bool
        Send log message to `grizli.utils.LOGFILE`
    
    gain : float
        Gain value passed to `sep.sum_circle`
        
    extract_pixstack : int
        See `sep.set_extract_pixstack`
    
    sub_object_limit : int
        See `sep.set_sub_object_limit`
    
    exposure_footprints : list, None
        An optional list of objects that can be parsed with `sregion.SRegion`.  If 
        specified, add a column ``nexp`` to the catalog corresponding to the number
        of entries in the list that overlap with a particular source position
    
    Returns
    -------
    tab : `~astropy.table.Table`
        Source catalog
    
    
    """
    # if log:
    #     frame = inspect.currentframe()
    #     utils.log_function_arguments(utils.LOGFILE, frame,
    #                                 'prep.make_SEP_catalog', verbose=True)

    import copy
    import astropy.units as u

    if detection_params is None:
        detection_params = prep.SEP_DETECT_PARAMS
    if phot_apertures is None:
        phot_apertures = prep.SEXTRACTOR_PHOT_APERTURES

    sep.set_extract_pixstack(extract_pixstack)
    sep.set_sub_object_limit(sub_object_limit)

    if in_dir is not None: in_dir = Path(in_dir)
    if out_dir is not None: out_dir = Path(out_dir)

    if sci is not None:
        drz_file = sci
    else:
        if in_dir is not None and in_dir.is_dir():
            drz_file = [str(f) for f in in_dir.glob(f'{root}_dr[zc]_sci.fits*')][0]
        else:
            drz_file = glob.glob('{root}_dr[zc]_sci.fits')[0]

    with pf.open(drz_file) as im:
        # im = pf.open(drz_file)

        # Filter
        drz_filter = utils.parse_filter_from_header(im[0].header)
        if 'PHOTPLAM' in im[0].header:
            drz_photplam = im[0].header['PHOTPLAM']
        else:
            drz_photplam = None

        # Get AB zeropoint
        ZP = utils.calc_header_zeropoint(im, ext=0)

        # logstr = 'sep: Image AB zeropoint =  {0:.3f}'.format(ZP)
        # utils.log_comment(utils.LOGFILE, logstr, verbose=verbose, show_date=True)

        # Scale fluxes to mico-Jy
        uJy_to_dn = 1/(3631*1e6*10**(-0.4*ZP))

    if wht is not None:
        weight_file = wht
    else:
        weight_file = drz_file.replace('_sci.fits', '_wht.fits')
        weight_file = weight_file.replace('_drz.fits', '_wht.fits')

    if (weight_file == drz_file) | (not os.path.exists(weight_file)):
        WEIGHT_TYPE = "NONE"
        weight_file = None
    else:
        WEIGHT_TYPE = "MAP_WEIGHT"

    with pf.open(drz_file) as drz_im:
    # drz_im = pf.open(drz_file)
        data = drz_im[0].data.byteswap().newbyteorder()

        try:
            wcs = pywcs.WCS(drz_im[0].header)
            wcs_header = utils.to_header(wcs)
            pixel_scale = utils.get_wcs_pscale(wcs)  # arcsec
        except:
            wcs = None
            wcs_header = drz_im[0].header.copy()
            pixel_scale = np.sqrt(wcs_header['CD1_1']**2+wcs_header['CD1_2']**2)
            pixel_scale *= 3600.  # arcsec

        # Add some header keywords to the wcs header
        for k in ['EXPSTART', 'EXPEND', 'EXPTIME']:
            if k in drz_im[0].header:
                wcs_header[k] = drz_im[0].header[k]

        if isinstance(phot_apertures, str):
            apertures = np.cast[float](phot_apertures.replace(',', '').split())
        else:
            apertures = []
            for ap in phot_apertures:
                if hasattr(ap, 'unit'):
                    apertures.append(ap.to(u.arcsec).value/pixel_scale)
                else:
                    apertures.append(ap)

        # Do we need to compute the error from the wht image?
        need_err = (not use_bkg_err) | (not get_background)
        if (weight_file is not None) & need_err:
            wht_im = pf.open(weight_file)
            wht_data = wht_im[0].data.byteswap().newbyteorder()

            err = 1/np.sqrt(wht_data)
            del(wht_data)

            # True mask pixels are masked with sep
            mask = (~np.isfinite(err)) | (err == 0) | (~np.isfinite(data))
            err[mask] = 0

            wht_im.close()
            del(wht_im)

        else:
            # True mask pixels are masked with sep
            mask = (data == 0) | (~np.isfinite(data))
            err = None

    try:
        # drz_im.close()
        del(drz_im)
    except:
        pass

    data_mask = np.cast[data.dtype](mask)

    if get_background | (err_scale < 0) | (use_bkg_err):

        # Account for pixel scale in bkg_params
        bkg_input = {}
        if 'pixel_scale' in bkg_params:
            bkg_pscale = bkg_params['pixel_scale']
        else:
            bkg_pscale = pixel_scale

        for k in bkg_params:
            if k in ['pixel_scale']:
                continue

            if k in ['bw', 'bh']:
                bkg_input[k] = bkg_params[k]*bkg_pscale/pixel_scale
            else:
                bkg_input[k] = bkg_params[k]

        logstr = 'SEP: Get background {0}'.format(bkg_input)
        utils.log_comment(utils.LOGFILE, logstr, verbose=verbose,
                        show_date=True)

        if bkg_mask is not None:
            bkg = sep.Background(data, mask=mask | bkg_mask, **bkg_input)
        else:
            bkg = sep.Background(data, mask=mask, **bkg_input)

        bkg_data = bkg.back()
        if bkg_only:
            return bkg_data
        
        if get_background == 2:
            if in_dir is not None and in_dir.is_dir():
                bkg_file = str(in_dir / f'{root}_bkg.fits')
            else:
                bkg_file = '{root}_bkg.fits'
            if os.path.exists(bkg_file):
                logstr = 'SEP: use background file {0}'.format(bkg_file)
                utils.log_comment(utils.LOGFILE, logstr, verbose=verbose,
                                show_date=True)
                
                bkg_im = pf.open(bkg_file)
                bkg_data = bkg_im[0].data*1
        # else:        
        #     pf.writeto('{0}_bkg.fits'.format(root), data=bkg_data,
        #             header=wcs_header, overwrite=True)

        if (err is None) | use_bkg_err:
            logstr = 'sep: Use bkg.rms() for error array'
            utils.log_comment(utils.LOGFILE, logstr, verbose=verbose,
                            show_date=True)

            err = bkg.rms()

        if err_scale == -np.inf:
            ratio = bkg.rms()/err
            err_scale = np.median(ratio[(~mask) & np.isfinite(ratio)])
        else:
            # Just return the error scale
            if err_scale < 0:
                ratio = bkg.rms()/err
                xerr_scale = np.median(ratio[(~mask) & np.isfinite(ratio)])
                del(bkg)
                return xerr_scale

        del(bkg)

    else:
        if err_scale is None:
            err_scale = 1.

    if not get_background:
        bkg_data = 0.
        data_bkg = data
    else:
        data_bkg = data - bkg_data

    if rescale_weight:
        if verbose:
            print('SEP: err_scale={:.3f}'.format(err_scale))

        err *= err_scale

    if source_xy is None:
        # Run the detection
        if verbose:
            print('SEP: Extract...')

        objects, seg = sep.extract(
            data_bkg, 
            threshold, 
            err=err,
            mask=mask, 
            segmentation_map=seg_image if seg_image is not None else True,
            **detection_params,
        )

        if verbose:
            print("\033[1A", end="\x1b[2K")
            print ("SEP: Extract... Done.")

        tab = utils.GTable(objects)
        tab.meta['VERSION'] = (sep.__version__, 'SEP version')

        # make unit-indexed like SExtractor
        tab['x_image'] = tab['x']+1
        tab['y_image'] = tab['y']+1

        # ID
        tab['number'] = np.arange(len(tab), dtype=np.int32)+1
        if seg_image is not None:
            ids = np.unique(seg_image)
            tab['number'] = ids[ids>0]

        tab['theta'] = np.clip(tab['theta'], -np.pi/2, np.pi/2)
        for c in ['a', 'b', 'x', 'y', 'x_image', 'y_image', 'theta']:
            tab = tab[np.isfinite(tab[c])]

        if seg_image is None:
            seg_image = seg
            
        # Segmentation
        seg_image[mask] = 0

        try:
            _out_path = out_dir / f"{root}_seg.fits"
        except:
            _out_path = f"{root}_seg.fits"
        pf.writeto(
            _out_path,
            data=seg_image,
            header=wcs_header, 
            overwrite=True,
        )

        # WCS coordinates
        if wcs is not None:
            tab['ra'], tab['dec'] = wcs.all_pix2world(tab['x'], tab['y'], 0)
            tab['ra'].unit = u.deg
            tab['dec'].unit = u.deg
            tab['x_world'], tab['y_world'] = tab['ra'], tab['dec']

        if 'minarea' in detection_params:
            tab.meta['MINAREA'] = (detection_params['minarea'],
                                'Minimum source area in pixels')
        else:
            tab.meta['MINAREA'] = (5, 'Minimum source area in pixels')

        if 'clean' in detection_params:
            tab.meta['CLEAN'] = (detection_params['clean'],
                                'Detection cleaning')
        else:
            tab.meta['CLEAN'] = (True, 'Detection cleaning')

        if 'deblend_cont' in detection_params:
            tab.meta['DEBCONT'] = (detection_params['deblend_cont'],
                                'Deblending contrast ratio')
        else:
            tab.meta['DEBCONT'] = (0.005, 'Deblending contrast ratio')

        if 'deblend_nthresh' in detection_params:
            tab.meta['DEBTHRSH'] = (detection_params['deblend_nthresh'],
                                'Number of deblending thresholds')
        else:
            tab.meta['DEBTHRSH'] = (32, 'Number of deblending thresholds')

        if 'filter_type' in detection_params:
            tab.meta['FILTER_TYPE'] = (detection_params['filter_type'],
                                    'Type of filter applied, conv or weight')
        else:
            tab.meta['FILTER_TYPE'] = ('conv',
                                    'Type of filter applied, conv or weight')

        tab.meta['THRESHOLD'] = (threshold, 'Detection threshold')

        # ISO fluxes (flux within segments)
        iso_flux, iso_fluxerr, iso_area = prep.get_seg_iso_flux(data_bkg, seg, tab,
                                                err=err, verbose=1)

        tab['flux_iso'] = iso_flux/uJy_to_dn*u.uJy
        tab['fluxerr_iso'] = iso_fluxerr/uJy_to_dn*u.uJy
        tab['area_iso'] = iso_area
        tab['mag_iso'] = 23.9-2.5*np.log10(tab['flux_iso'])

        # Compute FLUX_AUTO, FLUX_RADIUS
        if compute_auto_quantities:
            auto = prep.compute_SEP_auto_params(data, data_bkg, mask,
                                pixel_scale=pixel_scale,
                                err=err, segmap=seg, tab=tab,
                                autoparams=autoparams, flux_radii=flux_radii,
                                subpix=subpix, verbose=verbose)

            for k in auto.meta:
                tab.meta[k] = auto.meta[k]

            auto_flux_cols = ['flux_auto', 'fluxerr_auto', 'bkg_auto']
            for c in auto.colnames:
                if c in auto_flux_cols:
                    tab[c] = auto[c]/uJy_to_dn*u.uJy
                else:
                    tab[c] = auto[c]

            # Problematic sources
            # bad = (tab['flux_auto'] <= 0) | (tab['flux_radius'] <= 0)
            # bad |= (tab['kron_radius'] <= 0)
            # tab = tab[~bad]

            # Correction for flux outside Kron aperture
            tot_corr = prep.get_kron_tot_corr(tab, drz_filter, 
                                                pixel_scale=pixel_scale, 
                                                photplam=drz_photplam)

            tab['tot_corr'] = tot_corr
            tab.meta['TOTCFILT'] = (drz_filter, 'Filter for tot_corr')
            tab.meta['TOTCWAVE'] = (drz_photplam, 'PLAM for tot_corr')

            total_flux = tab['flux_auto']*tot_corr
            tab['mag_auto'] = 23.9-2.5*np.log10(total_flux)
            tab['magerr_auto'] = 2.5/np.log(10)*(tab['fluxerr_auto']/tab['flux_auto'])

        # More flux columns
        for c in ['cflux', 'flux', 'peak', 'cpeak']:
            tab[c] *= 1. / uJy_to_dn
            tab[c].unit = u.uJy

        source_x, source_y = tab['x'], tab['y']

        # Use segmentation image to mask aperture fluxes
        if aper_segmask:
            aseg = seg
            aseg_id = tab['number']
        else:
            aseg = aseg_id = None

        # Rename some columns to look like SExtractor
        for c in ['a', 'b', 'theta', 'cxx', 'cxy', 'cyy', 'x2', 'y2', 'xy']:
            tab.rename_column(c, c+'_image')

    else:
        if len(source_xy) == 2:
            source_x, source_y = source_xy
            aseg, aseg_id = None, None
            aper_segmask = False
        else:
            source_x, source_y, aseg, aseg_id = source_xy
            aper_segmask = True

        if hasattr(source_x, 'unit'):
            if source_x.unit == u.deg:
                # Input positions are ra/dec, convert with WCS
                ra, dec = source_x, source_y
                source_x, source_y = wcs.all_world2pix(ra, dec, 0)

        tab = utils.GTable()
        tab.meta['VERSION'] = (sep.__version__, 'SEP version')
    
    # Exposure footprints
    #--------------------
    if (exposure_footprints is not None) & ('ra' in tab.colnames):        
        tab['nexp'] =  prep.catalog_exposure_overlaps(tab['ra'], tab['dec'],
                                        exposure_footprints=exposure_footprints)
                                        
        tab['nexp'].description = 'Number of overlapping exposures'
    
    # Info
    tab.meta['ZP'] = (ZP, 'AB zeropoint')
    if 'PHOTPLAM' in im[0].header:
        tab.meta['PLAM'] = (im[0].header['PHOTPLAM'], 'Filter pivot wave')
        if 'PHOTFNU' in im[0].header:
            tab.meta['FNU'] = (im[0].header['PHOTFNU'], 'Scale to Jy')

        tab.meta['FLAM'] = (im[0].header['PHOTFLAM'], 'Scale to flam')

    tab.meta['uJy2dn'] = (uJy_to_dn, 'Convert uJy fluxes to image DN')

    tab.meta['DRZ_FILE'] = (drz_file[:36], 'SCI file')
    tab.meta['WHT_FILE'] = (weight_file[:36], 'WHT file')

    tab.meta['GET_BACK'] = (get_background, 'Background computed')
    for k in bkg_params:
        tab.meta[f'BACK_{k.upper()}'] = (bkg_params[k],
                                        f'Background param {k}')
        
    tab.meta['ERR_SCALE'] = (err_scale, 'Scale factor applied to weight image (like MAP_WEIGHT)')
    tab.meta['RESCALEW'] = (rescale_weight, 'Was the weight applied?')

    tab.meta['APERMASK'] = (aper_segmask, 'Mask apertures with seg image')

    # Photometry
    for iap, aper in enumerate(apertures):
        if sep.__version__ > '1.03':
            # Should work with the sep fork at gbrammer/sep and latest sep
            flux, fluxerr, flag = sep.sum_circle(data_bkg,
                                                source_x, source_y,
                                                aper/2, err=err,
                                                gain=gain, subpix=subpix,
                                                segmap=aseg, seg_id=aseg_id,
                                                mask=mask)
        else:
            tab.meta['APERMASK'] = (False, 'Mask apertures with seg image - Failed')
            flux, fluxerr, flag = sep.sum_circle(data_bkg,
                                                source_x, source_y,
                                                aper/2, err=err,
                                                gain=gain, subpix=subpix,
                                                mask=mask)

        tab.meta['GAIN'] = gain

        tab['flux_aper_{0}'.format(iap)] = flux/uJy_to_dn*u.uJy
        tab['fluxerr_aper_{0}'.format(iap)] = fluxerr/uJy_to_dn*u.uJy
        tab['flag_aper_{0}'.format(iap)] = flag

        if get_background:
            try:
                flux, fluxerr, flag = sep.sum_circle(bkg_data,
                                        source_x, source_y,
                                        aper/2, err=None, gain=1.0,
                                        segmap=aseg, seg_id=aseg_id,
                                        mask=mask)
            except:
                flux, fluxerr, flag = sep.sum_circle(bkg_data,
                                        source_x, source_y,
                                        aper/2, err=None, gain=1.0,
                                        mask=mask)

            tab['bkg_aper_{0}'.format(iap)] = flux/uJy_to_dn*u.uJy
        else:
            tab['bkg_aper_{0}'.format(iap)] = 0.*u.uJy

        # Count masked pixels in the aperture, not including segmask
        flux, fluxerr, flag = sep.sum_circle(data_mask,
                                    source_x, source_y,
                                    aper/2, err=err,
                                    gain=gain, subpix=subpix)

        tab['mask_aper_{0}'.format(iap)] = flux

        tab.meta['aper_{0}'.format(iap)] = (aper, 'Aperture diameter, pix')
        tab.meta['asec_{0}'.format(iap)] = (aper*pixel_scale,
                                            'Aperture diameter, arcsec')

    # # If blended, use largest aperture magnitude
    # if 'flag' in tab.colnames:
    #     last_flux = tab['flux_aper_{0}'.format(iap)]
    #     last_fluxerr = tab['fluxerr_aper_{0}'.format(iap)]
    #
    #     blended = (tab['flag'] & 1) > 0
    #
    #     total_corr = tab['flux_auto']/last_flux
    #     blended |= total_corr > max_total_corr
    #
    #     tab['flag'][blended] |= 1024
    #
    #     aper_mag = 23.9 - 2.5*np.log10(last_flux)
    #     aper_magerr = 2.5/np.log(10)*last_fluxerr/last_flux
    #
    #     tab['mag_auto'][blended] = aper_mag[blended]
    #     tab['magerr_auto'][blended] = aper_magerr[blended]
    #
    #     # "ISO" mag, integrated within the segment
    #     tab['mag_iso'] = 23.9-2.5*np.log10(tab['flux'])

    try:
        # Free memory objects explicitly
        del(data_mask)
        del(data)
        del(err)
    except:
        pass

    # if uppercase_columns:
    for c in tab.colnames:
        tab.rename_column(c, column_case(c))

    if save_to_fits:
        try:
            _out_path = out_dir / f"{root}.cat.fits"
        except:
            _out_path = f"{root}.cat.fits"
        tab.write(
            _out_path,
            format='fits', 
            overwrite=True,
        )

        if include_wcs_extension:
            try:
                hdul = pf.open(
                    _out_path,
                    mode='update',
                )
                wcs_hdu = pf.ImageHDU(header=wcs_header, data=None,
                                        name='WCS')
                hdul.append(wcs_hdu)
                hdul.flush()
            except:
                pass

    # logstr = '# SEP {0}.cat.fits: {1:d} objects'.format(root, len(tab))
    # utils.log_comment(utils.LOGFILE, logstr, verbose=verbose)

    return tab

def regen_multiband_catalogue(
    field_root,
    threshold=1.8, 
    detection_background=True, 
    photometry_background=True, 
    get_all_filters=False, 
    filters=None, 
    det_err_scale=-np.inf, 
    phot_err_scale=-np.inf, 
    rescale_weight=True, 
    run_detection=True, 
    detection_filter='ir', 
    detection_root=None, 
    output_root=None,
    use_psf_filter=True, 
    detection_params=None,  
    phot_apertures=None, 
    master_catalog=None,
    bkg_mask=None, 
    bkg_params={'bw': 64, 'bh': 64, 'fw': 3, 'fh': 3, 'pixel_scale': 0.06}, 
    use_bkg_err=False,
    aper_segmask=True, 
    sci_image=None, 
    clean_bkg=True,
    seg_image=None,
    in_dir=None,
    out_dir=None,
):
    if in_dir is not None: in_dir = Path(in_dir)
    if detection_params is None:
        detection_params=prep.SEP_DETECT_PARAMS
    if phot_apertures is None:
        phot_apertures=prep.SEXTRACTOR_PHOT_APERTURES_ARCSEC

    if detection_root is None:
        detection_root = '{0}-{1}'.format(field_root, detection_filter)

    if output_root is None:
        output_root = field_root

    if use_psf_filter:
        # psf_files = glob.glob('{0}*psf.fits'.format(self.field_root))
        if in_dir is not None and in_dir.is_dir():
            psf_files = [str(f) for f in in_dir.glob(f'{field_root}*psf.fits*')]
        else:
            psf_files = glob.glob('{0}*psf.fits'.format(field_root))

        if len(psf_files) > 0:
            psf_files.sort()
            psf_im = pf.open(psf_files[-1])

            msg = '# Generate PSF kernel from {0}\n'.format(psf_files[-1])
            utils.log_comment(utils.LOGFILE, msg, verbose=True)

            sh = psf_im['PSF', 'DRIZ1'].data.shape
            # Cut out center of PSF
            skip = (sh[0]-1-11)//2
            psf = psf_im['PSF', 'DRIZ1'].data[skip:-1-skip, skip:-1-skip]*1

            # Optimal filter is reversed PSF (i.e., PSF cross-correlation)
            # https://arxiv.org/pdf/1512.06872.pdf
            psf_kernel = psf[::-1, :][:, ::-1]
            psf_kernel /= psf_kernel.sum()

            detection_params['filter_kernel'] = psf_kernel

    tab = make_SEP_catalog(
        root=detection_root,
        sci=sci_image, 
        threshold=threshold,
        get_background=detection_background,
        save_to_fits=True,
        rescale_weight=rescale_weight,
        err_scale=det_err_scale,
        phot_apertures=phot_apertures,
        detection_params=detection_params,
        bkg_mask=bkg_mask,
        bkg_params=bkg_params,
        use_bkg_err=use_bkg_err,
        aper_segmask=aper_segmask,
        seg_image=seg_image,
        in_dir=in_dir,
        out_dir=out_dir,
    )

    cat_pixel_scale = tab.meta['asec_0'][0]/tab.meta['aper_0'][0]

    if aper_segmask:
        seg_data = pf.open(out_dir / f'{detection_root}_seg.fits')[0].data
        seg_data = np.cast[np.int32](seg_data)

        aseg, aseg_id = seg_data, tab['NUMBER']

        source_xy = tab['X_WORLD'], tab['Y_WORLD'], aseg, aseg_id
        aseg_half = None
    else:
        source_xy = tab['X_WORLD'], tab['Y_WORLD']

    if filters is None:
        #visits_file = '{0}_visits.yaml'.format(field_root)
        visits_file = auto_script.find_visit_file(root=field_root, path=str(in_dir))
        if visits_file is None:
            get_all_filters = True

        if get_all_filters:
            mq = '{0}-f*dr?_sci.fits*'
            mq = mq.format(field_root.replace('-100mas','-*mas'))
            # mosaic_files = glob.glob(mq)
            if in_dir is not None and in_dir.is_dir():
                mosaic_files = [str(f) for f in in_dir.glob(mq)]
            else:
                mosaic_files = glob.glob(mq)
            
            mq = '{0}-clear*dr?_sci.fits*'
            mq = mq.format(field_root.replace('-100mas','-*mas'))
            # mosaic_files += glob.glob(mq)
            if in_dir is not None and in_dir.is_dir():
                mosaic_files += [str(f) for f in in_dir.glob(mq)]
            else:
                mosaic_files += glob.glob(mq)
            
            mosaic_files.sort()
            
            filters = [Path(file).stem.split('_')[-3][len(field_root)+1:] 
                       for file in mosaic_files]
        else:
            #vfile = '{0}_visits.npy'.format(field_root)
            #visits, all_groups, info = np.load(vfile, allow_pickle=True)
            visits, all_groups, info = auto_script.load_visit_info(
                field_root, 
                path=in_dir,
                verbose=False,
            )

            # UVIS
            info_filters = [f for f in info['FILTER']]
            for i in range(len(info)):
                file_i = info['FILE'][i]
                if file_i.startswith('i') & ('_flc' in file_i):
                    info_filters[i] += 'U'

            info['FILTER'] = info_filters

            filters = [f.lower() for f in np.unique(info['FILTER'])]

    fq = '{0}-{1}_dr?_sci.fits*'

    utils.log_comment(
        utils.LOGFILE, 
        f"Filters: {filters}",
        verbose=True,
    ) 
    
    for ii, filt in enumerate(filters):
        utils.log_comment(
            utils.LOGFILE,
            f"Trying {filt} filter...",
            verbose=True,
        )
        if filt.startswith('g'):
            continue

        if filt not in ['g102', 'g141', 'g800l']:
            _fstr = fq.format(field_root.replace('-100mas','-*mas'), filt)
            # sci_files = glob.glob(_fstr)
            if in_dir is not None and in_dir.is_dir():
                sci_files = [str(f) for f in in_dir.glob(_fstr)]
            else:
                sci_files = glob.glob(_fstr)

            if len(sci_files) == 0:
                continue

            root = sci_files[0].split('{0}_dr'.format(filt))[0]+filt
            # root = '{0}-{1}'.format(field_root, filt)

            # Check for half-pixel optical images if using segmask
            if aper_segmask:
                sci = pf.open(sci_files[0])
                sci_shape = sci[0].data.shape
                sci.close()
                del(sci)

                if sci_shape[0] != aseg.shape[0]:
                    msg = '# filt={0}, need half-size segmentation image'
                    msg += ', shapes sci:{1} seg:{2}'
                    print(msg.format(filt, sci_shape, aseg.shape))
                    
                    if aseg_half is None:
                        aseg_half = np.zeros(sci_shape, dtype=aseg.dtype)
                        for i in [0, 1]:
                            for j in [0, 1]:
                                aseg_half[i::2, j::2] += aseg

                    source_xy = (tab['X_WORLD'], tab['Y_WORLD'],
                                 aseg_half, aseg_id)
                else:
                    source_xy = (tab['X_WORLD'], tab['Y_WORLD'],
                                 aseg, aseg_id)

            root = Path(root).relative_to(in_dir) if in_dir is not None else root
            filter_tab = make_SEP_catalog(
                root=root,
                threshold=threshold,
                rescale_weight=rescale_weight,
                err_scale=phot_err_scale,
                get_background=photometry_background,
                save_to_fits=False,
                source_xy=source_xy,
                phot_apertures=phot_apertures,
                bkg_mask=bkg_mask,
                bkg_params=bkg_params,
                use_bkg_err=use_bkg_err,
                sci=sci_image,
                seg_image=seg_image,
                in_dir=in_dir,
                out_dir=out_dir,
            )

            for k in filter_tab.meta:
                
                newk = '{0}_{1}'.format(filt.upper(), k)
                newk = newk.replace('-CLEAR','')
                tab.meta[newk] = filter_tab.meta[k]

            for c in filter_tab.colnames:
                newc = '{0}_{1}'.format(filt.upper(), c)
                tab[newc] = filter_tab[c]

            # Kron total correction from EE
            newk = '{0}_PLAM'.format(filt.upper())
            newk = newk.replace('-CLEAR','')
            # filt_plam = tab.meta[newk]

            tot_corr = prep.get_kron_tot_corr(tab, filt.lower(), 
                                                pixel_scale=cat_pixel_scale, 
                                                # photplam=filt_plam,
                                            )

            #ee_corr = prep.get_kron_tot_corr(tab, filter=filt.lower())
            tab['{0}_tot_corr'.format(filt.upper())] = tot_corr
            
            if clean_bkg:
                # bkg_files = glob.glob(f'{root}*{filt}*bkg.fits')

                if in_dir is not None and in_dir.is_dir():
                    bkg_files = [str(f) for f in in_dir.glob(f'{root}*{filt}*bkg.fits')]
                else:
                    bkg_files = glob.glob(f'{root}*{filt}*bkg.fits')

                for bfile in bkg_files:
                    print('# rm {bfile}')
                    os.remove(bfile)
        
        else:
            continue

    for c in tab.colnames:
        tab.rename_column(c, c.lower())

    idcol = utils.GTable.Column(data=tab['number'], name='id')
    tab.add_column(idcol, index=0)

    try:
        out_path = out_dir / f"{output_root}_phot.fits"
    except:
        out_path = f"{root}_phot.fits"
    tab.write(
        out_path,
        format='fits',
        overwrite=True,
    )

    return tab
    