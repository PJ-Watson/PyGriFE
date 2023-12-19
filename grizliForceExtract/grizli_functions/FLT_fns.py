import os
from grizli import model

def get_beams_with_spectrum(grp, id, size=10, center_rd=None, beam_id='A',
                min_overlap=0.1, min_valid_pix=10, min_mask=0.01,
                min_sens=0.08, mask_resid=True, get_slice_header=True,
                show_exception=False, spectrum1d=None):
    """Extract 2D spectra "beams" from the GroupFLT exposures.

    Parameters
    ----------
    id : int
        Catalog ID of the object to extract.

    size : int
        Half-size of the 2D spectrum to extract, along cross-dispersion
        axis.

    center_rd : optional, (float, float)
        Extract based on RA/Dec rather than catalog ID.

    beam_id : type
        Name of the order to extract.

    min_overlap : float
        Fraction of the spectrum along wavelength axis that has one
        or more valid pixels.

    min_valid_pix : int
        Minimum number of valid pixels (`beam.fit_mask == True`) in 2D
        spectrum.

    min_mask : float
        Minimum factor relative to the maximum pixel value of the flat
        f-lambda model where the 2D cutout data are considered good.
        Passed through to `~grizli.model.BeamCutout`.

    min_sens : float
        See `~grizli.model.BeamCutout`.

    get_slice_header : bool
        Passed to `~grizli.model.BeamCutout`.

    Returns
    -------
    beams : list
        List of `~grizli.model.BeamCutout` objects.

    """
    beams = grp.compute_single_model(
        id,
        center_rd=center_rd,
        size=size,
        store=False,
        get_beams=[beam_id],
        spectrum1d=spectrum1d,
    )

    out_beams = []
    for flt, beam in zip(grp.FLTs, beams):
        try:
            out_beam = model.BeamCutout(flt=flt, beam=beam[beam_id],
                                    conf=flt.conf, min_mask=min_mask,
                                    min_sens=min_sens,
                                    mask_resid=mask_resid,
                                    get_slice_header=get_slice_header)
        except:
            #print('Except: get_beams')
            if show_exception:
                utils.log_exception(utils.LOGFILE, traceback)
                
            continue

        valid = (out_beam.grism['SCI'] != 0)
        valid &= out_beam.fit_mask.reshape(out_beam.sh)
        hasdata = (valid.sum(axis=0) > 0).sum()
        if hasdata*1./out_beam.model.shape[1] < min_overlap:
            continue

        # Empty direct image?
        if out_beam.beam.total_flux == 0:
            continue

        if out_beam.fit_mask.sum() < min_valid_pix:
            continue

        out_beams.append(out_beam)

    return out_beams


def load_and_mod_FLT(grism_file, sci_extn, direct_file, pad, ref_file,
               ref_ext, seg_file, verbose, catalog, ix, use_jwst_crds):
    """Helper function for loading `.model.GrismFLT` objects with `multiprocessing`.

    TBD
    """
    import time
    try:
        import cPickle as pickle
    except:
        # Python 3
        import pickle
    
    # slight random delay to avoid synchronization problems
    # np.random.seed(ix)
    # sleeptime = ix*1
    # print '%s sleep %.3f %d' %(grism_file, sleeptime, ix)
    # time.sleep(sleeptime)

    # print grism_file, direct_file

    new_root = '.{0:02d}.GrismFLT.fits'.format(sci_extn)
    save_file = grism_file.replace('_flt.fits', new_root)
    save_file = save_file.replace('_flc.fits', new_root)
    save_file = save_file.replace('_cmb.fits', new_root)
    save_file = save_file.replace('_rate.fits', new_root)
    save_file = save_file.replace('_elec.fits', new_root)
    
    if (save_file == grism_file) & ('GrismFLT' not in grism_file):
        # couldn't build new filename based on the extensions
        # so just insert at the end
        save_file = grism_file.replace('.fits', new_root)
        
    if (grism_file.find('_') < 0) & ('GrismFLT' not in grism_file):
        save_file = 'xxxxxxxxxxxxxxxxxxx'

    if os.path.exists(save_file) & ('GrismFLT' in save_file):
        print('Load {0}!'.format(save_file))

        fp = open(save_file.replace('GrismFLT.fits', 'GrismFLT.pkl'), 'rb')
        flt = pickle.load(fp)
        fp.close()

        flt.conf_file = f"{os.environ['GRIZLI']}/CONF/{flt.conf_file.split('/')[-1]}"

        status = flt.load_from_fits(save_file)

        if seg_file is not None:
            flt.process_seg_file(seg_file)

    else:        
        flt = model.GrismFLT(grism_file=grism_file, sci_extn=sci_extn,
                         direct_file=direct_file, pad=pad,
                         ref_file=ref_file, ref_ext=ref_ext,
                         seg_file=seg_file, shrink_segimage=True,
                         verbose=verbose, use_jwst_crds=use_jwst_crds)

    if flt.direct.wcs.wcs.has_pc():
        for obj in [flt.grism, flt.direct]:
            obj.get_wcs()

    if catalog is not None:
        flt.catalog = flt.blot_catalog(catalog,
                                   sextractor=('X_WORLD' in catalog.colnames))
        flt.catalog_file = catalog

    else:
        flt.catalog = None
    
    if flt.grism.instrument in ['NIRCAM']:
        flt.apply_POM()

    if flt.grism.instrument in ['NIRISS', 'NIRCAM']:
        flt.transform_JWST_WFSS()
    
    if hasattr(flt, 'conf'):
        delattr(flt, 'conf')
    
    return flt  # , out_cat