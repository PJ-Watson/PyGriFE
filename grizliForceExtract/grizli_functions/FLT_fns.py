import os
from grizli import model, utils
import numpy as np

def get_beams_with_spectrum(grp, id, size=10, center_rd=None, beam_id='A',
                min_overlap=0.1, min_valid_pix=10, min_mask=0.01,
                min_sens=0.08, mask_resid=True, get_slice_header=True,
                show_exception=False, spectrum_1d=None, spectrum_is_cgs=False):
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
        # spectrum_1d=spectrum_1d,
        # is_cgs=True,
    )
    # print (beams)
    # print (dir(beams))
    # print (beams[0])
    # print (dir(beams[0]))
    # print (beams[0]["A"])
    # print (dir(beams[0]["A"]))
    # print ("1d_spec", beams[0]["A"].spectrum_1d)

    out_beams = []
    for flt, beam in zip(grp.FLTs, beams):
        # print (dir(flt))
        # for k in dir(flt):
        #     print (k, flt.k)
        print (flt.transform_JWST_WFSS)
        # print (flt[0])
        try:
            beam_in = beam[beam_id]
            beam_in.spectrum_1d = spectrum_1d
            beam_in.is_cgs = spectrum_is_cgs

            # Is this a modified segmentation map?
            if (spectrum_1d is not None) and hasattr(flt, "orig_seg"):

                # Find the most common id in the relevant area of the original map
                u, c = np.unique(flt.orig_seg[flt.orig_seg==id], return_counts=True)
                orig_id = u[c.argmax()]

                # scale the 1d spectrum by the fraction of the original flux 
                # contained in the new segmentation map
                # CAUTION - This assumes the new region is a subset of the original
                # print (beam_in.direct)
                # print ("FLT.DIRECT", dir(flt.direct))
                # print ("FLT", dir(flt))
                # print (flt.ref_file)
                # print (flt.direct.data["REF"].shape)
                # print (np.nansum(flt.direct.data["REF"]))
                # # print (flt.direct.data)
                # # print (flt.direct.data.shape)
                # print (beam_in.direct.shape)
                # print (flt.seg.shape)
                # print (beam_in.seg.shape)

                print (f'ORIG_ID: {orig_id}, OLD_FLUX: {flt.direct.data["REF"][flt.orig_seg == orig_id].sum()}, NEW_FLUX: {flt.direct.data["REF"][flt.seg == id].sum()} ')
                print ("SCALE FACTOR:", flt.direct.data["REF"][flt.seg == id].sum() / flt.direct.data["REF"][flt.orig_seg == orig_id].sum())
                spectrum_1d[1] = np.asarray(spectrum_1d[1])*flt.direct.data["REF"][flt.seg == id].sum() / flt.direct.data["REF"][flt.orig_seg == orig_id].sum()
                # print (spectrum_1d[1])
                # print (flt[0])
                # import matplotlib.pyplot as plt
                # fig, ax = plt.subplots(2,1)
                # ax[0].imshow(flt.orig_seg)
                # ax[1].imshow(flt.seg, cmap="plasma")
                # plt.show()

            out_beam = model.BeamCutout(flt=flt, beam=beam_in,
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
            if hasattr(flt, "seg"):
                flt.orig_seg = np.rot90(flt.seg, -1)
            if hasattr(flt, "seg_file"):
                flt.orig_seg_file = flt.seg_file
            flt.process_seg_file(seg_file)

        # print (dir(flt.object_dispersers))
        # print (dir(flt.object_dispersers[4149]))
        # for k, v in flt.object_dispersers[4149].items():
        #     print (f"key: {k}, value: {v}")
        # print (flt.object_dispersers["test"])
        # print (flt.test[0])
        # exit()

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
        
        # import matplotlib.pyplot as plt
        # fig, ax = plt.subplots(2,1)
        # ax[0].imshow(flt.orig_seg)
        # ax[1].imshow(flt.seg, cmap="plasma")
        # plt.show()

        if hasattr(flt, "orig_seg"):
            rot = 0
            if flt.grism.instrument == 'NIRISS':
                if flt.grism.filter == 'GR150C':
                    rot = 2
                else:
                    rot = -1

            elif flt.grism.instrument in ['NIRCAM', 'NIRCAMA']:
                if flt.grism.module == 'A':
                    #  Module A
                    if flt.grism.pupil == 'GRISMC':
                        rot = 1
                else:
                    # Module B
                    if flt.grism.pupil == 'GRISMC':
                        rot = 1
                    else:
                        rot = 2
                
            elif flt.grism.instrument == 'NIRCAMB':
                if flt.grism.pupil == 'GRISMC':
                    rot = 1
                else:
                    rot = 2

            if flt.is_rotated:
                rot *= -1

            is_rotated = not flt.is_rotated
            print('SECONDARY: Transform JWST WFSS: flip={0}'.format(is_rotated))

            rot += 1

            # Rotate segmentation image
            flt.orig_seg = np.rot90(flt.orig_seg, rot)

        flt.transform_JWST_WFSS()


        # import matplotlib.pyplot as plt
        # fig, ax = plt.subplots(2,1)
        # ax[0].imshow(flt.orig_seg)
        # ax[1].imshow(flt.seg, cmap="plasma")
        # plt.show()
    
    if hasattr(flt, 'conf'):
        delattr(flt, 'conf')
    
    return flt  # , out_cat