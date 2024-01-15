
import os
import glob
import numpy as np
import astropy.io.fits as pf
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from importlib.metadata import version
import astropy.units as u
from astropy.coordinates import SkyCoord 
from tqdm import tqdm
from pathlib import Path
import pickle
import sep
# from grizliForceExtract.utils_c import sep_fns
import math
import time
import multiprocessing
from functools import partial
import shutil
import logging

if version("sep").split(".")[1] < "3":
    raise ImportError("""
        Forced extraction requires SEP>=1.3.0. Please install the fork maintained
        at https://github.com/PJ-Watson/sep/tree/dual-image.
    """)

# grizli_dir = Path(
#     "/Path/to/out/dir" # Point this to the output directory
# ) / "grizli"

# grizli_dir.mkdir(exist_ok=True)
# (grizli_dir / "CONF").mkdir(exist_ok=True)
# (grizli_dir / "templates").mkdir(exist_ok=True)
# (grizli_dir / "iref").mkdir(exist_ok=True)
# (grizli_dir / "iref").mkdir(exist_ok=True)
# os.environ["GRIZLI"] = str(grizli_dir)
# os.environ["iref"] = str(grizli_dir / "iref")
# os.environ["jref"] = str(grizli_dir / "jref")
# cwd = Path.cwd()
# os.chdir(os.path.join(grizli.GRIZLI_PATH, 'CONF'))

# grizli.utils.fetch_default_calibs()
# grizli.utils.fetch_config_files(get_jwst=True)
# if not os.path.exists('GR150C.F115W.221215.conf'):
#     os.system('wget "https://zenodo.org/record/7628094/files/niriss_config_221215.tar.gz?download=1" -O niriss_config_221215.tar.gz')
#     os.system('tar xzvf niriss_config_221215.tar.gz')
# if not os.path.exists('niriss_sens_221215.tar.gz'):
#     os.system('wget "https://zenodo.org/record/7628094/files/niriss_sens_221215.tar.gz" -O niriss_sens_221215.tar.gz')
#     os.system('tar xzvf niriss_sens_221215.tar.gz')
# os.chdir(cwd)
# eazy.fetch_eazy_photoz()

# try:
#     grizli.utils.symlink_templates(force=True)
# except:
#     pass

try:
    Path(os.getenv("GRIZLI")).is_dir()
    Path(os.getenv("iref")).is_dir()
    Path(os.getenv("jref")).is_dir()
except:
    raise RuntimeError("""
        Either the grizli environment variables are not set correctly, 
        or the directories they point to do not yet exist. 
        Check that the environment is correctly configured, or edit and 
        uncomment the lines above.
    """)

import grizli
from grizli import utils, prep, jwst_utils, multifit, fitting
from grizli.pipeline import auto_script, photoz
import eazy

from grizli_functions import catalogue_fns, FLT_fns

extended_translate = {
    'f098m': 201, 'f105w': 202, 'f110w': 241, 'f125w': 203, 'f140w': 204, 
    'f160w': 205, 'f435w': 233, 'f475w': 234, 'f555w': 235, 'f606w': 236, 
    'f625w': 237, 'f775w': 238, 'f814w': 239, 'f850lp': 240, 'f702w': 15, 
    'f600lpu': 243, 'f225wu': 207, 'f275wu': 208, 'f336wu': 209, 'f350lpu': 339, 
    'f438wu': 211, 'f475wu': 212, 'f475xu': 242, 'f555wu': 213, 'f606wu': 214, 
    'f625wu': 215, 'f775wu': 216, 'f814wu': 217, 'f390wu': 210, 'ch1': 18, 
    'ch2': 19, 'f336w':209, 'f350lp':339, 'f115w': 309, 'f150w': 310, 'f200w': 311,
    "f277w": 375,
}

multifit._loadFLT = FLT_fns.load_and_mod_FLT

# photoz.FILTER_TRANS = extended_translate
# photoz.eazy_photoz = partial(photoz.eazy_photoz, filter_trans=extended_translate)

# print (photoz.eazy_photoz.__defaults__)

"""
Steps:
 - Recalculate photometry for affected objects - look at this function:

        auto_script.multiband_catalog(field_root=root,
            **multiband_catalog_args)

 - Edit contamination maps - need to modify `multifit.GroupFLT`

 - Permanently modify *GrismFLT.fits, or re-run modifications each time?

 - Force all functions here to save to a different folder, e.g. ForcedExtractions, and save all files there.

 - ImageViewer module can be as generic as possible, but `Extract Object` button can be customised to work only with grizli.

 - Add redshift range to search for (min/max entries).

 - [DONE] Sub-class imageviewer? Make one completely generic function, and one specifically for grizli.

 - Separate out installation (take inspiration from grizli). Make installing grizli components optional.

"""

# def _modify_conf_file_location(file):

#     if os.path.exists(file) & ('GrismFLT' in file):
#         print(f'Checking conf file location for {file}.')

#         with open(file.replace('GrismFLT.fits', 'GrismFLT.pkl'), 'rb') as fp:
#             flt = pickle.load(fp)

#         flt.conf_file = f"{os.environ['GRIZLI']}/CONF/{flt.conf_file.split('/')[-1]}"

#         with open(file.replace('GrismFLT.fits', 'GrismFLT.pkl'), 'wb') as fp:
#             pickle.dump(flt, fp)

class GrizliExtractor:

    def __init__(self, field_root, in_dir, out_dir):

        self.field_root = field_root
        self.in_dir = Path(in_dir)
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(exist_ok=True, parents=True)

        copy_patterns = ["*FLT.fits", "*FLT.pkl", "*01.wcs.fits"]
        for pattern in copy_patterns:
            orig_files = self.in_dir.glob(pattern)
            for o in orig_files:
                out_path = self.out_dir/o.name
                if not out_path.is_file():
                    utils.log_comment(
                        utils.LOGFILE,
                        f"Copying file {o}",
                        verbose=True,
                    )
                    shutil.copy(
                        src=o,
                        dst=self.out_dir/o.name,
                        follow_symlinks=True,
                    )
        
        link_patterns = ["*drz_sci.fits", "*drz_wht.fits"]
        for pattern in link_patterns:
            orig_files = self.in_dir.glob(self.field_root+pattern)
            for o in orig_files:
                out_path = self.out_dir/o.name
                try:
                    out_path.symlink_to(o)
                except:
                    utils.log_comment(
                        utils.LOGFILE,
                        f"File {out_path.name} exists already.",
                        verbose=True,
                    )
                # print (o.name)

        # print ([*flt_files])

        # grizli_dir = self.out_dir / "grizli"
        # grizli_dir.mkdir(exist_ok=True)
        # (grizli_dir / "CONF").mkdir(exist_ok=True)
        # (grizli_dir / "templates").mkdir(exist_ok=True)
        # (grizli_dir / "iref").mkdir(exist_ok=True)
        # (grizli_dir / "iref").mkdir(exist_ok=True)
        # os.environ["GRIZLI"] = str(grizli_dir)
        # os.environ["iref"] = str(grizli_dir / "iref")
        # os.environ["jref"] = str(grizli_dir / "jref")

        # print (os.environ["GRIZLI"])

        # print (grizli.GRIZLI_PATH)

        # grizli.GRIZLI_PATH = os.environ["GRIZLI"]

        # reload (grizli)
        # print (grizli.GRIZLI_PATH)

    def load_orig_seg_img(self, seg_img_path, ext=0):

        seg_img_path = Path(seg_img_path)

        self.seg_name = seg_img_path.name

        shutil.copy(
            src=seg_img_path,
            dst=self.out_dir/seg_img_path.with_stem(f"orig_{seg_img_path.stem}")
        )
        
        with pf.open(seg_img_path) as hdul:
            self.seg_img = hdul[ext].data.byteswap().newbyteorder()
            self.seg_hdr = hdul[ext].header
            self.seg_wcs = WCS(self.seg_hdr)

    def regen_multiband_catalogue(self, **kwargs):

        utils.log_comment(
            utils.LOGFILE,
            "Regenerating multiband catalogue...",
            verbose=True,
            show_date=True,
        )

        [p.unlink() for p in self.out_dir.glob("*phot_apcorr.fits")]

        kwargs["get_all_filters"] = kwargs.get("get_all_filters", True)

        self.catalogue = catalogue_fns.regen_multiband_catalogue(
            self.field_root, 
            seg_image=self.seg_img,
            in_dir=self.out_dir,
            out_dir=self.out_dir,
            **kwargs,
        )

        utils.log_comment(
            utils.LOGFILE,
            "Multiband catalogue complete.",
            verbose=True,
            show_date=True,
        )

    def _modify_FLT_seg(self, flt, seg_file=None):
        print (flt, seg_file)

        flt.process_seg_file(seg_file)


    def load_contamination_maps(self, detection_filter="ir", pad=800, cpu_count=4, **kwargs):

        """
        Be careful with the cpu_count - the memory footprint per process is extremely high 
        (e.g. with an 6P/8E CPU, and 32GB RAM, I typically limit this to <=6 cores).

        """

        grism_files = [str(p) for p in self.out_dir.glob("*GrismFLT.fits")]

        utils.log_comment(
            utils.LOGFILE,
            f"Loading {len(grism_files)} grism files, grizli version={grizli.__version__}",
            verbose=True,
            show_date=True,
        )

        catalog_path = Path(
            kwargs.get(
                "catalog",
                self.out_dir / f"{self.field_root}-{detection_filter}.cat.fits"
            ),
        )
        if not catalog_path.is_file():
            raise FileNotFoundError(
                f"Catalogue file not found at the specified location: {catalog_path}."
            )
        seg_path = Path(
            kwargs.get(
                "seg_file",
                self.out_dir / f"{self.field_root}-{detection_filter}_seg.fits"
            ),
        )
        if not seg_path.is_file():
            raise FileNotFoundError(
                f"Segmentation map not found at the specified location: {seg_path}."
            )

        self.grp = multifit.GroupFLT(
            cpu_count=cpu_count,
            grism_files=grism_files,
            pad=pad,
            seg_file=str(seg_path),
            catalog=str(catalog_path),
        )

    def extract_spectra(self, obj_id_list, z_range=[0.25,0.35], beams_kwargs=None, multibeam_kwargs=None, spectrum_1d=None):

        print("5. Extracting spectra...")
        os.chdir(self.out_dir)
        pline = {"kernel": "square", "pixfrac": 1.0, "pixscale": 0.03, "size": 8, "wcs": None}
        args = auto_script.generate_fit_params(
            pline=pline,
            field_root=self.field_root,
            min_sens=0.0,
            min_mask=0.0,
            # include_photometry=True,  # set both of these to True to include photometry in fitting
            # use_phot_obj=True,
        )  # set both of these to True to include photometry in fitting
        
        if beams_kwargs is None:
            beams_kwargs = {}
        beams_kwargs["size"] = beams_kwargs.get("size", 75)
        beams_kwargs["min_mask"] = beams_kwargs.get("min_mask", 0.)
        beams_kwargs["min_sens"] = beams_kwargs.get("min_sens", 0.)
        beams_kwargs["show_exception"] = beams_kwargs.get("show_exception", True)

        if multibeam_kwargs is None:
            multibeam_kwargs = {}
        multibeam_kwargs["fcontam"] = multibeam_kwargs.get("fcontam", 0.1)
        multibeam_kwargs["min_mask"] = multibeam_kwargs.get("min_mask", 0.)
        multibeam_kwargs["min_sens"] = multibeam_kwargs.get("min_sens", 0.)

        if not hasattr(obj_id_list, '__iter__'):
            obj_id_list = [obj_id_list]

        for obj_id in tqdm(obj_id_list):
            if spectrum_1d is not None:
                beams = FLT_fns.get_beams_with_spectrum(self.grp, obj_id, spectrum_1d=spectrum_1d, **beams_kwargs)
            else:
                beams = self.grp.get_beams(obj_id, **beams_kwargs)

                # with pf.open("/media/sharedData/data/2023_11_07_spectral_orders/Extractions/nis-wfss_02186.full.fits") as hdul:
                #     sp = utils.GTable(hdul['TEMPL'].data)
                #     dt = float
                #     wave = np.cast[dt](sp['wave'])  # .byteswap()
                #     flux = np.cast[dt](sp["full"])  # .byteswap()
                #     flux /= (np.nansum(hdul["DSCI"].data)*1e-19)
                #     print ("DIRECT_SUM:", np.nansum(hdul["DSCI"].data)*1e-19)

            mb = multifit.MultiBeam(beams, group_name=self.field_root, **multibeam_kwargs)
            #     _ = mb.oned_figure()
            #     _ = mb.drizzle_grisms_and_PAs(size=32, scale=0.5, diff=False)
            mb.write_master_fits()
            _ = fitting.run_all_parallel(
                obj_id, zr=z_range, verbose=True, get_output_data=True,
            )
        print("5. Extracting spectra...[COMPLETE]")

    def refine_contam_model_with_fits(self, spectrum="full", max_chinu=5, fit_files=None):
        """
        Refine the full-field grism models with the best fit spectra from 
        individual extractions.
        """

        if fit_files is None:
            fit_files = glob.glob('*full.fits')
            fit_files.sort()
        fit_files = np.atleast_1d(fit_files)
        N = len(fit_files)
        if N == 0:
            return False
        
        msg = 'Refine model ({0}/{1}): {2} / skip (chinu={3:.1f}, dof={4})'
        
        for i, file in enumerate(fit_files):
            try:
                hdu = pf.open(file)
                o_id = hdu[0].header['ID']

                fith = hdu['ZFIT_STACK'].header
                chinu = fith['CHIMIN']/fith['DOF']
                if (chinu > max_chinu) | (fith['DOF'] < 10):
                    print(msg.format(i, N, file, chinu, fith['DOF']))
                    continue

                sp = utils.GTable(hdu['TEMPL'].data)

                dt = float
                wave = np.cast[dt](sp['wave'])  # .byteswap()
                flux = np.cast[dt](sp[spectrum])  # .byteswap()
                self.grp.compute_single_model(int(o_id), mag=19, size=-1, store=False,
                                        spectrum_1d=[wave, flux], is_cgs=True,
                                        get_beams=None, in_place=True)
                print('Refine model ({0}/{1}): {2}'.format(i, N, file))
            except Exception as e:
                print('Refine model ({0}/{1}): {2} / failed {3}'.format(i, N, file, e))

        for f in self.grp.FLTs:
            f.orig_seg = f.seg

    def _create_circular_mask(self, x_c, y_c, radius):

        Y, X = np.ogrid[:self.seg_img.shape[0], :self.seg_img.shape[1]]

        sqrd_dist = (X - x_c)**2 + (Y-y_c)**2

        mask = sqrd_dist <= radius**2
        return mask

    
    def set_obj_circ(self, ra=None, dec=None, x=None, y=None, unit="deg", skycoords=None, radius=1, inner_radius=0):

        if not hasattr(self, "seg_img"):
            raise AttributeError("Segmentation map not set.")

        radius = np.atleast_1d(np.asarray(radius))
        inner_radius = np.atleast_1d(np.asarray(inner_radius))

        if (ra is not None) & (dec is not None):
            try:
                ra = np.atleast_1d(np.asarray(ra))
                dec = np.atleast_1d(np.asarray(dec))
                input_coords = SkyCoord(ra=ra, dec=dec, unit=unit)
            except Exception as e:
                raise Exception(f"Could not parse supplied ra, dec as on-sky coordinates. {e}")
            # try:

            print (self.seg_wcs.footprint_contains(input_coords))
            x_p, y_p = input_coords.to_pixel(self.seg_wcs)

        elif (x is not None) & (y is not None):
            try:
                x = np.asarray(x)
                y = np.asarray(y)
                input_coords = SkyCoord(x=x, y=y)
            except Exception as e:
                raise Exception(f"Could not parse supplied ra, dec as on-sky coordinates. {e}")
        else:
            raise Exception("Coordinate pair not supplied.")

        # print (np.ndim(radius))
        # print (radius.shape)
        # print (repr(radius))
        # print (x_p.shape)
        # print (len(radius))
        if radius.shape[0]<x_p.shape[0]:
            radius = np.array([radius[0]]*x_p.shape[0])
        if inner_radius.shape[0]<x_p.shape[0]:
            inner_radius = np.array([inner_radius[0]]*x_p.shape[0])
        for i, o in zip(inner_radius, radius):
            if i>=o:
                print (i, o)
                raise Exception("Inner radius cannot be greater than outer radius.")

        curr_max = np.nanmax(self.seg_img) + 1
        obj_ids = curr_max + np.arange(radius.shape[0])
        for o_i, x_i, y_i, r_i, i_i in zip(obj_ids, x_p, y_p, radius, inner_radius):
            mask = self._create_circular_mask(x_i,y_i,r_i)
            if i_i != 0:
                mask[self._create_circular_mask(x_i,y_i,i_i)] = 0
            self.seg_img[mask] = o_i
        # fig, ax = plt.subplots()
        # ax.imshow(self.seg_img, origin="lower", cmap="plasma")
        # plt.show()

        return obj_ids

        # try:
        #     print (in1.dtype)
        # except:
        #     print ("no fucking units")

        # try:
        #     input_coords = SkyCoord(ra, dec)
        # except Exception as e:
        #     raise Exception(f"Improper format for coordinates: {e}")

    def segment_split(self, ra=None, dec=None, x=None, y=None, unit="deg", skycoords=None, segments=4, angle=0, radius=0):

        if not hasattr(self, "seg_img"):
            raise AttributeError("Segmentation map not set.")

        radius = np.atleast_1d(np.asarray(radius))

        if (ra is not None) & (dec is not None):
            try:
                ra = np.atleast_1d(np.asarray(ra))
                dec = np.atleast_1d(np.asarray(dec))
                input_coords = SkyCoord(ra=ra, dec=dec, unit=unit)
            except Exception as e:
                raise Exception(f"Could not parse supplied ra, dec as on-sky coordinates. {e}")
            # try:

            print (self.seg_wcs.footprint_contains(input_coords))
            x_p, y_p = input_coords.to_pixel(self.seg_wcs)

        elif (x is not None) & (y is not None):
            try:
                x = np.asarray(x)
                y = np.asarray(y)
                input_coords = SkyCoord(x=x, y=y)
            except Exception as e:
                raise Exception(f"Could not parse supplied ra, dec as on-sky coordinates. {e}")
        else:
            raise Exception("Coordinate pair not supplied.")

        # print (np.ndim(radius))
        # print (radius.shape)
        # print (repr(radius))
        # print (x_p.shape)
        # print (len(radius))
        if radius.shape[0]<x_p.shape[0]:
            radius = np.array([radius[0]]*x_p.shape[0])

        curr_max = np.nanmax(self.seg_img) + 1
        obj_ids = curr_max + np.arange(radius.shape[0], step=int(segments+1))

        used_ids = []
        for o_i, x_i, y_i, r_i in zip(obj_ids, x_p, y_p, radius):

            orig_id = self.seg_img[int(y_i),int(x_i)]
            print (orig_id)

            y, x = np.indices(self.seg_img.shape)

            print (x_i, y_i)
            print (self.seg_img.shape)
            print (np.nanmax(x), np.nanmax(y))

            angle_arr = (np.rad2deg(np.arctan2(x_i-x, y-y_i))-angle)%360.

            # print (360/segments)
            for s in np.arange(segments):
                # print (np.nansum(np.where(
                #     (self.seg_img==orig_id)
                #     &(angle_arr>=s*360/segments)
                #     &(angle_arr<(s+1)*360/segments))))
                # print (o_i+s)
                # print (orig_id)
                self.seg_img[
                    np.where(
                        (self.seg_img==orig_id)
                        &(angle_arr>=s*360/segments)
                        &(angle_arr<(s+1)*360/segments)
                    )
                ] = o_i+s
                used_ids.append(o_i+s)

            # fig, ax = plt.subplots()
            # ax.imshow(self.seg_img[3000:3100,1650:1750], origin="lower")
            # # ax.imshow(angle_arr, origin="lower")
            # plt.show()

            if r_i != 0:
                mask = self._create_circular_mask(x_i,y_i,r_i)
            self.seg_img[mask] = o_i+segments
            used_ids.append(o_i+segments)
        # fig, ax = plt.subplots()
        # ax.imshow(self.seg_img, origin="lower", cmap="plasma")
        # plt.show()

        return used_ids

    def very_hacky_shit(self, old_id, new_ids):

        from copy import deepcopy

        for flt in self.grp.FLTs:
            for n in np.atleast_1d(new_ids).flatten():
                flt.object_dispersers[n] = deepcopy(flt.object_dispersers[old_id])
        
