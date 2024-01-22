
import os
import glob
import numpy as np
import astropy.io.fits as pf
import astropy
import warnings
import numpy.typing as npt
from typing import Optional, Any
from collections.abc import Iterable
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
from datetime import datetime, timezone 
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
from grizli import utils, prep, jwst_utils, multifit, fitting, model
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
model.GrismFLT.transform_JWST_WFSS = FLT_fns.mod_transform_JWST_WFSS
model.GrismFLT.compute_model_orders = FLT_fns.mod_compute_model_orders

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

    def __init__(self, field_root, in_dir, out_dir, seg_path=None):

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

        if seg_path is not None:
            self.load_orig_seg_map(seg_path)



    def load_orig_seg_map(self, seg_path, ext=0):

        seg_path = Path(seg_path)

        self.seg_name = seg_path.name

        # shutil.copy(
        #     src=seg_img_path,
        #     dst=self.out_dir/seg_img_path.with_stem(f"orig_{seg_img_path.stem}")
        # )
        
        with pf.open(seg_path) as hdul:
            self.seg_img = hdul[ext].data.byteswap().newbyteorder()
            self.seg_hdr = hdul[ext].header
            self.seg_wcs = WCS(self.seg_hdr)

    def regen_multiband_catalogue(self, **kwargs):
        """
        Regenerate the grizli multiband catalogue, for the current segmentation map.
        
        Parameters
        ----------
        **kwargs : dict, optional
            Keyword arguments are passed through to 
            `~grizli_functions.catalogue_fns.regen_multiband_catalogue()`.

        """

        utils.log_comment(
            utils.LOGFILE,
            "Regenerating multiband catalogue...",
            verbose=True,
            show_date=True,
        )

        [p.unlink() for p in self.out_dir.glob("*phot_apcorr.fits")]

        kwargs["get_all_filters"] = kwargs.get("get_all_filters", True)

        now = datetime.now(timezone.utc).strftime(r"%Y%m%dT%H%M%SZ")

        try:
            end_string = self.seg_name.split("-")[-1].split(".")[0]
            test = datetime.strptime(end_string, r"%Y%m%dT%H%M%SZ")
            if type(test) == datetime:
                self.seg_name = self.seg_name.replace(end_string,now)
        except:

            self.seg_name = self.seg_name.replace(".fits", f"-{now}.fits")
        
        seg_out_path = self.out_dir / self.seg_name

        self.catalogue = catalogue_fns.regen_multiband_catalogue(
            self.field_root, 
            seg_image=self.seg_img,
            in_dir=self.out_dir,
            out_dir=self.out_dir,
            seg_out_path=seg_out_path,
            **kwargs,
        )

        utils.log_comment(
            utils.LOGFILE,
            f"Multiband catalogue complete. Segmentation map saved to {self.seg_name}",
            verbose=True,
            show_date=True,
        )

    def load_contamination_maps(self, grism_files=None, detection_filter="ir", pad=800, cpu_count=4, **kwargs):

        '''
        Be careful with the cpu_count - the memory footprint per process is extremely high 
        (e.g. with a 6P/8E CPU, and 32GB RAM, I typically limit this to <=6 cores).

        '''

        if grism_files is None:
            grism_files = [str(p) for p in self.out_dir.glob("*GrismFLT.fits")]
        else:
            grism_files = np.atleast_1d(grism_files)
        if len(grism_files)==0:
            raise Exception("No grism files found.")

        utils.log_comment(
            utils.LOGFILE,
            f"Loading {len(grism_files)} grism files, grizli version={grizli.__version__}",
            verbose=True,
            show_date=True,
        )

        catalog_path = Path(
            kwargs.get(
                "catalog_path",
                self.out_dir / f"{self.field_root}-{detection_filter}.cat.fits"
            ),
        )
        if not catalog_path.is_file():
            raise FileNotFoundError(
                f"Catalogue file not found at the specified location: {catalog_path}."
            )
        seg_path = kwargs.get("seg_path")
        if seg_path is None and hasattr(self, "seg_name"):
            seg_path = self.out_dir / self.seg_name
        if seg_path is not None:
            if not Path(seg_path).is_file():
                raise FileNotFoundError(
                    f"Segmentation map not found at the specified location: {seg_path}."
                )
            seg_path = str(seg_path)

        self.grp = multifit.GroupFLT(
            cpu_count=cpu_count,
            grism_files=grism_files,
            pad=pad,
            seg_file=seg_path,
            catalog=str(catalog_path),
        )

    def extract_spectra(self, obj_id_list: npt.ArrayLike, z_range: npt.ArrayLike = [0., 0.5], beams_kwargs: dict[str, Any] | None = None, multibeam_kwargs: dict[str, Any] | None = None, spectrum_1d: npt.ArrayLike | None = None, is_cgs: bool = True):
        """
        Perform a full extraction of the specified objects.

        Parameters
        ----------
        obj_id_list : array_like
            The object ids in the segmentation map which will be extracted.
        z_range : array_like, optional
            The redshift range to consider for the extraction, by default 0 < z < 0.5
        beams_kwargs : dict, optional
            Keyword arguments to pass to grizli.multifit.GroupFLT.get_beams()
        multibeam_kwargs : dict, optional
            Keyword arguments to pass to grizli.multifit.MultiBeam()
        spectrum_1d : [wavelengths, flux], optional
            The flux spectrum and corresponding wavelengths of the object in the model.
            By default, this is calculated automatically from the stored object_dispersers.
        is_cgs : bool, optional
            The flux units of `spectrum_1d[1]` are cgs f_lambda flux densities,
            rather than normalised in the detection band, by default True.
        """

        if not hasattr(self, "grp"):
            raise Exception("GrismFLT files not loaded. Run `load_contamination_maps()' first.")
        if Path(self.grp.FLTs[0].seg_file).name != self.seg_name:
            raise Exception(
                f"The current segmentation map ({self.seg_name}) does not match the"
                f"name stored in the GrismFLT files ({Path(self.grp.FLTs[0].seg_file).name})."
                "Run `load_contamination_maps()' before extracting any spectra, or load "
                "the correct segmentation map."
            )
            
        utils.log_comment(
            utils.LOGFILE,
            f"Generating fit parameters.",
            verbose=True,
            show_date=True,
        )
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
        beams_kwargs["size"] = beams_kwargs.get("size", -1)
        beams_kwargs["min_mask"] = beams_kwargs.get("min_mask", 0.)
        beams_kwargs["min_sens"] = beams_kwargs.get("min_sens", 0.)
        beams_kwargs["show_exception"] = beams_kwargs.get("show_exception", True)

        if multibeam_kwargs is None:
            multibeam_kwargs = {}
        multibeam_kwargs["fcontam"] = multibeam_kwargs.get("fcontam", 0.1)
        multibeam_kwargs["min_mask"] = multibeam_kwargs.get("min_mask", 0.)
        multibeam_kwargs["min_sens"] = multibeam_kwargs.get("min_sens", 0.)

        # if not hasattr(obj_id_list, '__iter__'):
        #     obj_id_list = [obj_id_list]
        obj_id_arr = np.atleast_1d(obj_id_list).flatten()

        for obj_id in tqdm(obj_id_arr):
            beams = FLT_fns.get_beams_with_spectrum(self.grp, obj_id, spectrum_1d=spectrum_1d, is_cgs=is_cgs, **beams_kwargs)

            mb = multifit.MultiBeam(beams, group_name=self.field_root, **multibeam_kwargs)
            mb.write_master_fits()
            _ = fitting.run_all_parallel(
                obj_id, zr=z_range, verbose=True, get_output_data=True,
            )        
        utils.log_comment(
            utils.LOGFILE,
            f"Finished extracting spectra.",
            verbose=True,
            show_date=True,
        )

    def refine_contam_model_with_fits(self, spectrum: str = "full", max_chinu: int | float = 5, 
        fit_files: list[str] | list[os.PathLike] | None = None) -> bool:
        """
        Refine the full-field grism models with the best fit spectra from 
        individual extractions. [Modified version of a grizli function]

        Parameters
        ----------
        spectrum : str, optional
            The component of the best-fit spectrum to use, either `full' or `continuum'.
        max_chinu : int | float, optional
            The maximum reduced chi-squared value of the fit to accept, in order to 
            refine the contamination model with the resulting spectrum, by default 5.
        fit_files : list[str] | list[os.PathLike] | None, optional
            An explicit list of the best-fit files to use. By default, all 
            `*full.fits' files in the current directory will be used.

        Returns
        -------
        status : bool
            Returns False if the contamination maps are not modified.
        """

        if fit_files is None:
            fit_files = glob.glob('*full.fits')
            fit_files.sort()
        fit_files_arr = np.atleast_1d(np.asarray(fit_files))
        N = fit_files_arr.shape[0]
        if N == 0:
            return False
        
        msg = 'Refine model ({0}/{1}): {2} / skip (chinu={3:.1f}, dof={4})'
        
        for i, file in enumerate(fit_files_arr):
            try:
                hdu = pf.open(file)
                o_id = hdu[0].header['ID']

                fith = hdu['ZFIT_STACK'].header
                chinu = fith['CHIMIN']/fith['DOF']
                if (chinu > max_chinu) | (fith['DOF'] < 10):
                    print(msg.format(i, N, file, chinu, fith['DOF']))
                    continue

                sp = utils.GTable(hdu['TEMPL'].data)

                wave = np.cast[float](sp['wave'])  # .byteswap()
                flux = np.cast[float](sp[spectrum])  # .byteswap()
                for flt in self.grp.FLTs:
                    if int(o_id) not in flt.object_dispersers:
                        old_obj_ids = np.unique(flt.orig_seg[flt.seg==o_id])
                        old_obj_ids = old_obj_ids.ravel()[np.flatnonzero(old_obj_ids)].astype(int)
                self.grp.compute_single_model(int(o_id), mag=19, size=-1, store=False,
                                        spectrum_1d=[wave, flux], is_cgs=True,
                                        get_beams=None, in_place=True)
                print('Refine model ({0}/{1}): {2}'.format(i, N, file))
            except Exception as e:
                print('Refine model ({0}/{1}): {2} / failed {3}'.format(i, N, file, e))

        for f in self.grp.FLTs:
            f.orig_seg = f.seg
            f.orig_seg_file = f.seg_file

        return True


    def _create_circular_mask(self, x_c: float, y_c: float, radius: float) -> npt.NDArray[np.bool_]:
        """
        Create a boolean mask of all elements in the segmentation map within 
        a specified distance of a point.

        Parameters
        ----------
        x_c : float
            The x-coordinate of the reference point.
        y_c : float
            The y-coordinate of the reference point.
        radius : float
            The maximum radius allowed.

        Returns
        -------
        mask : ndarray, bool
            The mask, where elements are True if the distance to 
            (x_c, y_c) is less than or equal to radius.
        """

        Y, X = np.ogrid[:self.seg_img.shape[0], :self.seg_img.shape[1]]

        sqrd_dist = (X - x_c)**2 + (Y-y_c)**2

        mask = sqrd_dist <= radius**2
        return mask

    def add_circ_obj(self, radius: astropy.units.Quantity | float | Iterable[astropy.units.Quantity | float] = 3*u.arcsec, 
        inner_radius: astropy.units.Quantity | float | Iterable[astropy.units.Quantity | float] = 0, 
        centre: astropy.coordinates.SkyCoord = None, **skycoord_kwargs):

        """
        Add one or more circular objects to the segmentation map.

        Parameters
        ----------
        radius : astropy.units.Quantity | float | Iterable[astropy.units.Quantity  |  float], optional
            The outer radius of the aperture, by default 3 arcseconds.
        inner_radius : astropy.units.Quantity | float | Iterable[astropy.units.Quantity  |  float], optional
            The inner radius, if specified, creates an annulus instead of an aperture.
        centre : astropy.coordinates.SkyCoord, optional
            The centre of the aperture.
        **skycoord_kwargs : dict, optional
            Any inputs accepted by astropy.coordinates.SkyCoord, if skycoords is None.

        Returns
        -------
        new_obj_id : int
            The id corresponding to the new object in the segmentation map.
        """

        if not hasattr(self, "seg_img"):
            raise AttributeError("Segmentation map not set.")

        if centre is None:
            try:
                centre = SkyCoord(**skycoord_kwargs)
            except Exception as e:
                raise Exception(f"Could not parse supplied arguments as on-sky coordinates: {e}")

        centre = np.atleast_1d(centre).flatten()
        outer_radius = np.atleast_1d(radius).flatten()
        inner_radius = np.atleast_1d(inner_radius).flatten()

        for r in [inner_radius, outer_radius]:
            if (r.shape[0]>1) and (r.shape[0] != centre.shape[0]):
                raise ValueError(
                    f"Size of inputs do not match. {centre.shape[0]}"
                    f" coordinates and {r.shape[0]} radii have been supplied."
                )

        if (inner_radius.shape[0]==1) and (centre.shape[0]>1):
            try:
                inner_radius = np.tile(inner_radius.value,centre.shape[0])*inner_radius.unit
            except:
                inner_radius = np.tile(inner_radius, centre.shape[0])
        if (outer_radius.shape[0]==1) and (centre.shape[0]>1):
            try:
                outer_radius = np.tile(outer_radius.value,centre.shape[0])*outer_radius.unit
            except:
                outer_radius = np.tile(outer_radius, centre.shape[0])

        for i, o in zip(inner_radius, outer_radius):
            if i>=o:
                raise ValueError("Inner radius cannot be greater than outer radius.")

        contained = np.asarray([self.seg_wcs.footprint_contains(c) for c in centre])
        if not all(contained):
            warnings.warn(
                f"The following coordinates are outside the segmentation map footprint, "
                f"and will be skipped: {centre[~contained].to_string()}",
            )
        xs, ys = centre[contained].to_pixel(self.seg_wcs)
        
        inner_radius, outer_radius = inner_radius[contained], outer_radius[contained]

        radii = np.zeros((xs.shape[0], 2))
        pix_scale = np.nanmean(
            [d.value for d in self.seg_wcs.proj_plane_pixel_scales()]
        ) * self.seg_wcs.proj_plane_pixel_scales()[0].unit /u.pix
        for n in range(xs.shape[0]):
            for i, r in enumerate([inner_radius[n], outer_radius[n]]):
                if isinstance(r, u.Quantity):
                    if u.get_physical_type(r)=="angle":
                        radii[n, i] = (r / pix_scale).to(u.pix).value
                    else:
                        radii[n, i] = r.value
                else:
                    radii[n, i] = r

        curr_max = np.nanmax(self.seg_img) + 1
        new_obj_ids = curr_max + np.arange(xs.shape[0])
        for new_id, x_c, y_c, rads in zip(new_obj_ids, xs, ys, radii):
            mask = self._create_circular_mask(x_c, y_c, rads[1])
            if rads[0] != 0:
                mask[self._create_circular_mask(x_c, y_c, rads[0])] = 0
            self.seg_img[mask] = new_id

        return new_obj_ids
            

    
    # def add_circ_obj(self, ra=None, dec=None, x=None, y=None, unit="deg", skycoords=None, radius=1, inner_radius=0):

    #     if not hasattr(self, "seg_img"):
    #         raise AttributeError("Segmentation map not set.")

    #     radius = np.atleast_1d(np.asarray(radius))
    #     inner_radius = np.atleast_1d(np.asarray(inner_radius))

    #     if (ra is not None) & (dec is not None):
    #         try:
    #             ra = np.atleast_1d(np.asarray(ra))
    #             dec = np.atleast_1d(np.asarray(dec))
    #             input_coords = SkyCoord(ra=ra, dec=dec, unit=unit)
    #         except Exception as e:
    #             raise Exception(f"Could not parse supplied ra, dec as on-sky coordinates. {e}")
    #         if not self.seg_wcs.footprint_contains(input_coords):
    #             raise Exception(f"Supplied coordinates are outside the segmentation map footprint.")
    #         x_p, y_p = input_coords.to_pixel(self.seg_wcs)

    #     elif (x is not None) & (y is not None):
    #         try:
    #             x = np.asarray(x)
    #             y = np.asarray(y)
    #             input_coords = SkyCoord(x=x, y=y)
    #         except Exception as e:
    #             raise Exception(f"Could not parse supplied ra, dec as on-sky coordinates. {e}")
    #     else:
    #         raise Exception("Coordinate pair not supplied.")

    #     if radius.shape[0]<x_p.shape[0]:
    #         radius = np.array([radius[0]]*x_p.shape[0])
    #     if inner_radius.shape[0]<x_p.shape[0]:
    #         inner_radius = np.array([inner_radius[0]]*x_p.shape[0])
    #     for i, o in zip(inner_radius, radius):
    #         if i>=o:
    #             print (i, o)
    #             raise Exception("Inner radius cannot be greater than outer radius.")

    #     curr_max = np.nanmax(self.seg_img) + 1
    #     obj_ids = curr_max + np.arange(radius.shape[0])
    #     for o_i, x_i, y_i, r_i, i_i in zip(obj_ids, x_p, y_p, radius, inner_radius):
    #         mask = self._create_circular_mask(x_i,y_i,r_i)
    #         if i_i != 0:
    #             mask[self._create_circular_mask(x_i,y_i,i_i)] = 0
    #         self.seg_img[mask] = o_i

    #     return obj_ids

    def add_segment_obj(self, ra=None, dec=None, x=None, y=None, unit="deg", skycoords=None, segments=4, angle=0, radius=0):

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

        if radius.shape[0]<x_p.shape[0]:
            radius = np.array([radius[0]]*x_p.shape[0])

        curr_max = np.nanmax(self.seg_img) + 1
        obj_ids = curr_max + np.arange(radius.shape[0], step=int(segments+1))

        used_ids = []
        for o_i, x_i, y_i, r_i in zip(obj_ids, x_p, y_p, radius):

            orig_id = self.seg_img[int(y_i),int(x_i)]

            y, x = np.indices(self.seg_img.shape)

            angle_arr = (np.rad2deg(np.arctan2(x_i-x, y-y_i))-angle)%360.

            for s in np.arange(segments):
                self.seg_img[
                    np.where(
                        (self.seg_img==orig_id)
                        &(angle_arr>=s*360/segments)
                        &(angle_arr<(s+1)*360/segments)
                    )
                ] = o_i+s
                used_ids.append(o_i+s)
            if r_i != 0:
                mask = self._create_circular_mask(x_i,y_i,r_i)
            self.seg_img[mask] = o_i+segments
            used_ids.append(o_i+segments)

        return used_ids
        
    def add_reg_obj(self, reg_path: str | os.PathLike, format: str | None = None, reg_wcs: astropy.wcs.WCS | None = None) -> int:

        """
        Read a regions file, and add the relevant region to the 
        segmentation map as a new object.

        Parameters
        ----------
        reg_path : str | os.PathLike
            The path pointing to the region
        format : str | None, optional
            The file format specifier. If None, the format is 
            automatically inferred from the file extension.
        reg_wcs : astropy.wcs.WCS | None, optional
            The WCS to use to convert pixels to world coordinates. 
            By default, the segmentation map WCS will be used.

        Returns
        -------
        new_obj_id : int
            The id corresponding to the new object in the segmentation map.

        Raises
        ------
        AttributeError
            If the segmentation map has not been loaded.
        ImportError
            If the regions package is not available.
        Exception
            If the supplied regions file cannot be read.
        """

        if not hasattr(self, "seg_img"):
            raise AttributeError("Segmentation map not set.")

        try:
            from regions import Regions, SkyRegion
        except:
            raise ImportError(
                "Astropy Regions is required to import region files."
                "Please read https://astropy-regions.readthedocs.io for more information."
            )
        
        try:
            region = Regions.read(reg_path, format=format).regions[0]
        except Exception as e:
            raise Exception(
                f"Could not parse regions file, with the following error: {e}"
            )
        
        if issubclass(type(region), SkyRegion):
            pixel_region = region.to_pixel(self.seg_wcs if reg_wcs is None else reg_wcs)
        else:
            pixel_region = region
        
        mask = pixel_region.to_mask(mode="subpixels")

        matched_mask = mask.to_image(self.seg_img.shape) > 0.5

        new_obj_id = np.nanmax(self.seg_img)+1
        self.seg_img[matched_mask] = new_obj_id

        return new_obj_id