
import os
import glob
import numpy as np
import astropy.io.fits as pf
import astropy.wcs as pywcs
import matplotlib.pyplot as plt
# from grizli import utils, prep, jwst_utils, multifit
# from grizli.pipeline import auto_script
# from grizli.pipeline.auto_script import get_yml_parameters
# import grizli
from importlib import reload
# from grizli import model, multifit, grismconf
# from grizli import fitting
import astropy.units as u
from tqdm import tqdm
from pathlib import Path
import pickle
import sep
import math
# import eazy

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

class GrizliExtractor:

    def __init__(self, field_name, in_dir, out_dir):

        self.field_name = field_name
        self.in_dir = Path(in_dir)
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(exist_ok=True, parents=True)

        grizli_dir = self.out_dir / "grizli"
        grizli_dir.mkdir(exist_ok=True)
        (grizli_dir / "CONF").mkdir(exist_ok=True)
        (grizli_dir / "templates").mkdir(exist_ok=True)
        (grizli_dir / "iref").mkdir(exist_ok=True)
        (grizli_dir / "iref").mkdir(exist_ok=True)
        os.environ["GRIZLI"] = str(grizli_dir)
        os.environ["iref"] = str(grizli_dir / "iref")
        os.environ["jref"] = str(grizli_dir / "jref")

        # print (os.environ["GRIZLI"])

        # print (grizli.GRIZLI_PATH)

        # grizli.GRIZLI_PATH = os.environ["GRIZLI"]

        # reload (grizli)
        # print (grizli.GRIZLI_PATH)

        import grizli
        from grizli import utils, prep, jwst_utils, multifit, fitting
        from grizli.pipeline import auto_script
        import eazy
        global grizli, utils, prep, jwst_utils, multifit, eazy, auto_script, fitting


        os.chdir(os.path.join(grizli.GRIZLI_PATH, 'CONF'))
        if not os.path.exists('GR150C.F115W.221215.conf'):
            os.system('wget "https://zenodo.org/record/7628094/files/niriss_config_221215.tar.gz?download=1" -O niriss_config_221215.tar.gz')
            os.system('tar xzvf niriss_config_221215.tar.gz')

        if not os.path.exists('niriss_sens_221215.tar.gz'):
            os.system('wget "https://zenodo.org/record/7628094/files/niriss_sens_221215.tar.gz" -O niriss_sens_221215.tar.gz')
            os.system('tar xzvf niriss_sens_221215.tar.gz')

        os.chdir(self.in_dir)

        # grizli.utils.fetch_default_calibs()
        # print ("done1")
        # grizli.utils.fetch_config_files(get_jwst=True)

        # eazy.fetch_eazy_photoz()

        # try:
        #     grizli.utils.symlink_templates(force=True)
        # except:
        #     pass

    def load_contamination_maps(self):


        print (grizli.__version__)

        grism_files = [str(p) for p in self.in_dir.glob("*GrismFLT.fits")]

        for g in grism_files:

            if os.path.exists(g) & ('GrismFLT' in g):
                print(f'Checking conf file location for {g}.')

                # fp = open(g.replace('GrismFLT.fits', 'GrismFLT.pkl'), 'rb')
                # print (g.replace('GrismFLT.fits', 'GrismFLT.pkl'))
                # flt = pickle.load(fp)
                # fp.close()

                with open(g.replace('GrismFLT.fits', 'GrismFLT.pkl'), 'rb') as fp:
                    flt = pickle.load(fp)
                    # keys = [*flt.grism.header.keys()]
                    # keys.sort()
                    # print (keys)
                # flt = pickle.load(fp)
                # fp.close()
                # print (dir(flt))
                # print (flt.conf_file)

                flt.conf_file = f"{os.environ['GRIZLI']}/CONF/{flt.conf_file.split('/')[-1]}"
                # print (flt.conf_file)

                with open(g.replace('GrismFLT.fits', 'GrismFLT.pkl'), 'wb') as fp:
                    pickle.dump(flt, fp)

                # status = flt.load_from_fits(save_file)

            # with pf.open(g) as g_hdul:
            #     g_hdul.info()


        self.grp = multifit.GroupFLT(
            grism_files=grism_files,
            catalog=str(self.in_dir / f"{self.field_name}-ir.cat.fits"),
            cpu_count=-1,
            sci_extn=1,
            pad=800,
        )
        # print (dir(self.grp))

    def extract_spectra(self, obj_id_list):

        print("5. Extracting spectra...")
        os.chdir(self.out_dir)
        pline = {"kernel": "square", "pixfrac": 1.0, "pixscale": 0.03, "size": 8, "wcs": None}
        args = auto_script.generate_fit_params(
            pline=pline,
            field_root=self.field_name,
            min_sens=0.0,
            min_mask=0.0,
            include_photometry=True,  # set both of these to True to include photometry in fitting
            use_phot_obj=True
        )  # set both of these to True to include photometry in fitting
        for obj_id in tqdm(obj_id_list):
            beams = self.grp.get_beams(obj_id, size=32, min_mask=0, min_sens=0)
            mb = multifit.MultiBeam(beams, fcontam=0.1, min_sens=0, min_mask=0, group_name=self.field_name)
            #     _ = mb.oned_figure()
            #     _ = mb.drizzle_grisms_and_PAs(size=32, scale=0.5, diff=False)
            mb.write_master_fits()
            _ = fitting.run_all_parallel(
                obj_id, zr=[0.29,0.30], verbose=True, get_output_data=True, skip_complete=False
            )
        print("5. Extracting spectra...[COMPLETE]")


    def extract_sep(self, obj_id, threshold=1.8):

        x = [0,6573]
        y = [0,6429]

        # x = [2750,3250]
        # y = [1450,2000]

        # x = [2916,2960]
        # y = [1881,1926]

        # x = [2818,2909]
        # y = [2208,2305]

        with pf.open("/media/sharedData/data/2023_11_07_spectral_orders/Prep/nis-wfss-ir_seg.fits") as hdul:

            seg_img = hdul[0].data[x[0]:x[1],y[0]:y[1]].byteswap().newbyteorder()

        with pf.open("/media/sharedData/data/2023_11_07_spectral_orders/Prep/nis-wfss-ir_bkg.fits") as hdul:

            bkg_img = hdul[0].data[x[0]:x[1],y[0]:y[1]].byteswap().newbyteorder()
            
        with pf.open("/media/sharedData/data/2023_11_07_spectral_orders/Prep/nis-wfss-ir_drz_sci.fits") as hdul:

            det_img = hdul[0].data[x[0]:x[1],y[0]:y[1]].byteswap().newbyteorder()
            
        with pf.open("/media/sharedData/data/2023_11_07_spectral_orders/Prep/nis-wfss-ir_drz_wht.fits") as hdul:

            wht_img = hdul[0].data[x[0]:x[1],y[0]:y[1]].byteswap().newbyteorder()

        import matplotlib.pyplot as plt

        # fig, ax = plt.subplots()
        # ax.imshow(seg_img)
        # plt.show()

        # seg_img[seg_img!=2125] = 0

        # area = np.where(seg_img==obj_id)

        # details = {}

        # details["NPIX"] = len(area[0])
        # details["XMIN"] = np.nanmin(area[1])
        # details["XMAX"] = np.nanmax(area[1])
        # details["YMIN"] = np.nanmin(area[0])
        # details["YMAX"] = np.nanmax(area[0])

        # print (np.nanmean(area[0]))
        # print (np.nanmin(det_img[area]))

        # print (np.nanmin((wht_img*(det_img-bkg_img)**2)[area]))

        # non_sky_pix = (det_img-bkg_img)

        # wht_img = np.sqrt(wht_img)

        # print ("Weight", np.nanmean(wht_img[area]))
        # print ("Detection", np.nanmean(det_img[area]))

        # print ("bkg", np.nanmean(np))


        # print (np.nanmin((det_img*wht_img)[area]))

        err = 1/np.sqrt(wht_img)
        # True mask pixels are masked with sep
        mask = (~np.isfinite(err)) | (err == 0) | (~np.isfinite(det_img))
        err[mask] = 0

        # fig, ax = plt.subplots()
        # ax.imshow(mask)
        # plt.show()

        bkg_input={'bw': 128, 'bh': 128, 'fw': 3, 'fh': 3}

        bkg = sep.Background(det_img, mask=mask, **bkg_input)
        bkg_data = bkg.back()
        data_bkg = det_img - bkg_data
        det_img_sub = data_bkg
        # err = bkg.rms()
        ratio = bkg.rms()/err
        err_scale = np.median(ratio[(~mask) & np.isfinite(ratio)])
        # print (err_scale)

        err *= err_scale

        # print ("ERROR_STATS: ", err_scale, np.nanmin(err), np.nanmedian(err), np.nanstd(err))
        # print ("ERROR_STATS: ", np.nanmean(data_bkg), err_scale, np.nanmin(data_bkg), np.nanmedian(data_bkg), np.nanstd(data_bkg))
        # exit

        # sep.set_extract_pixstack(int(3e7))
        # sep.set_sub_object_limit(4096)
        # objects, seg = sep.extract(
        #     data_bkg, 
        #     threshold, 
        #     err=err,
        #     mask=mask, 
        #     segmentation_map=True,
        #     **prep.SEP_DETECT_PARAMS
        # )

        print (seg_img.dtype)

        # Look at updating info[co] rather than setting luflag

        # exit()

        from grizliForceExtract.utils_c import sep_fns

        # data_bkg = np.zeros()

        # prep.SEP_DETECT_PARAMS['minarea'] = 1
        # prep.SEP_DETECT_PARAMS['clean'] = False
        # prep.SEP_DETECT_PARAMS['filter_kernel'] = None

        # print (prep.SEP_DETECT_PARAMS)



        # data_bkg = np.array(
        #     # [
        #     #     [0,0,0,1,0],
        #     #     [0,0,0,1,0],
        #     #     [0,1,0,1,0],
        #     #     [1,1,0,0,0],
        #     #     [0,1,0,0,0],
        #     # ],
        #     [
        #         [0,1,1,1,0],
        #         [0,0,1,1,0],
        #         [0,0,0,1,0],
        #         [0,1,0,1,0],
        #         [0,1,0,0,0],
        #     ],
        #     dtype=float,
        # )#.byteswap().newbyteorder()

        # seg_img = np.array(
        #     # [
        #     #     [0,0,0,1,0],
        #     #     [0,0,0,1,0],
        #     #     [0,1,0,1,0],
        #     #     [1,1,0,0,0],
        #     #     [0,1,0,0,0],
        #     # ],
        #     [
        #         [0,1,1,1,0],
        #         [0,0,1,1,0],
        #         [0,0,0,1,0],
        #         [0,2,0,1,0],
        #         [0,2,0,0,0],
        #     ],
        #     dtype=float,
        # )#.byteswap().newbyteorder()

        # err = np.full_like(data_bkg, 0.01, dtype=float)#.byteswap().newbyteorder()
        # mask = np.zeros_like(data_bkg, dtype=float)#.byteswap().newbyteorder()

        sep_fns.set_extract_pixstack(int(3e7))
        sep_fns.set_sub_object_limit(4096)
        objects, seg = sep_fns.extract(
            data_bkg, 
            threshold, 
            err=err,
            mask=mask, 
            segmentation_map=seg_img,
            # segmentation_map=True,
            **prep.SEP_DETECT_PARAMS
        )

        # print (len(objects))
        print ("Objects", objects)
        print (seg)

        output_markers = np.array(
            [
                [" "," "," ","S","F"],
                [" "," "," ","S","F"],
                [" ","S","F","S","F"],
                ["S"," ","F"," "," "],
                [" ","S","F"," "," "],
            ],
            # dtype=float,
        )#.byteswap().newbyteorder()

        # # from astropy.table import Table
        # # test_tab = Table.read("info.txt", format="ascii.no_header")
        # # # print (test_tab)
        # # test_array = np.zeros_like(seg_img)
        # # # print (np.max(test_tab["col1"]))
        # # # print (np.max(test_tab["col2"]))
        # # # print (test_array.shape)
        # # for row in test_tab:
        # #     try:
        # #         test_array[row[1], row[0]] = row[2]
        # #     except IndexError:
        # #         pass

        fig, ax = plt.subplots()
        seg = seg.astype(float)
        seg[seg==0] = np.nan
        # im = ax.imshow(test_array, origin="lower")
        # plt.colorbar(im)
        # plt.show()
        ax.set_facecolor("0.7")
        ax.imshow(seg, cmap="plasma")
        plt.show()

        # print (err[3000:3005, 3000:30005])
        # from astropy.convolution import convolve_fft as convolve
#         from astropy.convolution import convolve
#         print (prep.GAUSS_3_7x7)

        # self.convolve_matched(det_img_sub, err, prep.GAUSS_3_7x7)

#         det_img_orig = np.copy(det_img_sub)
#         det_img_sub = convolve(det_img_sub, prep.GAUSS_3_7x7)
#         err_conv = convolve(err, prep.GAUSS_3_7x7)

#         print (np.nanmedian(err[area])*1.8)
#         print (np.nanmean(err[area])*1.8)
#         thresholds = err*2
#         print (sum(mask[area]))
#         # thresh_mask = (det_img_sub > thresholds) & (seg_img==obj_id)

#         # seg_img[~thresh_mask] = 0
#         # new_area = np.where(seg_img==obj_id)

#         new_area = area
# # ERROR_STATS:  0.0 0.00046226996 0.0034473338
#         rv = np.nansum(det_img_sub[new_area])

#         # det_img_sub = det_img
#         # det_img_sub += thresholds
#         print (np.nanmax(det_img_sub[new_area]))

        
#         mx = my = tv = 0.0
#         mx2 = my2 = mxy = 0.0

#         mean_thresh = np.nanmean(err[area])*1.8
#         print ("THRESHOLDS?")
#         print (np.nanmin(err[area]*1.8))
#         print (np.nanmedian(err[area]*1.8))
#         print (np.nanmean(err_conv[area]*1.8))
#         print (np.nanmin(err_conv[area]*1.8))
#         print (np.nanmedian(err_conv[area]*1.8))
#         print (np.nanmean(err_conv[area]*1.8))
#         print ("THRESHOLDS over")

#         for y, x in zip(new_area[0], new_area[1]):
#             # if det_img_orig[y,x]>1.8*err[y,x]:
#             # if det_img_sub[y,x]>1.8*err_conv[y,x]:
#             if det_img_sub[y,x]> 0.004350837785750628:

#             # if det_img_orig[y,x]>thresholds[y,x]:
#             # if det_img_sub[y,x]>mean_thresh:
#             # if det_img_sub[y,x]>np.nanmean(err[area])*1.8:
#                 tv+=1
#                 # print (det_img_sub[y,x], thresholds[y,x])
#                 # print (a)
#             cval = det_img_sub[y,x]
#             x -= details["XMIN"]
#             y -= details["YMIN"]

#             # mx += x*det_img_sub[y,x]
#             # my += y*det_img_sub[y,x]
#             # mx2 += x*x*det_img_sub[y,x]
#             # my2 += y*y*det_img_sub[y,x]
#             # mxy += x*y*det_img_sub[y,x]

#             mx += cval * x
#             my += cval * y
#             mx2 += cval * x*x
#             my2 += cval * y*y
#             mxy += cval * x*y

#         xm = mx/rv
#         ym = my/rv

#         xm2 = mx2 / rv - xm * xm
#         ym2 = my2 / rv - ym * ym
#         xym = mxy / rv - xm * ym

#         print (xm+details["XMIN"], ym+details["YMIN"], xm2, ym2, xym)
#         print (xm2*ym2-xym*xym)
#         print (details, tv)
#         print (rv)
#         # print (rv/np.nansum(err[new_area]))
#         # print (np.nanmean(new_area[1]*det_img_sub))
#         # print (np.nanmean(area[1]))
#         # print (1.8*err[area])


#         # print (bkg.rms()[3000:3005, 3000:30005])

#         # # det_test = 
#         # bkg = sep.Background(det_img)
#         # # print (bkg.back()[3000:3005, 3000:30005])
#         # # print ((bkg.back()**2)[3000:3005, 3000:30005])
#         # # print (np.nanmax(bkg))
#         # # bkg_rms = bkg.rms()
#         # print (bkg.rms()[3000:3005, 3000:30005])
#         # print (np.nanmax(bkg_rms))


#         # details["thresh"] = 

    @staticmethod
    def convolve_matched(img, noise, conv):

        print ("\nBEGIN MATCHED FILTERING\n")
        # sum(conv_i * f_i / n_i^2) / sqrt(sum(conv_i^2 / n_i^2))

        convw, convh = conv.shape
        convw2 = math.trunc(convw/2)

        out = np.zeros_like(img)

        print (convw)

        # for y in range(img.shape[1]):

        for y in [7]:

            y0 = y-math.trunc(convh/2)

        
            # /* Cut off top of kernel if it extends beyond image */
            if y0 + convh > img.shape[1]:
                convh = img.shape[1] - y0

            # /* cut off bottom of kernel if it extends beyond image */
            if (y0 < 0):
                convh = convh + y0
                conv += convw*(-y0)
                y0 = 0

            convn = convw * convh

            print (convn)

            for i in range (convn):
                cx = i % convw #  /* x index in conv kernel */
                cy = math.trunc(i / convw) #  /* y index in conv kernel */
                print (y0, y, cy)
                print (cx, cy)
                imline = img[:,y0 - y + cy]
                nline = noise[:,y0 - y + cy]

                # /* get start and end positions in the source and target line */
                dcx = cx - convw2 # /* offset of conv pixel from conv center;
                                  #     determines offset between in and out line */
                if (dcx >= 0):
                
                    src_im = imline + dcx
                    src_n = nline + dcx
                    dst_num = out
                    dst_denom = work
                    dst_num_end = outend - dcx
                
                else:
                    src_im = imline
                    src_n = nline
                    dst_num = out - dcx
                    dst_denom = work - dcx
                    dst_num_end = outend

                print (src_im)

                # # /* actually calculate values */
                # while (dst_num < dst_num_end):
                #     imval = src_im
                #     varval = src_n**2
                #     if (varval != 0.0)
                #         {
                #         *dst_num   += conv[i] * imval / varval;
                #         *dst_denom += conv[i] * conv[i] / varval;
                #         }
                #     src_im+;
                #     src_n++;
                #     dst_num++;
                #     dst_denom++;
                #     }



# print("Grizli version: ", grizli.__version__)

# print("1. Importing Python packages...[COMPLETE]")
# ############################################
# ### Define necessary variables and paths ###
# ############################################
# print("2. Defining variables and paths...")

# HOME_PATH = "/media/sharedData/data/GLASS_owncloud/NIRISS/ABELL2744/v3"  # change this
# root = "nis-wfss"

# os.chdir(HOME_PATH + "/Prep")

# print("2. Defining variables and paths...[COMPLETE]")

if __name__=="__main__":

    ge = GrizliExtractor(
        "nis-wfss", 
        # "/media/sharedData/data/2023_11_07_spectral_orders/Prep",
        # "/media/sharedData/data/2023_11_07_spectral_orders/ForcedExtractions",
        # "/media/sharedData/data/GLASS_owncloud/NIRISS/ABELL2744/v3/Prep",
        # "/media/sharedData/data/GLASS_owncloud/NIRISS/ABELL2744/v3/ForcedExtractions",
        "/media/sharedData/data/2023_11_07_spectral_orders/Prep_testing_v3",
        "/media/sharedData/data/2023_11_07_spectral_orders/ForcedExtractions",
    )

    ge.extract_sep(4)

    # ge.load_contamination_maps()
    # ge.extract_spectra([2125,2126,2127,2128,2129,2130,2132,2133,2134])

    # import eazy
    # import grizli
    # import os

    # print (os.environ["GRIZLI"])
    # os.environ["GRIZLI"] = 