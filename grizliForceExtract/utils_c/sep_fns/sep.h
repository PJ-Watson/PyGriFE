/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
* This file is part of SEP
*
* Copyright 1993-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
* Copyright 2014 SEP developers
*
* SEP is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SEP is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with SEP.  If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef _MSC_VER
#define SEP_API __declspec(dllexport)
#else
#define SEP_API __attribute__((visibility("default")))
#endif

/* datatype codes */
#define SEP_TBYTE        11  /* 8-bit unsigned byte */
#define SEP_TINT         31  /* native int type */
#define SEP_TFLOAT       42
#define SEP_TDOUBLE      82

/* object & aperture flags */
#define SEP_OBJ_MERGED       0x0001  /* object is result of deblending */
#define SEP_OBJ_TRUNC        0x0002  /* object truncated at image boundary */
#define SEP_OBJ_DOVERFLOW    0x0004  /* not currently used, but could be */
#define SEP_OBJ_SINGU        0x0008  /* x,y fully correlated */
#define SEP_APER_TRUNC       0x0010
#define SEP_APER_HASMASKED   0x0020
#define SEP_APER_ALLMASKED   0x0040
#define SEP_APER_NONPOSITIVE 0x0080

/* noise_type values in sep_image */
#define SEP_NOISE_NONE   0
#define SEP_NOISE_STDDEV 1
#define SEP_NOISE_VAR    2

/* input flags for aperture photometry */
#define SEP_MASK_IGNORE      0x0004

/* threshold interpretation for sep_extract */
#define SEP_THRESH_REL 0  /* in units of standard deviations (sigma) */
#define SEP_THRESH_ABS 1  /* absolute data values */

/* filter types for sep_extract */
#define SEP_FILTER_CONV    0
#define SEP_FILTER_MATCHED 1

/* structs ------------------------------------------------------------------*/

/* sep_image
 *
 * Represents an image, including data, noise and mask arrays, and
 * gain.
 */
typedef struct {
  const void *data;  /* data array                */
  const void *noise; /* noise array (can be NULL) */
  const void *mask;  /* mask array (can be NULL)  */
  const void *segmap;/* segmap array (can be NULL)  */
  int dtype;         /* element type of image     */
  int ndtype;        /* element type of noise     */
  int mdtype;        /* element type of mask      */
  int sdtype;        /* element type of segmap    */
  int *segids;       /* array of unique ids in segmap  */
  int *idcounts;     /* counts of unique ids in segmap  */
  int numids;        /* number of unique ids in segmap */
  int w;             /* array width               */
  int h;             /* array height              */
  double noiseval;   /* scalar noise value; used only if noise == NULL */
  short noise_type;  /* interpretation of noise value                  */
  double gain;       /* (poisson counts / data unit)                   */
  double maskthresh; /* pixel considered masked if mask > maskthresh   */
} sep_image;

/* sep_catalog
 *
 * The result of sep_extract(). This is a struct of arrays. Each array has
 * one entry per detected object.
 */
typedef struct {
  int    nobj;                 /* number of objects (length of all arrays) */
  float	 *thresh;              /* threshold (ADU)                          */
  int	 *npix;                 /* # pixels extracted (size of pix array)   */
  int    *tnpix;                /* # pixels above thresh (unconvolved)      */
  int	 *xmin, *xmax;
  int    *ymin, *ymax;
  double *x, *y;                 /* barycenter (first moments)               */
  double *x2, *y2, *xy;		 /* second moments                           */
  double *errx2, *erry2, *errxy;      /* second moment errors            */
  float	 *a, *b, *theta;    /* ellipse parameters                       */
  float	 *cxx, *cyy, *cxy;  /* ellipse parameters (alternative)         */
  float	 *cflux;                /* total flux of pixels (convolved im)      */
  float	 *flux;      		 /* total flux of pixels (unconvolved)       */
  float  *cpeak;                /* peak intensity (ADU) (convolved)         */
  float  *peak;                 /* peak intensity (ADU) (unconvolved)       */
  int    *xcpeak, *ycpeak;       /* x, y coords of peak (convolved) pixel    */
  int    *xpeak, *ypeak;         /* x, y coords of peak (unconvolved) pixel  */
  short	 *flag;                 /* extraction flags                         */
  int    **pix;             /* array giving indicies of object's pixels in   */
                            /* image (linearly indexed). Length is `npix`.  */
                            /* (pointer to within the `objectspix` buffer)  */
  int    *objectspix;      /* buffer holding pixel indicies for all objects */
} sep_catalog;

/*-------------------------- source extraction ------------------------------*/

/* sep_extract()
 *
 * Extract sources from an image. Source Extractor defaults are shown
 * in [ ] above.
 *
 * Notes
 * -----
 * `dtype` and `ndtype` indicate the data type (float, int, double) of the
 * image and noise arrays, respectively.
 *
 * If `noise` is NULL, thresh is interpreted as an absolute threshold.
 * If `noise` is not null, thresh is interpreted as a relative threshold
 * (the absolute threshold will be thresh*noise[i,j]).
 *
 */
SEP_API int sep_extract(const sep_image *image,
		float thresh,         /* detection threshold           [1.5] */
                int thresh_type,      /* threshold units    [SEP_THRESH_REL] */
		int minarea,          /* minimum area in pixels          [5] */
		const float *conv,    /* convolution array (can be NULL)     */
                                      /*               [{1 2 1 2 4 2 1 2 1}] */
		int convw, int convh, /* w, h of convolution array     [3,3] */
                int filter_type,      /* convolution (0) or matched (1)  [0] */
		int deblend_nthresh,  /* deblending thresholds          [32] */
		double deblend_cont,  /* min. deblending contrast    [0.005] */
		int clean_flag,       /* perform cleaning?               [1] */
		double clean_param,   /* clean parameter               [1.0] */
                sep_catalog **catalog); /* OUTPUT catalog                    */



/* set and get the size of the pixel stack used in extract() */
SEP_API void sep_set_extract_pixstack(size_t val);
SEP_API size_t sep_get_extract_pixstack(void);

/* set and get the number of sub-objects limit when deblending in extract() */
SEP_API void sep_set_sub_object_limit(int val);
SEP_API int sep_get_sub_object_limit(void);

/* free memory associated with a catalog */
SEP_API void sep_catalog_free(sep_catalog *catalog);
