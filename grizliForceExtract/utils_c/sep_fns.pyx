import numpy as np
cimport numpy as np
from libc cimport limits
from libc.math cimport sqrt
cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.version cimport PY_MAJOR_VERSION

np.import_array()  # To access the numpy C-API.

# -----------------------------------------------------------------------------
# Definitions from the SEP C library

# macro definitions from sep.h
DEF SEP_TBYTE = 11
DEF SEP_TINT = 31
DEF SEP_TFLOAT = 42
DEF SEP_TDOUBLE = 82

# input flag values (C macros)
DEF SEP_NOISE_NONE = 0
DEF SEP_NOISE_STDDEV = 1
DEF SEP_NOISE_VAR = 2

# filter types for sep_extract
DEF SEP_FILTER_CONV = 0
DEF SEP_FILTER_MATCHED = 1

# Threshold types
DEF SEP_THRESH_REL = 0
DEF SEP_THRESH_ABS = 1

# input flags for aperture photometry
DEF SEP_MASK_IGNORE = 0x0004

# Output flag values accessible from python
OBJ_MERGED = np.short(0x0001)
OBJ_TRUNC = np.short(0x0002)
OBJ_DOVERFLOW = np.short(0x0004)
OBJ_SINGU = np.short(0x0008)
APER_TRUNC = np.short(0x0010)
APER_HASMASKED = np.short(0x0020)
APER_ALLMASKED = np.short(0x0040)
APER_NONPOSITIVE = np.short(0x0080)

# macro defintion from sepcore.h
# This is not part of the SEP API, but we pull it out because want to
# explicitly detect memory errors so that we can raise MemoryError().
DEF MEMORY_ALLOC_ERROR = 1

# header definitions
cdef extern from "sep.h":

    ctypedef struct sep_image:
        const void *data
        const void *noise
        const void *mask
        const void *segmap
        int dtype
        int ndtype
        int mdtype
        int sdtype
        int *segids
        int *idcounts
        int numids
        int w
        int h
        double noiseval
        short noise_type
        double gain
        double maskthresh

    ctypedef struct sep_bkg:
        int w
        int h
        float globalback
        float globalrms

    ctypedef struct sep_catalog:
        int    nobj
        float  *thresh
        int    *npix
        int    *tnpix
        int    *xmin
        int    *xmax
        int    *ymin
        int    *ymax
        double *x
        double *y
        double *x2
        double *y2
        double *xy
        double *errx2
        double *erry2
        double *errxy
        float  *a
        float  *b
        float  *theta
        float  *cxx
        float  *cyy
        float  *cxy
        float  *cflux
        float  *flux
        float  *cpeak
        float  *peak
        int    *xcpeak
        int    *ycpeak
        int    *xpeak
        int    *ypeak
        short  *flag
        int    **pix
        int    *objectspix

    int sep_extract(const sep_image *image,
                    float thresh,
                    int thresh_type,
                    int minarea,
                    float *conv,
                    int convw, int convh,
                    int filter_type,
                    int deblend_nthresh,
                    double deblend_cont,
                    int clean_flag,
                    double clean_param,
                    sep_catalog **catalog)

    void sep_catalog_free(sep_catalog *catalog)

    int sep_flux_radius(const sep_image *image,
                        double x, double y, double rmax, int id, int subpix,
                        short inflag,
                        double *fluxtot, double *fluxfrac, int n,
                        double *r, short *flag)

    int sep_kron_radius(const sep_image *image,
                        double x, double y, double cxx, double cyy,
                        double cxy, double r, int id,
                        double *kronrad, short *flag)

    void sep_set_extract_pixstack(size_t val)
    size_t sep_get_extract_pixstack()

    void sep_set_sub_object_limit(int val)
    int sep_get_sub_object_limit()

    void sep_get_errmsg(int status, char *errtext)
    void sep_get_errdetail(char *errtext)

# -----------------------------------------------------------------------------
# Utility functions

cdef int _get_sep_dtype(dtype) except -1:
    """Convert a numpy dtype to the corresponding SEP dtype integer code."""
    if not dtype.isnative:
        raise ValueError(
            "Input array with dtype '{0}' has non-native byte order. "
            "Only native byte order arrays are supported. "
            "To change the byte order of the array 'data', do "
            "'data = data.byteswap().newbyteorder()'".format(dtype))
    t = dtype.type
    if t is np.single:
        return SEP_TFLOAT
    elif t is np.bool_ or t is np.ubyte:
        return SEP_TBYTE
    elif dtype == np.double:
        return SEP_TDOUBLE
    elif dtype == np.intc:
        return SEP_TINT
    raise ValueError('input array dtype not supported: {0}'.format(dtype))


cdef int _check_array_get_dims(np.ndarray arr, int *w, int *h) except -1:
    """Check some things about an array and return dimensions"""

    # Raise an informative message if array is not C-contiguous
    if not arr.flags["C_CONTIGUOUS"]:
        raise ValueError("array is not C-contiguous")

    # Check that there are exactly 2 dimensions
    if arr.ndim != 2:
        raise ValueError("array must be 2-d")

    # ensure that arr dimensions are not too large for C ints.
    if arr.shape[0] <= <Py_ssize_t> limits.INT_MAX:
        h[0] = arr.shape[0]
    else:
        raise ValueError("array height  ({0:d}) greater than INT_MAX ({1:d})"
                         .format(arr.shape[0], limits.INT_MAX))
    if arr.shape[1] <= <Py_ssize_t> limits.INT_MAX:
       w[0] = arr.shape[1]
    else:
        raise ValueError("array width ({0:d}) greater than INT_MAX ({1:d})"
                         .format(arr.shape[1], limits.INT_MAX))
    return 0

cdef int _assert_ok(int status) except -1:
    """Get the SEP error message corresponding to status code"""
    cdef char *errmsg
    cdef char *errdetail

    if status == 0:
        return 0

    # First check if we have an out-of-memory error, so we don't try to
    # allocate more memory to hold the error message.
    if status == MEMORY_ALLOC_ERROR:
        raise MemoryError

    # Otherwise, get error message.
    errmsg = <char *>PyMem_Malloc(61 * sizeof(char))
    sep_get_errmsg(status, errmsg)
    pyerrmsg = <bytes> errmsg
    PyMem_Free(errmsg)

    # Get error detail.
    errdetail = <char *>PyMem_Malloc(512 * sizeof(char))
    sep_get_errdetail(errdetail)
    pyerrdetail = <bytes> errdetail
    PyMem_Free(errdetail)

    # If error detail is present, append it to the message.
    if pyerrdetail != b"":
        pyerrmsg = pyerrmsg + b": " + pyerrdetail

    # Convert string to unicode if on python 3
    if PY_MAJOR_VERSION == 3:
        msg = pyerrmsg.decode()
    else:
        msg = pyerrmsg

    raise Exception(msg)


cdef int _parse_arrays(np.ndarray data, err, var, mask, segmap,
                       sep_image *im) except -1:
    """Helper function for functions accepting data, error, mask & segmap arrays.
    Fills in an sep_image struct."""

    cdef int ew, eh, mw, mh, sw, sh
    cdef np.uint8_t[:,:] buf, ebuf, mbuf, sbuf
    cdef np.uint8_t[:] idbuf, countbuf

    # Clear im fields we might not touch (everything besides data, dtype, w, h)
    im.noise = NULL
    im.mask = NULL
    im.segmap = NULL
    im.numids = 0
    im.ndtype = 0
    im.mdtype = 0
    im.noiseval = 0.0
    im.noise_type = SEP_NOISE_NONE
    im.gain = 0.0
    im.maskthresh = 0.0

    # Get main image info
    _check_array_get_dims(data, &(im.w), &(im.h))
    im.dtype = _get_sep_dtype(data.dtype)
    buf = data.view(dtype=np.uint8)
    im.data = <void*>&buf[0, 0]

    # Check if noise is error or variance.
    noise = None  # will point to either error or variance.
    if err is not None:
        if var is not None:
            raise ValueError("Cannot specify both err and var")
        noise = err
        im.noise_type = SEP_NOISE_STDDEV
    elif var is not None:
        noise = var
        im.noise_type = SEP_NOISE_VAR

    # parse noise
    if noise is None:
        im.noise = NULL
        im.noise_type = SEP_NOISE_NONE
        im.noiseval = 0.0
    elif isinstance(noise, np.ndarray):
        if noise.ndim == 0:
            im.noise = NULL
            im.noiseval = noise
        elif noise.ndim == 2:
            _check_array_get_dims(noise, &ew, &eh)
            if ew != im.w or eh != im.h:
                raise ValueError("size of error/variance array must match"
                                 " data")
            im.ndtype = _get_sep_dtype(noise.dtype)
            ebuf = noise.view(dtype=np.uint8)
            im.noise = <void*>&ebuf[0, 0]
        else:
            raise ValueError("error/variance array must be 0-d or 2-d")
    else:
        im.noise = NULL
        im.noiseval = noise

    # Optional input: mask
    if mask is None:
        im.mask = NULL
    else:
        _check_array_get_dims(mask, &mw, &mh)
        if mw != im.w or mh != im.h:
            raise ValueError("size of mask array must match data")
        im.mdtype = _get_sep_dtype(mask.dtype)
        mbuf = mask.view(dtype=np.uint8)
        im.mask = <void*>&mbuf[0, 0]

    # Optional input: segmap
    if segmap is None:
        im.segmap = NULL
    else:
        _check_array_get_dims(segmap, &sw, &sh)
        if sw != im.w or sh != im.h:
            raise ValueError("size of segmap array must match data")
        im.sdtype = _get_sep_dtype(segmap.dtype)
        # ids_buf = np.ascontiguousarray(np.unique(segmap), dtype=ctypes.c_int)
        # print (ids_buf)
        # idbuf = np.unique(segmap).view(dtype=np.uint8)
        # im.segids = <void*>&idbuf[0]
        # im.segids = np.unique(segmap).ascontiguousarray(dtype=ctypes.c_int)
        # from cython.view cimport array as cy_array
        ids, counts = np.unique(segmap, return_counts=True)
        print (ids, counts)
        segids = np.ascontiguousarray(ids[1:], dtype=np.dtype("i"))
        idcounts = np.ascontiguousarray(counts[1:], dtype=np.dtype("i"))
        idbuf = segids.view(dtype=np.uint8)
        countbuf = idcounts.view(dtype=np.uint8)
        im.segids = <int*>&idbuf[0]
        im.idcounts = <int*>&countbuf[0]

        
        # print (segids)
        # # im.segids = 
        # sg = cy_array(shape=len(segids), itemsize=sizeof(int), format="i")
        # for s in segids:
        #     sg[i] = s
        # im.segids = <int*>PyMem_Malloc(len(segids)*sizeof(int))
        # print (im.segids)
    # im.segids = np.zeros([], dtype=int)
        # <int*>segids.data
        im.numids = len(ids[1:])
        print (f"Unique ids: {im.numids}.\n")
        sbuf = segmap.view(dtype=np.uint8)
        im.segmap = <void*>&sbuf[0, 0]

# -----------------------------------------------------------------------------
# Source Extraction

# This needs to match the result from extract
cdef packed struct Object:
    np.float64_t thresh
    np.int_t npix
    np.int_t tnpix
    np.int_t xmin
    np.int_t xmax
    np.int_t ymin
    np.int_t ymax
    np.float64_t x
    np.float64_t y
    np.float64_t x2
    np.float64_t y2
    np.float64_t xy
    np.float64_t a
    np.float64_t b
    np.float64_t theta
    np.float64_t cxx
    np.float64_t cyy
    np.float64_t cxy
    np.float64_t cflux
    np.float64_t flux
    np.float64_t cpeak
    np.float64_t peak
    np.float64_t errx2
    np.float64_t erry2
    np.float64_t errxy
    np.int_t xcpeak
    np.int_t ycpeak
    np.int_t xpeak
    np.int_t ypeak
    np.int_t flag

default_kernel = np.array([[1.0, 2.0, 1.0],
                           [2.0, 4.0, 2.0],
                           [1.0, 2.0, 1.0]], dtype=np.float32)

@cython.boundscheck(True)
def extract(np.ndarray data not None, float thresh, err=None, var=None,
            gain=None, np.ndarray mask=None, double maskthresh=0.0,
            int minarea=5,
            np.ndarray filter_kernel=default_kernel, filter_type='matched',
            int deblend_nthresh=32, double deblend_cont=0.005,
            bint clean=True, double clean_param=1.0,
            segmentation_map=None):
    """extract(data, thresh, err=None, mask=None, minarea=5,
               filter_kernel=default_kernel, filter_type='matched',
               deblend_nthresh=32, deblend_cont=0.005, clean=True,
               clean_param=1.0, segmentation_map=False)

    Extract sources from an image.

    Parameters
    ----------
    data : `~numpy.ndarray`
        Data array (2-d).
    thresh : float
        Threshold pixel value for detection. If an ``err`` or ``var`` array
        is not given, this is interpreted as an absolute threshold. If ``err``
        or ``var`` is given, this is interpreted as a relative threshold: the
        absolute threshold at pixel (j, i) will be ``thresh * err[j, i]`` or
        ``thresh * sqrt(var[j, i])``.
    err, var : float or `~numpy.ndarray`, optional
        Error *or* variance (specify at most one). This can be used to
        specify a pixel-by-pixel detection threshold; see "thresh" argument.
    gain : float, optional
        Conversion factor between data array units and poisson counts. This
        does not affect detection; it is used only in calculating Poisson
        noise contribution to uncertainty parameters such as ``errx2``. If
        not given, no Poisson noise will be added.
    mask : `~numpy.ndarray`, optional
        Mask array. ``True`` values, or numeric values greater than
        ``maskthresh``, are considered masked. Masking a pixel is equivalent
        to setting data to zero and noise (if present) to infinity.
    maskthresh : float, optional
        Threshold for a pixel to be masked. Default is ``0.0``.
    minarea : int, optional
        Minimum number of pixels required for an object. Default is 5.
    filter_kernel : `~numpy.ndarray` or None, optional
        Filter kernel used for on-the-fly filtering (used to
        enhance detection). Default is a 3x3 array:
        [[1,2,1], [2,4,2], [1,2,1]]. Set to ``None`` to skip
        convolution.
    filter_type : {'matched', 'conv'}, optional
        Filter treatment. This affects filtering behavior when a noise
        array is supplied. ``'matched'`` (default) accounts for
        pixel-to-pixel noise in the filter kernel. ``'conv'`` is
        simple convolution of the data array, ignoring pixel-to-pixel
        noise across the kernel.  ``'matched'`` should yield better
        detection of faint sources in areas of rapidly varying noise
        (such as found in coadded images made from semi-overlapping
        exposures).  The two options are equivalent when noise is
        constant.
    deblend_nthresh : int, optional
        Number of thresholds used for object deblending. Default is 32.
    deblend_cont : float, optional
        Minimum contrast ratio used for object deblending. Default is 0.005.
        To entirely disable deblending, set to 1.0.
    clean : bool, optional
        Perform cleaning? Default is True.
    clean_param : float, optional
        Cleaning parameter (see SExtractor manual). Default is 1.0.
    segmentation_map : bool, optional
        If True, also return a "segmentation map" giving the member
        pixels of each object. Default is False.

    Returns
    -------
    objects : `~numpy.ndarray`
        Extracted object parameters (structured array). Available fields are:

        * ``thresh`` (float) Threshold at object location.
        * ``npix`` (int) Number of pixels belonging to the object.
        * ``tnpix`` (int) Number of pixels above threshold (unconvolved data).
        * ``xmin``, ``xmax`` (int) Minimum, maximum x coordinates of pixels.
        * ``ymin``, ``ymax`` (int) Minimum, maximum y coordinates of pixels.
        * ``x``, ``y`` (float) object barycenter (first moments).
        * ``x2``, ``y2``, ``xy`` (float) Second moments.
        * ``errx2``, ``erry2``, ``errxy`` (float) Second moment errors.
          Note that these will be zero if error is not given.
        * ``a``, ``b``, ``theta`` (float) Ellipse parameters, scaled as
          described by Section 8.4.2 in "The Source Extractor Guide" or
          Section 10.1.5-6 of v2.13 of SExtractor's User Manual.
        * ``cxx``, ``cyy``, ``cxy`` (float) Alternative ellipse parameters.
        * ``cflux`` (float) Sum of member pixels in convolved data.
        * ``flux`` (float) Sum of member pixels in unconvolved data.
        * ``cpeak`` (float) Peak value in convolved data.
        * ``peak`` (float) Peak value in unconvolved data.
        * ``xcpeak``, ``ycpeak`` (int) Coordinate of convolved peak pixel.
        * ``xpeak``, ``ypeak`` (int) Coordinate of unconvolved peak pixel.
        * ``flag`` (int) Extraction flags.

    segmap : `~numpy.ndarray`, optional
        Array of integers with same shape as data. Pixels not belonging to
        any object have value 0. All pixels belonging to the ``i``-th object
        (e.g., ``objects[i]``) have value ``i+1``. Only returned if
        ``segmentation_map=True``.
    """

    cdef int kernelw, kernelh, status, i, j
    cdef int filter_typecode, thresh_type
    cdef sep_catalog *catalog = NULL
    cdef np.ndarray[Object] result
    cdef float[:, :] kernelflt
    cdef float *kernelptr
    cdef np.int32_t[:, :] segmap_buf
    cdef np.int32_t *segmap_ptr
    cdef int *objpix
    cdef sep_image im

    if type(segmentation_map) is np.ndarray:
        _parse_arrays(data, err, var, mask, segmentation_map, &im)
    else:
        _parse_arrays(data, err, var, mask, None, &im)
    im.maskthresh = maskthresh
    if gain is not None:
        im.gain = gain

    # Parse filter input
    if filter_kernel is None:
        kernelptr = NULL
        kernelw = 0
        kernelh = 0
    else:
        kernelflt = filter_kernel.astype(np.float32)
        kernelptr = &kernelflt[0, 0]
        kernelw = kernelflt.shape[1]
        kernelh = kernelflt.shape[0]

    if filter_type == 'matched':
        filter_typecode = SEP_FILTER_MATCHED
    elif filter_type == 'conv':
        filter_typecode = SEP_FILTER_CONV
    else:
        raise ValueError("unknown filter_type: {!r}".format(filter_type))

    # If image has error info, the threshold is relative, otherwise
    # it is absolute.
    if im.noise_type == SEP_NOISE_NONE:
        thresh_type = SEP_THRESH_ABS
    else:
        thresh_type = SEP_THRESH_REL

    status = sep_extract(&im,
                         thresh, thresh_type, minarea,
                         kernelptr, kernelw, kernelh, filter_typecode,
                         deblend_nthresh, deblend_cont, clean, clean_param,
                         &catalog)
    _assert_ok(status)

    # Allocate result record array and fill it
    result = np.empty(catalog.nobj,
                      dtype=np.dtype([('thresh', np.float64),
                                      ('npix', np.int_),
                                      ('tnpix', np.int_),
                                      ('xmin', np.int_),
                                      ('xmax', np.int_),
                                      ('ymin', np.int_),
                                      ('ymax', np.int_),
                                      ('x', np.float64),
                                      ('y', np.float64),
                                      ('x2', np.float64),
                                      ('y2', np.float64),
                                      ('xy', np.float64),
                                      ('errx2', np.float64),
                                      ('erry2', np.float64),
                                      ('errxy', np.float64),
                                      ('a', np.float64),
                                      ('b', np.float64),
                                      ('theta', np.float64),
                                      ('cxx', np.float64),
                                      ('cyy', np.float64),
                                      ('cxy', np.float64),
                                      ('cflux', np.float64),
                                      ('flux', np.float64),
                                      ('cpeak', np.float64),
                                      ('peak', np.float64),
                                      ('xcpeak', np.int_),
                                      ('ycpeak', np.int_),
                                      ('xpeak', np.int_),
                                      ('ypeak', np.int_),
                                      ('flag', np.int_)]))

    for i in range(catalog.nobj):
        result['thresh'][i] = catalog.thresh[i]
        result['npix'][i] = catalog.npix[i]
        result['tnpix'][i] = catalog.tnpix[i]
        result['xmin'][i] = catalog.xmin[i]
        result['xmax'][i] = catalog.xmax[i]
        result['ymin'][i] = catalog.ymin[i]
        result['ymax'][i] = catalog.ymax[i]
        result['x'][i] = catalog.x[i]
        result['y'][i] = catalog.y[i]
        result['x2'][i] = catalog.x2[i]
        result['y2'][i] = catalog.y2[i]
        result['xy'][i] = catalog.xy[i]
        result['errx2'][i] = catalog.errx2[i]
        result['erry2'][i] = catalog.erry2[i]
        result['errxy'][i] = catalog.errxy[i]
        result['a'][i] = catalog.a[i]
        result['b'][i] = catalog.b[i]
        result['theta'][i] = catalog.theta[i]
        result['cxx'][i] = catalog.cxx[i]
        result['cyy'][i] = catalog.cyy[i]
        result['cxy'][i] = catalog.cxy[i]
        result['cflux'][i] = catalog.cflux[i]
        result['flux'][i] = catalog.flux[i]
        result['cpeak'][i] = catalog.cpeak[i]
        result['peak'][i] = catalog.peak[i]
        result['xcpeak'][i] = catalog.xcpeak[i]
        result['ycpeak'][i] = catalog.ycpeak[i]
        result['xpeak'][i] = catalog.xpeak[i]
        result['ypeak'][i] = catalog.ypeak[i]
        result['flag'][i] = catalog.flag[i]

    # construct a segmentation map, if it was requested.
    if type(segmentation_map) is np.ndarray or segmentation_map:
        # Note: We have to write out `(data.shape[0], data.shape[1])` because
        # because Cython turns `data.shape` later into an int pointer when
        # the function argument is typed as np.ndarray.
        segmap = np.zeros((data.shape[0], data.shape[1]), dtype=np.int32)
        segmap_buf = segmap
        segmap_ptr = &segmap_buf[0, 0]
        for i in range(catalog.nobj):
            objpix = catalog.pix[i]
            for j in range(catalog.npix[i]):
                segmap_ptr[objpix[j]] = i + 1

    # Free the C catalog
    sep_catalog_free(catalog)

    if type(segmentation_map) is np.ndarray or segmentation_map:
        return result, segmap
    else:
        return result

# -----------------------------------------------------------------------------
# Utility functions

def set_extract_pixstack(size_t size):
    """set_extract_pixstack(size)

    Set the size in pixels of the internal pixel buffer used in extract().

    The current value can be retrieved with get_extract_pixstack. The
    initial default is 300000.
    """
    sep_set_extract_pixstack(size)

def get_extract_pixstack():
    """get_extract_pixstack()

    Get the size in pixels of the internal pixel buffer used in extract().
    """
    return sep_get_extract_pixstack()


def set_sub_object_limit(int limit):
    """set_sub_object_limit(limit)

    Set the limit on the number of sub-objects when deblending in extract().

    The current value can be retrieved with get_sub_object_limit. The
    initial default is 1024.
    """
    sep_set_sub_object_limit(limit)

def get_sub_object_limit():
    """get_sub_object_limit()

    Get the limit on the number of sub-objects when deblending in extract().
    """
    return sep_get_sub_object_limit()