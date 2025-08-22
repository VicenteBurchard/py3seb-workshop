import numpy as np
import skimage.morphology as morph
from tseb_input_creator import gdal_utils as gu
import logging

log = logging.getLogger(__name__)


def confidence_level(values, t1, t2):
    scale = 1. / (t2 - t1)
    prob = np.full_like(values, np.nan)
    valid = np.isfinite(values)

    # probability = 0
    case = np.logical_and(valid, values < t1)
    prob[case] = 0
    # probability = 1
    case = np.logical_and(valid, values > t2)
    prob[case] = 1
    # Intermediate probabilities
    case = np.logical_and.reduce((valid,
                                  values >= t1,
                                  values <= t2))
    prob[case] = (values[case] - t1) * scale
    return prob


def sen2cor_cloud_prob(rhos):
    # Step 1a - Brightness thresholds on red (Band 4)
    s1a = confidence_level(rhos[2], 0.06, 0.25)
    # Step 1b – Normalized Difference Snow Index (NDSI)
    ndsi = _veg_index(rhos[1], rhos[8])
    s1b = confidence_level(ndsi, -0.24, -0.16)
    cloud = s1a * s1b
    return cloud


def sen2cor_mask(rhos, threshold=0, buffer_clouds=3):
    rhos = rhos.astype(float)
    no_data = rhos == -9999
    rhos[no_data] = np.nan
    rhos[~no_data] = rhos[~no_data] / 1e4
    cloud_prob = sen2cor_cloud_prob(rhos)
    clear = cloud_prob <= threshold
    
    if buffer_clouds > 0:
        se = morph.disk(buffer_clouds)
        clear = morph.binary_erosion(clear, se)

    clear = np.logical_and(np.isfinite(cloud_prob), clear)
    return clear


def reclassify_qai(lc_array,
                   best_aerosol=False,
                   best_illumination=True,
                   best_wvp=False):
    """Reclassifies FORCE output QAI such that all the pixels which should be masked
    have a value of zero, and pixels which should not be masked have a value of 1.

    Parameters
    ----------
    lc_array : array
        Numpy array of the mask.
    best_aerosol : bool
        True if want to discard interpolated aerosols.
    best_illumination : bool
        True if want to discard medium illumination conditions.
    best_wvp : bool
        True if want to discard interpolated water vapour.

    Returns
    -------
    valid : array
        Boolean array with valid pixels (True or 1 is valid observation)
    """
    # Bitwise operators
    # if not valid (bit 0) == 1 then no valid
    no_valid = (lc_array & (1 << 0)) > 0
    # if cloud confidence (bits 1 & 2) != 00 then clouds
    clouds = (lc_array & ((1 << 1) | (1 << 2))) > 0
    # if shadows (bit 3) == 1 then shadows
    shadow = (lc_array & (1 << 3)) > 0
    # if snow (bit 4 ) == 1 then snow
    snow = (lc_array & (1 << 4)) > 0
    # if water (bit 5) == 1 then water
    water = (lc_array & (1 << 5)) > 0
    # if high aerosol load (bit 7) == 1 then high_aerosol
    high_aerosol = (lc_array & (1 << 7)) > 0
    # if interpolated aerosol load (bit 6) == 1 int_aerosol
    int_aerosol = (lc_array & (1 << 6)) > 0
    # if poor/shadow incidence angle (bit 12) == 1 then poor_incidence
    poor_incidence = (lc_array & (1 << 12)) > 0
    # if medium incidence angle (bit 11) == 1 then med_incidence
    med_incidence = (lc_array & (1 << 11)) > 0
    # if fill wvp (bit 14) == 1 then wvp
    wvp = (lc_array & (1 << 14)) > 0
    # Final mask is without clouds nor snow, with sensor observations and succesful AOT
    valid = np.logical_and.reduce((~no_valid, ~clouds, ~shadow, ~snow,
                                   ~water, ~high_aerosol, ~poor_incidence))
    if best_aerosol is True:
        valid = np.logical_and(valid, ~int_aerosol)
    if best_illumination is True:
        valid = np.logical_and(valid, ~med_incidence)
    if best_wvp is True:
        valid = np.logical_and(valid, ~wvp)
    valid[lc_array == 1] = False
    return valid


def s2l2a_mask(qai_file):
    date, level, satellite, product = qai_file.stem.split("_")
    scene = f"{date}_{level}_{satellite}"
    log.info(f'Creating binary QA mask from {qai_file}')
    prj, gt, *_ = gu.raster_info(qai_file)
    valid = gu.get_raster_data(qai_file)
    clear_qai = reclassify_qai(valid, best_illumination=False)
    boa_file = qai_file.parent / f"{scene}_BOA.tif"
    image_array = gu.get_raster_data(boa_file)
    clear_sen2cor = sen2cor_mask(image_array)
    valid = np.logical_and(clear_qai, clear_sen2cor)
    output_file = qai_file.parent / f"{scene}_MASK.tif"
    gu.save_image(valid, gt, prj, output_file)
    return valid


def reclassify_syn_flags(lc_array,
                         best_aerosol=False,
                         all_observations=False,
                         buffer_clouds=3):
    """Reclassifies FORCE output QAI such that all the pixels which should be masked
    have a value of zero, and pixels which should not be masked have a value of 1.

    Parameters
    ----------
    lc_array : array
        Numpy array of the mask.
    best_aerosol : bool
        True if want to discard interpolated aerosols.
    best_illumination : bool
        True if want to discard medium illumination conditions.
    best_wvp : bool
        True if want to discard interpolated water vapour.

    Returns
    -------
    valid : array
        Boolean array with valid pixels (True or 1 is valid observation)
    """
    lc_array = lc_array.astype(int)
    no_data = lc_array < 0
    # Bitwise operators
    # if clouod (bit 0) then cloud
    clouds = (lc_array & (1 << 0)) > 0
    # areas that are not clouds, but features similar to clouds (bit 1)
    snow_risk = (lc_array & (1 << 1)) > 0
    if buffer_clouds > 0:
        se = morph.disk(buffer_clouds)
        clouds = morph.binary_dilation(clouds, se)
        snow_risk = morph.binary_dilation(snow_risk, se)

    # Land (bit 4)
    land = (lc_array & (1 << 4)) > 0
    # No OLCI radiances (bit 5)
    no_olci = (lc_array & (1 << 5)) > 0
    # No SLSTR nadir radiance (bit 6)
    no_sln = (lc_array & (1 << 6)) > 0
    # No SLSTR oblique radiance (bit 7)
    no_slo = (lc_array & (1 << 7)) > 0
    # No aerosol filled (bit 11)
    aod_filled = (lc_array & (1 << 11)) > 0
    # Aerosol success (bit 12)
    aod_success = (lc_array & (1 << 12)) > 0
    aod_valid = np.logical_or(aod_filled, aod_success)
    # Final mask is without clouds nor snow,
    # with sensor observations and succesful AOT
    valid = np.logical_and.reduce((~no_data, ~clouds, ~snow_risk, land,
                                   ~no_olci, ~no_sln, aod_valid))

    if best_aerosol:
        valid = np.logical_and(valid, ~aod_filled)

    if all_observations:
        valid = np.logical_and(valid, ~no_slo)

    return valid


def reclassify_qa_pixel(lc_array,
                        buffer_clouds=1,
                        buffer_shadow=0):
    """Reclassifies fmask output map such that all the pixels which should be
    masked have a value of zero, and pixels which should not be masked have a
    value of 1.
    Parameters:
        qa_pixel_file (str):
            Path to fmask output map.
        buffer_clouds (int):
            Buffer in pixel size for clouds.
        buffer_shadow (int):
            Buffer in pixel size for shadows.
    """

    # if clouds (bit 3) and low/medium/high probability (bit 8 and 9) then clouds
    clouds = ((lc_array & (1 << 3)) > 1) & ((lc_array & ((1 << 8) | (1 << 9))) > 1)
    # if cirrus (bit 2) and low/medium/high probability shadows (bit 14 and 15) then shadows
    cirrus = ((lc_array & (1 << 2)) > 1) & ((lc_array & ((1 << 14) | (1 << 15))) > 1)
    # if shadows (bit 4) and low/medium/high probability shadows (bit 10 and 11) then shadows
    shadow = ((lc_array & (1 << 4)) > 1) & ((lc_array & ((1 << 10) | (1 << 11))) > 1)

    if buffer_clouds > 0:
        se = morph.disk(buffer_clouds)
        clouds = morph.binary_dilation(clouds, se)
        cirrus = morph.binary_dilation(cirrus, se)

    # Find shadow flags
    if buffer_shadow > 0:
        se = morph.disk(buffer_shadow)
        shadow = morph.binary_dilation(shadow, se)

    # Final mask is without clouds nor snow, with sensor observations and succesful AOT
    valid = np.logical_and.reduce((~clouds, ~shadow, ~cirrus))

    valid[lc_array == 1] = False
    return valid


def _veg_index(band1, band2):
    vi = np.full_like(band1, np.nan)
    valid = np.logical_and.reduce((np.isfinite(band1),
                                   np.isfinite(band2),
                                   band1 + band2 != 0))
    vi[valid] = (band1[valid] - band2[valid]) \
                / (band1[valid] + band2[valid])
    return vi


