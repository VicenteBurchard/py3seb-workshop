import os
import logging
from pathlib import Path
from typing import Union, List, Optional, Tuple

import numpy as np
import rasterio
import xarray as xr

from functions.eomaji.utils.general_utils import load_lut

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()],
)

lut = load_lut()


def get_biopar(
    connection, product: str, date: Union[str, List[str]], aoi: dict
) -> xr.Dataset:
    """
    Retrieve BIOPAR product from an openEO connection.

    Args:
        connection: openEO connection instance.
        product (str): Product type name (e.g. 'LAI', 'FAPAR').
        date (Union[str, List[str]]): Single date or date range [start, end].
        aoi (dict): GeoJSON-like polygon area of interest.

    Returns:
        xr.Dataset: BIOPAR dataset.
    """
    if isinstance(date, str):
        date = [date, date]

    biopar = connection.datacube_from_process(
        "BIOPAR",
        namespace=(
            "https://openeo.dataspace.copernicus.eu/openeo/1.1/processes/"
            "u:3e24e251-2e9a-438f-90a9-d4500e576574/BIOPAR"
        ),
        date=date,
        polygon=aoi,
        biopar_type=product,
    )
    return biopar


def calc_canopy(
    lai_path: Union[str, Path],
    worldcover_path: Union[str, Path],
    fg_path: Union[str, Path],
) -> None:
    """
    Estimate canopy height based on LAI, land cover, and green fraction.

    Args:
        lai_path (str | Path): Path to LAI GeoTIFF file.
        worldcover_path (str | Path): Path to land cover GeoTIFF file.
        fg_path (str | Path): Path to fraction green GeoTIFF file.

    Returns:
        None. Saves an H_C GeoTIFF to disk.
    """
    with rasterio.open(lai_path) as lai_src:
        lai = lai_src.read(1).astype(np.float32)
        profile = lai_src.profile

    with rasterio.open(worldcover_path) as worldcover_src:
        landcover = worldcover_src.read(1).astype(np.int32)
        landcover = 10 * (landcover // 10)

    with rasterio.open(fg_path) as fg_src:
        fg = fg_src.read(1).astype(np.float32)

    param_value = np.full_like(landcover, np.nan, dtype=np.float32)

    for lc_class in np.unique(landcover[~np.isnan(landcover)]):
        lc_pixels = np.where(landcover == lc_class)
        lc_index = lut[lut["landcover_class"] == lc_class].index[0]
        param_value[lc_pixels] = lut["veg_height"][lc_index]

        if lut["is_herbaceous"][lc_index] == 1:
            pai = lai / fg
            pai = pai[lc_pixels]
            param_value[lc_pixels] = 0.1 * param_value[lc_pixels] + 0.9 * param_value[
                lc_pixels
            ] * np.minimum((pai / lut["veg_height"][lc_index]) ** 3.0, 1.0)

    output_path = str(lai_path).replace("LAI", "H_C")
    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(param_value, 1)

    logging.info(f"Saved H_C to {output_path}")


def calc_fg(
    fapar_path: Union[str, Path], lai_path: Union[str, Path], sza_path: Union[str, Path]
) -> None:
    """
    Estimate green fraction (F_G) using FAPAR, LAI, and solar zenith angle.

    Args:
        fapar_path (Union[str, Path]): Path to FAPAR GeoTIFF.
        lai_path (Union[str, Path]): Path to LAI GeoTIFF.
        sza_path (Union[str, Path]): Path to solar zenith angle GeoTIFF.

    Returns:
        None. Saves F_G as GeoTIFF.
    """
    from pyTSEB import TSEB

    with rasterio.open(fapar_path) as fapar_src:
        fapar = fapar_src.read(1).astype(np.float32)
        profile = fapar_src.profile

    with rasterio.open(lai_path) as lai_src:
        lai = lai_src.read(1).astype(np.float32)

    with rasterio.open(sza_path) as sza_src:
        sza = sza_src.read(1).astype(np.float32)

    f_g = np.ones(lai.shape, dtype=np.float32)
    converged = np.zeros(lai.shape, dtype=bool)
    converged[np.logical_or(lai <= 0.2, fapar <= 0.1)] = True
    min_frac_green = 0.01

    for _ in range(50):
        f_g_old = f_g.copy()
        fipar = TSEB.calc_F_theta_campbell(
            sza[~converged], lai[~converged] / f_g[~converged], w_C=1, Omega0=1, x_LAD=1
        )
        f_g[~converged] = fapar[~converged] / fipar
        f_g = np.clip(f_g, min_frac_green, 1.0)
        converged = np.logical_or(np.isnan(f_g), np.abs(f_g - f_g_old) < 0.02)
        if np.all(converged):
            break

    profile.update(dtype=rasterio.float32, count=1)
    output_path = str(lai_path).replace("LAI", "F_G")
    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(f_g, 1)

    logging.info(f"Saved frac_green to {output_path}")


def split_tifs(nc_file: Union[str, Path], date_str: str) -> None:
    """
    Splits NetCDF bands into separate GeoTIFF files per variable.

    Args:
        nc_file (Union[str, Path]): Path to input NetCDF.
        date_str (str): Date string to include in output filenames.

    Returns:
        None. GeoTIFFs are written to disk.
    """
    out_dir = Path(nc_file).parent
    data = xr.open_dataset(nc_file)
    s2_bands = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]

    if Path(nc_file).stem == "s2_data":
        refl = data[s2_bands]
        refl = refl.rio.write_crs(data.crs.crs_wkt)
        output_file = out_dir / f"{date_str}_REFL.tif"
        refl.rio.to_raster(output_file)

    for var_name in data.data_vars:
        if var_name in s2_bands + ["crs"]:
            continue

        band = data[var_name]
        band.attrs.pop("grid_mapping", None)
        band = band.rio.write_crs(data.crs.crs_wkt)

        if var_name == "sunZenithAngles":
            output_file = out_dir / f"{date_str}_SZA.tif"
        else:
            output_file = out_dir / f"{date_str}_{var_name}.tif"

        band.rio.to_raster(output_file)
        logging.info(f"Saved {var_name} to {output_file}")


def _estimate_param_value(
    worldcover_path: Union[str, Path],
    lut: xr.Dataset,
    band: str,
    output_path: Union[str, Path],
) -> np.ndarray:
    """
    Estimate parameter value (e.g., leaf width) from land cover using a LUT.

    Args:
        worldcover_path (Union[str, Path]): Path to land cover GeoTIFF.
        lut (xr.Dataset): Lookup table mapping land cover class to param value.
        band (str): Parameter name (column in LUT).
        output_path (Union[str, Path]): Output GeoTIFF path.

    Returns:
        np.ndarray: Array of estimated values.
    """
    with rasterio.open(worldcover_path) as worldcover_src:
        landcover = worldcover_src.read(1).astype(np.int32)
        landcover = 10 * (landcover // 10)
        profile = worldcover_src.profile

    param_value = np.full(landcover.shape, np.nan, dtype=np.float32)

    for lc_class in np.unique(landcover[~np.isnan(landcover)]):
        lc_pixels = np.where(landcover == lc_class)
        lc_index = lut[lut["landcover_class"] == lc_class].index[0]
        param_value[lc_pixels] = lut[band][lc_index]

    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(param_value, 1)

    logging.info(f"Saved {band} to {output_path}")
    return param_value


def watercloud_model(param, a, b, c):
    result = a + b * (1.0 - np.exp(c * param))

    return result


def cab_to_vis_spectrum(
    cab,
    coeffs_wc_rho_vis=[0.14096573, -0.09648072, -0.06328343],
    coeffs_wc_tau_vis=[0.08543707, -0.08072709, -0.06562554],
):
    rho_leaf_vis = watercloud_model(cab, *coeffs_wc_rho_vis)
    tau_leaf_vis = watercloud_model(cab, *coeffs_wc_tau_vis)

    rho_leaf_vis = np.clip(rho_leaf_vis, 0, 1)
    tau_leaf_vis = np.clip(tau_leaf_vis, 0, 1)

    return rho_leaf_vis, tau_leaf_vis


def cw_to_nir_spectrum(
    cw,
    coeffs_wc_rho_nir=[0.38976106, -0.17260689, -65.7445699],
    coeffs_wc_tau_nir=[0.36187620, -0.18374560, -65.3125878],
):
    rho_leaf_nir = watercloud_model(cw, *coeffs_wc_rho_nir)
    tau_leaf_nir = watercloud_model(cw, *coeffs_wc_tau_nir)

    rho_leaf_nir = np.clip(rho_leaf_nir, 0, 1)
    tau_leaf_nir = np.clip(rho_leaf_nir, 0, 1)

    return rho_leaf_nir, tau_leaf_nir


def process_lai_to_vis(lai_path):
    """Processes a LAI raster file to generate visible spectrum reflectance and transmittance TIFFs."""
    try:
        if not Path(lai_path).exists():
            raise FileNotFoundError(f"LAI file not found: {lai_path}")

        with rasterio.open(lai_path) as src:
            meta = src.meta.copy()
            meta.update(dtype="float32")

            lai = src.read(1)
            cab = np.clip(np.array(lai), 0.0, 140.0)
            refl_vis, trans_vis = cab_to_vis_spectrum(cab)  # Function assumed to exist

            rho_vis_path = str(lai_path).replace("LAI", "RHO_VIS_C")
            tau_vis_path = str(lai_path).replace("LAI", "TAU_VIS_C")

            save_raster(rho_vis_path, refl_vis, meta)
            save_raster(tau_vis_path, trans_vis, meta)

            logging.info(f"Processed LAI to VIS: {rho_vis_path}, {tau_vis_path}")

    except Exception as e:
        logging.error(f"Error processing LAI: {e}")
        raise  # Ensure function fails on error


def process_cwc_to_nir(cw_path):
    """Processes a CWC raster file to generate NIR reflectance and transmittance TIFFs."""
    try:
        if not Path(cw_path).exists():
            raise FileNotFoundError(f"CWC file not found: {cw_path}")

        with rasterio.open(cw_path) as src:
            meta = src.meta.copy()
            meta.update(dtype="float32")

            cw = src.read(1)
            cw = np.clip(np.array(cw), 0.0, 0.1)
            refl_nir, trans_nir = cw_to_nir_spectrum(cw)  # Function assumed to exist

            rho_nir_path = str(cw_path).replace("CWC", "RHO_NIR_C")
            tau_nir_path = str(cw_path).replace("CWC", "TAU_NIR_C")

            save_raster(rho_nir_path, refl_nir, meta)
            save_raster(tau_nir_path, trans_nir, meta)

            logging.info(f"Processed CWC to NIR: {rho_nir_path}, {tau_nir_path}")

    except Exception as e:
        logging.error(f"Error processing CWC: {e}")
        raise  # Ensure function fails on error


def save_raster(output_path, data, meta):
    """Saves an array as a GeoTIFF using Rasterio."""
    try:
        with rasterio.open(output_path, "w", **meta) as dst:
            dst.write(data.astype("float32"), 1)
        logging.info(f"Saved raster: {output_path}")
    except Exception as e:
        logging.error(f"Failed to save raster {output_path}: {e}")
        raise  # Ensure function fails on error


def process_lai_and_cwc(lai_path, cw_path):
    process_lai_to_vis(lai_path)
    process_cwc_to_nir(cw_path)


def save_lat_lon_as_tifs(nc_file, out_dir, date):
    data = xr.open_dataset(nc_file)

    lat, lon = xr.broadcast(data["y"], data["x"])

    lat = lat.rio.write_crs(data.crs.crs_wkt)
    lat.rio.to_raster(f"{out_dir}/{date}_LAT.tif")

    lon = lon.rio.write_crs(data.crs.crs_wkt)
    lon.rio.to_raster(f"{out_dir}/{date}_LON.tif")


def split_datasets_to_tiffs(
    s2_path, s3_path, worldcover_path, date, out_dir: str | Path = None
):
    if out_dir:
        base_dir = Path(out_dir)
    else:
        base_dir = Path(s2_path).parent
    datestr = str(date).replace("-", "")

    split_tifs(s2_path, datestr)
    split_tifs(s3_path, datestr)

    lai_path = base_dir / f"{datestr}_LAI.tif"
    fapar_path = base_dir / f"{datestr}_FAPAR.tif"
    sza_path = base_dir / f"{datestr}_SZA.tif"
    fg_path = base_dir / f"{datestr}_F_G.tif"

    calc_fg(fapar_path, lai_path, sza_path)
    calc_canopy(lai_path, worldcover_path, fg_path)

    out_path = base_dir / f"{datestr}_W_C.tif"
    _ = _estimate_param_value(worldcover_path, lut, "veg_height_width_ratio", out_path)
    out_path = base_dir / f"{datestr}_LEAF_WIDTH.tif"
    _ = _estimate_param_value(worldcover_path, lut, "veg_leaf_width", out_path)

    cw_path = base_dir / f"{datestr}_CWC.tif"
    process_lai_and_cwc(lai_path, cw_path)

    return base_dir
