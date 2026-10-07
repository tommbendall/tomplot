"""
Routines for computing zonally-averaged fields on spherical domains.
"""

__all__ = ['zonal_average']

import numpy as np
import pandas as pd
from .data_extraction import reshape_gusto_data


def zonal_average(field_data, coords_lon, coords_lat, coords_height,
                  num_bins=18, lat_bins=None):
    """
    Computes a zonally-averaged field on a spherical domain, returning a
    lat-height slice.

    The field is averaged over longitude within a series of latitude bands,
    for each vertical (model) level. Data is assumed to be unstructured in
    the horizontal but may have varying heights per column (e.g. due to
    orography), so the returned heights are the mean height at each level.

    Args:
        field_data (`numpy.ndarray`): 1D unstructured array of field values.
        coords_lon (`numpy.ndarray`): 1D unstructured array of longitude
            coordinates.
        coords_lat (`numpy.ndarray`): 1D unstructured array of latitude
            coordinates.
        coords_height (`numpy.ndarray`): 1D unstructured array of height (or
            other vertical) coordinates, corresponding to each point in
            field_data.
        num_bins (int, optional): number of equal-width latitude bins to
            average over. Ignored if lat_bins is specified. Defaults to 18
            (10 degree bins).
        lat_bins (`numpy.ndarray`, optional): array of latitude bin edges to
            use instead of automatically generating equal-width bins.
            Defaults to None.

    Returns:
        tuple of `numpy.ndarray`: (zonal_mean, lat_bin_centres, level_heights).
            zonal_mean: 2D array of shape (num_lat_bins, num_levels) of the
                zonally-averaged field. Bins containing no data points are
                omitted entirely (rather than being filled with NaN), so
                num_lat_bins may be smaller than the requested num_bins.
            lat_bin_centres: 1D array of the centres of the latitude bins that
                contain data.
            level_heights: 1D array of the mean height at each level.
    """

    if len(np.shape(field_data)) != 1:
        raise ValueError('zonal_average: field_data must be 1D unstructured data')
    if len(np.shape(coords_lon)) != 1:
        raise ValueError('zonal_average: coords_lon must be 1D unstructured data')
    if len(np.shape(coords_lat)) != 1:
        raise ValueError('zonal_average: coords_lat must be 1D unstructured data')
    if len(np.shape(coords_height)) != 1:
        raise ValueError('zonal_average: coords_height must be 1D unstructured data')

    # ------------------------------------------------------------------------ #
    # Reshape data to be structured in the vertical
    # ------------------------------------------------------------------------ #
    field_2d, _, lat_2d, height_2d = reshape_gusto_data(
        field_data, coords_lon, coords_lat, coords_Z_3d=coords_height)

    num_levels = np.shape(field_2d)[1]

    # ------------------------------------------------------------------------ #
    # Determine latitude bins
    # ------------------------------------------------------------------------ #
    if lat_bins is None:
        # Guess whether latitudes are in degrees or radians, based on their
        # range (as done elsewhere in the repo, e.g. regrid_horizontal_slice)
        lat_range = np.max(coords_lat) - np.min(coords_lat)
        if lat_range > 10:
            min_lat, max_lat = -90.0, 90.0
        else:
            min_lat, max_lat = -np.pi/2, np.pi/2
        lat_bins = np.linspace(min_lat, max_lat, num_bins+1)

    lat_bin_centres = 0.5*(lat_bins[:-1] + lat_bins[1:])

    # ------------------------------------------------------------------------ #
    # Determine which bins actually contain data points
    # ------------------------------------------------------------------------ #
    # Latitude (and hence bin membership) is assumed to be the same for every
    # column regardless of level, so this only needs to be checked once
    col_bins = pd.cut(lat_2d[:, 0], bins=lat_bins, labels=False,
                      include_lowest=True)
    valid_bins = sorted(int(b) for b in pd.unique(col_bins) if not np.isnan(b))
    num_lat_bins = len(valid_bins)
    bin_position = {bin_idx: pos for pos, bin_idx in enumerate(valid_bins)}
    lat_bin_centres = lat_bin_centres[valid_bins]

    # ------------------------------------------------------------------------ #
    # Loop through levels, computing zonal mean for each
    # ------------------------------------------------------------------------ #
    zonal_mean = np.full((num_lat_bins, num_levels), np.nan)
    level_heights = np.zeros(num_levels)

    for lev_idx in range(num_levels):
        df = pd.DataFrame({'lat': lat_2d[:, lev_idx],
                           'field': field_2d[:, lev_idx]})
        df['lat_bin'] = pd.cut(df['lat'], bins=lat_bins, labels=False,
                               include_lowest=True)
        bin_means = df.groupby('lat_bin')['field'].mean()

        for bin_idx in bin_means.index:
            if int(bin_idx) in bin_position:
                zonal_mean[bin_position[int(bin_idx)], lev_idx] = bin_means[bin_idx]

        level_heights[lev_idx] = np.mean(height_2d[:, lev_idx])

    return zonal_mean, lat_bin_centres, level_heights
