"""
Routines for interpreting information to do with a spatial domain, and applying
it to a figure.
"""

import numpy as np

__all__ = ["apply_gusto_domain"]


def apply_gusto_domain(ax, dataset, slice_along=None, units=None, xlabel=True,
                       ylabel=True, xlabelpad=-10, ylabelpad=-10):
    """
    Edits the limits and applies labels to an axes, based on metadata relating
    to the domain that is stored in the Gusto netCDF data file.

    This is not intended to be used for making neat plots, but is useful for
    quickly making plots.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object to alter.
        dataset (:class:`Dataset`): the netCDF dataset containing the data.
        slice_along (str, optional): which coordinate axis to take this plot
            represents a slice along. Valid values are 'x', 'y', 'lon', 'lat',
            'z', 'alpha' or 'beta'. Defaults to None.
        xlabel (str, optional): x-axis label, to override that generated by the
            domain data. Defaults to True, in which case it is generated from
            the domain data.
        ylabel (str, optional): y-axis label, to override that generated by the
            domain data. Defaults to True, in which case it is generated from
            the domain data.
        xlabelpad (float, optional): padding to apply to the x-axis label.
            Defaults to None.
        ylabelpad (float, optional): padding to apply to the y-axis label.
            Defaults to None.
    """

    assert slice_along in [None, 'x', 'y', 'lon', 'lat', 'alpha', 'beta', 'z'], \
        f'slice_along variable must correspond to a coordinate. {slice_along}' \
        + ' is not a valid value'

    stored_xlabel = xlabel
    stored_ylabel = ylabel

    # ------------------------------------------------------------------------ #
    # Checks on units argument and set default for this domain
    # ------------------------------------------------------------------------ #
    domain = dataset.variables['domain_type'][:]

    if domain == 'spherical_shell':
        assert units in [None, 'deg', 'rad'], 'apply_gusto_domain: for a ' \
            + f'spherical shell domain units must be "deg" or "rad" not {units}'
        if units is None:
            units = 'deg'
    elif domain == 'vertical_slice':
        assert units in [None, 'm', 'km'], 'apply_gusto_domain: for a ' \
            + f'vertical slice domain units must be "m" or "km" not {units}'
        if units is None:
            units = 'km'
    elif domain == 'extruded_spherical_shell':
        assert units in [None, 'deg', 'rad'], 'apply_gusto_domain: for an ' \
            + f'extruded sphere, units must be "deg" or "rad" not {units}'
        if units is None:
            units = 'deg'
    else:
        raise NotImplementedError(f'apply_gusto_domain: domain {domain} '
                                  +' either not implemented or recognised')

    # ------------------------------------------------------------------------ #
    # Work out domain extent based on the domain type
    # ------------------------------------------------------------------------ #
    if domain == 'spherical_shell':
        if units == 'deg':
            xlims = [-180, 180]
            ylims = [-90, 90]
            xlabel = r'$\lambda \ / $ deg'
            ylabel = r'$\phi \ / $ deg'
        elif units == 'rad':
            xlims = [-np.pi, np.pi]
            ylims = [-np.pi/2, np.pi/2]
            xlabel = r'$\lambda \ / $ rad'
            ylabel = r'$\phi \ / $ rad'

    elif domain == 'vertical_slice':
        if units == 'm':
            unit_factor = 1
            xlabel = r'$x \ / $ m'
            ylabel = r'$z \ / $ m'
        elif units == 'km':
            unit_factor = 0.001
            xlabel = r'$x \ / $ km'
            ylabel = r'$z \ / $ km'

        xlims = [0, unit_factor*dataset['domain_extent_x'][:].data]
        ylims = [0, unit_factor*dataset['domain_extent_z'][:].data]

    elif domain == 'extruded_spherical_shell':
        if units in 'deg':
            xlims = [-180, 180]
            ylims = [-90, 90]
            xlabel = r'$\lambda \ / $ deg'
            ylabel = r'$\phi \ / $ deg'
        elif units == 'rad':
            xlims = [-np.pi, np.pi]
            ylims = [-np.pi/2, np.pi/2]
            xlabel = r'$\lambda \ / $ rad'
            ylabel = r'$\phi \ / $ rad'

        zlabel = r'$z \ / $ m'
        # round this value as it may be a bit odd
        zlims = [0, round(1.0*dataset['domain_extent_z'][:].data, 3)]

        if slice_along is None:
            raise ValueError('apply_gusto_domain: slice_along variable needs '
                             + 'providing for extruded spherical shell data')

        if slice_along in ['y', 'lat']:
            ylabel = zlabel
            ylims = zlims
        elif slice_along in ['x', 'lon']:
            xlabel = ylabel
            xlims = ylims
            ylabel = zlabel
            ylims = zlims
        elif slice_along == 'z':
            # Everything already correct, pass
            pass

        else:
            raise ValueError('apply_gusto_domain: slice_along variable '
                             + f'{slice_along} not valid for extruded sphere')

    else:
        raise NotImplementedError(f'apply_gusto_domain: domain {domain} '
                                  +' either not implemented or recognised')

    # ------------------------------------------------------------------------ #
    # Apply domain limits and labels
    # ------------------------------------------------------------------------ #

    ax.set_xlim(xlims)
    ax.set_xticks(xlims)
    ax.set_xticklabels(xlims)
    ax.set_ylim(ylims)
    ax.set_yticks(ylims)
    ax.set_yticklabels(ylims)

    # Apply labels. Possibly these have been specified by the user.
    if stored_xlabel == True:
        ax.set_xlabel(xlabel, labelpad=xlabelpad)
    elif stored_xlabel == False:
        pass
    else:
        ax.set_xlabel(stored_xlabel, labelpad=xlabelpad)

    if stored_ylabel == True:
        ax.set_ylabel(ylabel, labelpad=ylabelpad)
    elif stored_ylabel == False:
        pass
    else:
        ax.set_ylabel(stored_ylabel, labelpad=ylabelpad)