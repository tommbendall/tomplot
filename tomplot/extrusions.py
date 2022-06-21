"""
Makes some LFRic extrusions.
"""

import numpy as np

def generate_extrusion(extrusion_details, vert_placement, num_points):

    extrusion = extrusion_details['extrusion']
    num_levels = num_points if vert_placement == 'full_levels' else num_points + 1

    if  extrusion == 'linear':

        full_level_heights = np.linspace(extrusion_details['zmin'],
                                         extrusion_details['zmax'],
                                         num_levels)

    else:
        raise NotImplementedError(f'Extrusion {extrusion} not implemented')

    
    if vert_placement == 'full_levels':

        heights_1d = full_level_heights

    elif vert_placement == 'half_levels':

        heights_1d = [0.5*(full_level_heights[i]+full_level_heights[i+1])
                      for i in range(len(full_level_heights)-1)]

    else:
        raise ValueError(f'Vertical placement {vert_placement} not recognised')

    return np.array(heights_1d)