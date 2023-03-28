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
    
    elif extrusion == 'um_L70_50t_20s_80km':

        assert num_levels == 71, 'Need to have 71 full levels for 70-layer UM extrusion'
        
        full_level_heights = extrusion_details['zmax']*np.array([
            0.0000000,  0.0002500,  0.0006667,  0.0012500,
            0.0020000,  0.0029167,  0.0040000,  0.0052500,
            0.0066667,  0.0082500,  0.0100000,  0.0119167,
            0.0140000,  0.0162500,  0.0186667,  0.0212500,
            0.0240000,  0.0269167,  0.0300000,  0.0332500,
            0.0366667,  0.0402500,  0.0440000,  0.0479167,
            0.0520000,  0.0562500,  0.0606667,  0.0652500,
            0.0700000,  0.0749167,  0.0800000,  0.0852500,
            0.0906668,  0.0962505,  0.1020017,  0.1079213,
            0.1140113,  0.1202745,  0.1267154,  0.1333406,
            0.1401592,  0.1471838,  0.1544313,  0.1619238,
            0.1696895,  0.1777643,  0.1861929,  0.1950307,
            0.2043451,  0.2142178,  0.2247466,  0.2360480,
            0.2482597,  0.2615432,  0.2760868,  0.2921094,
            0.3098631,  0.3296378,  0.3517651,  0.3766222,
            0.4046373,  0.4362943,  0.4721379,  0.5127798,
            0.5589045,  0.6112759,  0.6707432,  0.7382500,
            0.8148403,  0.9016668,  1.0000000 ])

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