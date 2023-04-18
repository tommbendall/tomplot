from tomplot import compute_field_measures, DiagnosticInfo
import numpy as np


target_dir = '/data/users/tbendall/results/moist_gw_all_convergence'
base_data_dir = '/data/users/tbendall/moist_gw_data_all'
true_source_dir = f'{base_data_dir}/D1200_P1200'
variable_name = 'theta_e'
prognostic_variables = ['theta', 'exner', 'm_v']
measures = 'L2_error'

Lx = 300000.0
extrusion_details = {'domain':'plane', 'extrusion':'linear',
                     'zmin':0.0, 'zmax':10000, 'topological_dimension':3}

# First number dynamics, second is physics
# res_pairs = [[100, 100], [100, 200], [100, 600], [120, 120], [120, 60],
#              [150, 150], [150, 300], [150, 600], [150, 75], [200, 100],
#              [200, 200], [200, 400], [200, 600], [300, 150], [300, 300],
#              [300, 600], [400, 200], [400, 400], [600, 100], [600, 150],
#              [600, 200], [600, 300], [600, 60], [600, 75], [60, 120],
#              [60, 60], [60, 600], [75, 150], [75, 600], [75, 75]]

res_pairs = [[60, 60], [75, 75], [100, 100], [120, 120], [150, 150],
             [200, 200], [300, 300], [400, 400], [600,600],
             [100, 50], [100, 200], [100, 300], [100, 400], [100, 600],
             [200, 50], [200, 100], [200, 400], [200, 600], [200, 800],
             [300, 60], [300, 75], [300, 100], [300, 150], [300, 600], [300, 900], [300, 1200],
             [400, 50], [400, 80], [400, 100], [400, 200], [400, 800], [400, 1200],
             [600, 120], [600, 150], [600, 200], [600, 300] ]


coarse_dirs = [f'D{nxs[0]}_P{nxs[1]}' for nxs in res_pairs]
dyn_dxs = [Lx/nxs[0] for nxs in res_pairs]
phys_dxs = [Lx/nxs[1] for nxs in res_pairs]
run_params = {'dyn_dx':dyn_dxs, 'phys_dx':phys_dxs}
source_dirs = [f'{base_data_dir}/{coarse_dir}' for coarse_dir in coarse_dirs]

def theta_e_evaluator(data):

    Lv = 2.501e6
    cp = 1005.0

    m_v = data['m_v']
    exner = data['exner']
    theta = data['theta']

    # Unfortunately have to map exner to Wth as diagnostic exner_in_wth
    # is on the physics mesh and not the dynamics mesh
    data_shape = np.shape(exner)
    exner_in_wth = np.zeros((data_shape[0]+1, data_shape[1]))

    # Assume uniform extrusion
    exner_in_wth[0,:] = 2.0*exner[0,:] - exner[1,:]
    exner_in_wth[-1,:] = 2.0*exner[-1,:] - exner[-2,:]
    for j in range(1,data_shape[0]):
        exner_in_wth[j,:] = 0.5*(exner[j,:]+exner[j-1,:])

    T = theta * exner_in_wth
    exp_arg = Lv * m_v / (cp * T)
    theta_e = theta * np.exp(exp_arg)

    return theta_e

diagnostic = DiagnosticInfo(variable_name, prognostic_variables, theta_e_evaluator)

compute_field_measures(target_dir, true_source_dir, source_dirs, diagnostic,
                       measures, extrusion_details, run_params=run_params)