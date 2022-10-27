from tomplot import convert_global_output
import numpy as np

# base_dir = '/data/users/tbendall/results/pdc_idealised_paper'

# # ---------------------------------------------------------------------------- #
# # Convert convergence tests
# # ---------------------------------------------------------------------------- #

# radius = 6371229.0
# n_pairs = [[48,24],[24,12],[64,32],[32,16],[96,48],
#            [96,12],[96,16],[96,24],[96,32],]
# fine_dxs = [2*np.pi*radius / n[0] for n in n_pairs]
# coarse_dxs = [2*np.pi*radius / n[1] for n in n_pairs]
# run_params = {'fine_dx':fine_dxs, 'coarse_dx':coarse_dxs}
# test_dirs = [f'stretchy_sphere_conv_C{n[0]}_C{n[1]}' for n in n_pairs]
# source_dirs = [f'{base_dir}/{opt}' for opt in test_dirs]
# target_dir = f'{base_dir}/stretchy_sphere_conv_all'

# convert_global_output(target_dir, source_dirs, mode='transport_stats', run_params=run_params)

# test_dirs = [f'stretchy_sphere_conv_adv_C{n[0]}_C{n[1]}' for n in n_pairs]
# source_dirs = [f'{base_dir}/{opt}' for opt in test_dirs]
# target_dir = f'{base_dir}/stretchy_sphere_conv_adv_all'

# convert_global_output(target_dir, source_dirs, mode='transport_stats', run_params=run_params)


# # ---------------------------------------------------------------------------- #
# # Convert consistency tests
# # ---------------------------------------------------------------------------- #

# source_dir = f'{base_dir}/stretchy_sphere_cons'
# target_dir = source_dir

# convert_global_output(target_dir, source_dir, mode='transport_stats')

# ---------------------------------------------------------------------------- #
# Convert consistency tests
# ---------------------------------------------------------------------------- #

source_dir = '/data/users/tbendall/results/pdc_schar/schar_pdc'
target_dir = source_dir

convert_global_output(target_dir, source_dir, mode='transport_stats')