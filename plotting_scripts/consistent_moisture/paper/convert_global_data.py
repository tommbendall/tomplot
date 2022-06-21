from tomplot import convert_global_output

# ---------------------------------------------------------------------------- #
# Convert Bryan-Fritsch bubble logs
# ---------------------------------------------------------------------------- #

base_dir = 'results/consistent_moisture_paper'
test_dirs = ['bryan_fritsch-adv', 'bryan_fritsch-consist']
source_dirs = [f'{base_dir}/{opt}' for opt in test_dirs]
target_dirs = source_dirs.copy()

for source_dir, target_dir in zip(source_dirs, target_dirs):
    convert_global_output(target_dir, source_dir, mode='gungho_mass', dt=2.0)

# ---------------------------------------------------------------------------- #
# Convert convergence tests
# ---------------------------------------------------------------------------- #

Lx = 2000.0
ns = [120,140,160,180,200,220,240]
dxs = [Lx / n for n in ns]
base_dir = 'results/consistent_moisture_paper'
test_dirs = [f'stretchy-conv-adv-{n}' for n in ns]
source_dirs = [f'{base_dir}/{opt}' for opt in test_dirs]
target_dir = source_dirs[0]
run_params = {'dx':dxs}

convert_global_output(target_dir, source_dirs, mode='transport_stats', run_params=run_params)

test_dirs = [f'stretchy-conv-cons-{n}' for n in ns]
source_dirs = [f'{base_dir}/{opt}' for opt in test_dirs]
target_dir = source_dirs[0]

convert_global_output(target_dir, source_dirs, mode='transport_stats', run_params=run_params)

# ---------------------------------------------------------------------------- #
# Convert consistency tests
# ---------------------------------------------------------------------------- #

base_dir = 'results/consistent_moisture_paper'
test_dirs = ['stretchy-cons-cons-stretched', 'stretchy-cons-adv-stretched',
             'stretchy-cons-cons', 'stretchy-cons-adv']
source_dirs = [f'{base_dir}/{opt}' for opt in test_dirs]
target_dirs = source_dirs.copy()

for source_dir, target_dir in zip(source_dirs, target_dirs):
    convert_global_output(target_dir, source_dir, mode='transport_stats')