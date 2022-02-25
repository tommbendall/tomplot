from tomplot import convert_global_output

# ---------------------------------------------------------------------------- #
# Convert Bryan-Fritsch bubble logs
# ---------------------------------------------------------------------------- #

base_dir = 'results/consistent_moisture_paper'
test_dirs = ['bryan-fritsch-adv', 'bryan-fritsch-consist']
source_dirs = [f'{base_dir}/{opt}' for opt in test_dirs]
target_dirs = source_dirs.copy()

for source_dir, target_dir in zip(source_dirs, target_dirs):
    convert_global_output(target_dir, source_dir, mode='gungho_mass')

# ---------------------------------------------------------------------------- #
# Convert convergence tests
# ---------------------------------------------------------------------------- #

Lx = 2000.0
ns = [160,180,200,220,240]
dxs = [Lx / n for n in ns]
base_dir = 'results/consistent_moisture_paper'
test_dirs = [f'stretchy-conv-adv-{n}' for n in ns]
source_dirs = [f'{base_dir}/{opt}' for opt in test_dirs]
target_dir = source_dirs[0]

convert_global_output(target_dir, source_dirs, dxs=dxs, mode='transport_stats')