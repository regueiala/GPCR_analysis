import os
import numpy as np

# Get all .dat files in the current directory
dat_files = [f for f in os.listdir() if f.endswith('.dat')]

for dat in dat_files:
    x_list, y_list, z_list = [], [], []
    with open(dat, 'r') as file:
        for line in file:
            # Skip lines starting with '0'
            if not line.startswith('0'):
                parts = line.strip().split()
                if len(parts) >= 4:  # Ensure there are at least 4 columns
                    x_list.append(float(parts[1]))
                    y_list.append(float(parts[2]))
                    z_list.append(float(parts[3]))

# Compute mean coordinates with offsets
x_mean = np.mean(x_list) - 2
y_mean = np.mean(y_list)
z_mean = np.mean(z_list) + 5


for dat in dat_files:
    # Create output .cfg file
    cfg_filename = f"{dat.split('_10us')[0]}_config.cfg"
    with open(cfg_filename, 'w') as out_file:
        out_file.write(f"""[DEFAULT]
grid_spacing    = 0.5

[pock_{dat.split('_10us')[0]}]
include_sphere  = {x_mean} {y_mean} {z_mean} 12
include_sphere  = 5.889100154576584 68.55149318469131 30.771209711176557 10
""")

    print(f"{dat.split('_10us')[0]}: {x_mean}, {y_mean}, {z_mean}")
