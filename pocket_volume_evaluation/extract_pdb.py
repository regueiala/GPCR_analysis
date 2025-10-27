import MDAnalysis as mda
import sys
# --- First trajectory (XTC) ---
topo = sys.argv[1]
traj = sys.argv[2]
u1 = mda.Universe(topo, traj)
selection = u1.select_atoms("backbone")

u1.trajectory[1]
output = topo.split('.pdb')[0]

print("Writing last frame from XTC to PDB...")
with mda.Writer(f"{output}_f1.pdb", selection.n_atoms) as W:
    W.write(selection)
