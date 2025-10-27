#!/usr/bin/env python3

import sys
from MDAnalysis import Universe, Writer
from MDAnalysis.coordinates.XTC import XTCWriter



pdbfile_list = [
    "pth1r-Gq_rep1_10us_stripped_skip10_aligned.pdb",
    "pth1r_pthrp_Gq_rep1_10us_stripped_skip10_aligned.pdb",
    "pth1r-Gq_rep2_10us_stripped_skip10_aligned.pdb",
    "pth1r_pthrp_Gq_rep2_10us_stripped_skip10_aligned.pdb",
    "pth1r-Gq_rep3_10us_stripped_skip10_aligned.pdb",
    "pth1r_pthrp_Gq_rep3_10us_stripped_skip10_aligned.pdb",
    "pth1r-Gs_rep1_10us_stripped_skip10_aligned.pdb",
    "pth-8flq_rep1_10us_stripped_skip10_aligned.pdb",
    "pth1r-Gs_rep2_10us_stripped_skip10_aligned.pdb",
    "pth-8flq_rep2_10us_stripped_skip10_aligned.pdb",
    "pth1r-Gs_rep3_10us_stripped_skip10_aligned.pdb",
    "pth-8flq_rep3_10us_stripped_skip10_aligned.pdb",
    "pth1r_pth_Gs_rep1_10us_stripped_skip10_aligned.pdb",
    "pthrp-8flr_rep1_10us_stripped_skip10_aligned.pdb",
    "pth1r_pth_Gs_rep2_10us_stripped_skip10_aligned.pdb",
    "pthrp-8flr_rep2_10us_stripped_skip10_aligned.pdb",
    "pth1r_pth_Gs_rep3_10us_stripped_skip10_aligned.pdb",
    "pthrp-8flr_rep3_10us_stripped_skip10_aligned.pdb",
]

dcdfile_list = [
    "pth1r-Gq_rep1_10us_stripped_skip10_aligned.dcd",
    "pth1r_pthrp_Gq_rep1_10us_stripped_skip10_aligned.dcd",
    "pth1r-Gq_rep2_10us_stripped_skip10_aligned.dcd",
    "pth1r_pthrp_Gq_rep2_10us_stripped_skip10_aligned.dcd",
    "pth1r-Gq_rep3_10us_stripped_skip10_aligned.dcd",
    "pth1r_pthrp_Gq_rep3_10us_stripped_skip10_aligned.dcd",
    "pth1r-Gs_rep1_10us_stripped_skip10_aligned.dcd",
    "pth-8flq_rep1_10us_stripped_skip10_aligned.dcd",
    "pth1r-Gs_rep2_10us_stripped_skip10_aligned.dcd",
    "pth-8flq_rep2_10us_stripped_skip10_aligned.dcd",
    "pth1r-Gs_rep3_10us_stripped_skip10_aligned.dcd",
    "pth-8flq_rep3_10us_stripped_skip10_aligned.dcd",
    "pth1r_pth_Gs_rep1_10us_stripped_skip10_aligned.dcd",
    "pthrp-8flr_rep1_10us_stripped_skip10_aligned.dcd",
    "pth1r_pth_Gs_rep2_10us_stripped_skip10_aligned.dcd",
    "pthrp-8flr_rep2_10us_stripped_skip10_aligned.dcd",
    "pth1r_pth_Gs_rep3_10us_stripped_skip10_aligned.dcd",
    "pthrp-8flr_rep3_10us_stripped_skip10_aligned.dcd",
]


for pdb , dcd  in zip(pdbfile_list, dcdfile_list):
    xtc = dcd.replace(".dcd", ".xtc")
    xtc = xtc

    # Create a Universe with structure and trajectory
    u = Universe(pdb, dcd)
    selection = u.select_atoms("resindex 0:397")
    # Initialize XTC writer
    with XTCWriter(xtc, n_atoms=selection.n_atoms) as w:
        for ts in u.trajectory:
            w.write(selection)
            print(f"Converted frame {ts.frame}")


    tmp_pdb = pdb.replace("_aligned", "_aligned_sel")

    u.trajectory[1]           # se placer au premier frame
    selection.write(tmp_pdb)  # écrire uniquement ce frame


    print(f"Converted '{dcd}' --> '{xtc}'")
