#!/usr/bin/env bash
set -euo pipefail

pdbfile_list=(
    pth1r-Gq_rep1_10us_stripped_skip10_aligned_sel.pdb
    pth1r_pthrp_Gq_rep1_10us_stripped_skip10_aligned_sel.pdb
    pth1r-Gq_rep2_10us_stripped_skip10_aligned_sel.pdb
    pth1r_pthrp_Gq_rep2_10us_stripped_skip10_aligned_sel.pdb
    pth1r-Gq_rep3_10us_stripped_skip10_aligned_sel.pdb
    pth1r_pthrp_Gq_rep3_10us_stripped_skip10_aligned_sel.pdb
    pth1r-Gs_rep1_10us_stripped_skip10_aligned_sel.pdb
    pth-8flq_rep1_10us_stripped_skip10_aligned_sel.pdb
    pth1r-Gs_rep2_10us_stripped_skip10_aligned_sel.pdb
    pth-8flq_rep2_10us_stripped_skip10_aligned_sel.pdb
    pth1r-Gs_rep3_10us_stripped_skip10_aligned_sel.pdb
    pth-8flq_rep3_10us_stripped_skip10_aligned_sel.pdb
    pth1r_pth_Gs_rep1_10us_stripped_skip10_aligned_sel.pdb
    pthrp-8flr_rep1_10us_stripped_skip10_aligned_sel.pdb
    pth1r_pth_Gs_rep2_10us_stripped_skip10_aligned_sel.pdb
    pthrp-8flr_rep2_10us_stripped_skip10_aligned_sel.pdb
    pth1r_pth_Gs_rep3_10us_stripped_skip10_aligned_sel.pdb
    pthrp-8flr_rep3_10us_stripped_skip10_aligned_sel.pdb
)

xtcfile_list=(
    pth1r-Gq_rep1_10us_stripped_skip10_aligned.xtc
    pth1r_pthrp_Gq_rep1_10us_stripped_skip10_aligned.xtc
    pth1r-Gq_rep2_10us_stripped_skip10_aligned.xtc
    pth1r_pthrp_Gq_rep2_10us_stripped_skip10_aligned.xtc
    pth1r-Gq_rep3_10us_stripped_skip10_aligned.xtc
    pth1r_pthrp_Gq_rep3_10us_stripped_skip10_aligned.xtc
    pth1r-Gs_rep1_10us_stripped_skip10_aligned.xtc
    pth-8flq_rep1_10us_stripped_skip10_aligned.xtc
    pth1r-Gs_rep2_10us_stripped_skip10_aligned.xtc
    pth-8flq_rep2_10us_stripped_skip10_aligned.xtc
    pth1r-Gs_rep3_10us_stripped_skip10_aligned.xtc
    pth-8flq_rep3_10us_stripped_skip10_aligned.xtc
    pth1r_pth_Gs_rep1_10us_stripped_skip10_aligned.xtc
    pthrp-8flr_rep1_10us_stripped_skip10_aligned.xtc
    pth1r_pth_Gs_rep2_10us_stripped_skip10_aligned.xtc
    pthrp-8flr_rep2_10us_stripped_skip10_aligned.xtc
    pth1r_pth_Gs_rep3_10us_stripped_skip10_aligned.xtc
    pthrp-8flr_rep3_10us_stripped_skip10_aligned.xtc
)

config_list=(
    pth1r-Gq_rep1_config.cfg
    pth1r_pthrp_Gq_rep1_config.cfg
    pth1r-Gq_rep2_config.cfg
    pth1r_pthrp_Gq_rep2_config.cfg
    pth1r-Gq_rep3_config.cfg
    pth1r_pthrp_Gq_rep3_config.cfg
    pth1r-Gs_rep1_config.cfg
    pth-8flq_rep1_config.cfg
    pth1r-Gs_rep2_config.cfg
    pth-8flq_rep2_config.cfg
    pth1r-Gs_rep3_config.cfg
    pth-8flq_rep3_config.cfg
    pth1r_pth_Gs_rep1_config.cfg
    pthrp-8flr_rep1_config.cfg
    pth1r_pth_Gs_rep2_config.cfg
    pthrp-8flr_rep2_config.cfg
    pth1r_pth_Gs_rep3_config.cfg
    pthrp-8flr_rep3_config.cfg
)

# Limit tsp to 20 concurrent jobs
tsp -S 20

# Create central folder for output volumes
mkdir -p volume_analysis

# Loop through all files
for ((i=0; i<${#pdbfile_list[@]}; i++)); do
    pdb=${pdbfile_list[$i]}
    traj=${xtcfile_list[$i]}
    config=${config_list[$i]}

    folder="analysis_${pdb%.pdb}"
    mkdir -p "$folder"

    outfile="${config%.cfg}_volume.dat"

    echo "Queueing: $config → $outfile"

    # Submit job safely using tsp
    tsp bash -c "
        cd '$folder' &&
        epock -s '../$pdb' -f '../$traj' -c '../$config' --ox -o '$outfile' -v &&
        cp '$outfile' '../volume_analysis/'
    "
done
