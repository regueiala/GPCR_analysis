def rmsd(u_complex, u_ref, sel):
    '''Calculate RMSD between a trajectory and a reference structure based on a selection.'''
    rmsd_analysis = rms.RMSD(u_complex, u_ref, select=sel) #, ref_frame=0
    rmsd_analysis.run()
    rmsd_df = pd.DataFrame(rmsd_analysis.results.rmsd, columns=["Frame", "Time", "RMSD"])
    return rmsd_df

def rmsd_sel ( df, trajs, names, res_sele, name_sele, ref_structure ,plot = False) : 
    '''Calculate RMSD for multiple trajectories and add results to a DataFrame. Optionally plot the RMSD over time with smoothing.'''
    # df contains Frames time System sim_name REP
    # Distinct colors per sim_name    
    plt.figure(figsize=(12, 7))
    avg_win = 100
    
    for traj, name in zip(trajs, names):
        rmsd_df = rmsd(traj, ref_structure, f"{res_sele}")

        mask = df["System"] == name
        n_rows = mask.sum()
        assert n_rows == len(rmsd_df)

        df.loc[mask, f"{name_sele}"] = rmsd_df["RMSD"].values
        
        if plot == True:
            local_df = df.loc[mask].copy()
            local_df["avg"] = gaussian_filter1d(local_df[f"{name_sele}"], avg_win)

            sim = local_df["sim_name"].iloc[0]
            color = palette.get(sim, "gray")

            sns.lineplot(data=local_df, x="$Time\;(\mu s)$", y=f"{name_sele}", color=color, lw=0.1, alpha=0.3)
            sns.lineplot(data=local_df, x="$Time\;(\mu s)$", y="avg", color=color, lw=2, label=name)

    if plot == True:
        plt.xlabel('$Time\;(\mu s)$')
        plt.ylabel('RMSD ($\\AA$)')
        plt.legend(title="Replicas", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.ylim(0, 15)
        plt.title(f'RMSD of {name_sele}')
        plt.tight_layout()
        plt.show()

    return df

def plot_rmsd_distribution (df,name_sele) :
    '''Plot the distribution of RMSD values for a given selection across different systems.'''
    plt.figure(figsize=(12, 6))
    # List of unique systems
    sys_list = np.unique(df["System"])
    
    # Plot a smoothed distribution instead of a simple histogram for better clarity
    for i, sys in enumerate(sys_list):
        data = df[df["System"] == sys][f"{name_sele}"].dropna()
        sim = df[df["System"] == sys]["sim_name"].iloc[0]
        color = palette.get(sim, "gray")  # gray if sim_name unknown
    
        if len(data) == 0:
            continue  # Ignore systems without RMSD data
        sns.kdeplot(data, label=sys, color=color, fill=True, alpha=0.4)
    
    # Axes and title
    plt.xlabel("RMSD values ($\\AA$)", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.title(f"Distribution of {name_sele} RMSD", fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.show()

def compute_volume (df, dat_list):
    '''Compute and plot volume over time from .dat files using the epock tool.'''
    plt.figure(figsize=(12, 6))
    avg_win = 100

    names = [
    "PTH1R-Gs1", "PTH1R-Gs2", "PTH1R-Gs3", "pth1r_pth_gs1", "pth1r_pth_gs2", "pth1r_pth_gs3",
    "PTH1R-Gq1", "PTH1R-Gq2", "PTH1R-Gq3", "pth1r_pthrp_gq1", "pth1r_pthrp_gq2", "pth1r_pthrp_gq3",
    "pth_8flq1", "pth_8flq2", "pth_8flq3", "pthrp_8flr1", "pthrp_8flr2", "pthrp_8flr3"]

    for name in names :
        for dat in dat_list:
            df_vol = pd.read_csv(dat, delim_whitespace=True)
            col = df_vol.columns[-1]
            col = col.replace("_rep", "")
            col = col.replace("pock_", "")
            if col.find("_Gs") != -1 : 
                col= col.replace("_Gs","_gs" )
            if col.find("_Gq") != -1 : 
                col= col.replace("_Gq","_gq" )
            if col in ['pth-8flq1','pth-8flq2','pth-8flq3','pthrp-8flr1','pthrp-8flr2','pthrp-8flr3'] : 
                col= col.replace("-8","_8" )
            if col in ['pth1r-Gq1','pth1r-Gq2','pth1r-Gq3','pth1r-Gs1','pth1r-Gs2','pth1r-Gs3'] : 
                col= col.replace("pth1r","PTH1R" )
            if name == col :
                mask = df["System"] == name
                n_rows = mask.sum()
                assert n_rows == df_vol.shape[0]
                
                # Inject the calculated RMSD into the main dataframe
                df.loc[mask, "volume"] = df_vol[df_vol.columns[-1]].values
            
                local_df = df.loc[mask].copy()
                local_df["avg"] = gaussian_filter1d(local_df["volume"], avg_win)
            
                # Determine color via sim_name
                sim = local_df["sim_name"].iloc[0]
                color = palette.get(sim, "gray")  # gray if sim_name unknown
            
                # Plot: raw RMSD
                sns.lineplot(data=local_df,x="$Time\;(\mu s)$",y="volume",color=color,lw=0.1,alpha=0.3)
            
                # Plot: smoothed RMSD, with replica legend
                sns.lineplot(data=local_df,x="$Time\;(\mu s)$",y="avg",color=color,lw=2,label=name)
    
    # Graph formatting
    plt.xlabel('$Time\;(\mu s)$')
    plt.ylabel('volume ($\mathrm{\AA}³$)')
    plt.legend(title="Replicas", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('volume = f(t)')
    plt.tight_layout()
    plt.show()
   
    return df 

def compute_tilt(traj, name, tmd_res_sel):
    '''Compute tilt angle of a helix in a trajectory over time along the z-axis.'''
    u = traj
    helix = u.select_atoms(tmd_res_sel)
    z_axis = np.array([0, 0, 1])
    v_membrane = z_axis / np.linalg.norm(z_axis)
    frames = []
    tilts = []

    for ts in u.trajectory:
        # Get helix vector (first–last Cα positions)
        v_helix = helix.positions[-1] - helix.positions[0]
        v_helix /= np.linalg.norm(v_helix)

        # Compute tilt angle (degrees)
        tilt_angle = np.degrees(np.arccos(np.dot(v_helix, v_membrane)))

        frames.append(ts.frame)
        tilts.append(tilt_angle)

    tilt_df = pd.DataFrame({"Frame": frames,"Tilt_Angle": tilts})

    return tilt_df

def compute_kappa(traj, part1_sel, part2_sel):
    '''Compute kappa angle (kink) between two parts of a molecule over time in a trajectory.'''
    kappa = []
    frames = []
    for ts in traj.trajectory:
        part1 = traj.select_atoms(part1_sel)
        part2  = traj.select_atoms(part2_sel)
        
        v1 = part1.positions[-1] - part1.positions[0]
        v1 /= np.linalg.norm(v1)
        
        v2 = part2.positions[-1] - part2.positions[0]
        v2 /= np.linalg.norm(v2)
        kappa_angle = np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)))
    
        frames.append(ts.frame)
        kappa.append(kappa_angle)
    kappa_df = pd.DataFrame({"Frame": frames,"kappa_Angle": kappa})
    return kappa_df
