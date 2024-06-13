# Carlos Angel
# Python 3.10.12
# Date Created: 2023-11-27

# Testing CTCF 

## pcSAC modules
import globals 
import utils
import reconstruct


## General modules
import sys
import numpy as np
import scipy as scp
import pandas as pd
from pathlib import Path
import os



## Hi-C modules
import cooler

## Plotting modules
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

# Global variables
microC = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/data/MicroC/"
microC_path= microC + "GSE130275_mESC_WT_combined_1.3B.mcool"

rcmc_path = "/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_02_28_spliC/2023_02_28_v0_dataGathering/RCMC_data/GSM6281849_RCMC_BR1_merged_allCap_WT_mm39.merged.50.mcool"
clamp_n = True


def plot_matrices(LR,ctcf_matrix,merged_matrix,HR,start,end):
    extents = (start, end, end, start)
    norm_balanced = LogNorm(vmax=0.12,vmin=0.00007)
    matrices = [LR,ctcf_matrix,merged_matrix,HR]

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    titles = ['Hi-C LR', 'CTCF', 'Hi-C LR + CTCF (Input)', 'RCMC']

    for i, ax in enumerate(axes.flat):
    
        cax = ax.matshow(
            matrices[i],
            cmap='rocket_r',
            norm=norm_balanced,
            extent=extents
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(titles[i])

    fig.colorbar(cax, ax=axes.ravel().tolist(), fraction=0.036, pad=0.03)

    plt.show()

def check_shapes(matrix_ctcf,matrix_LR):
    print("Matrices have different shapes")
    if abs(matrix_ctcf.shape[0] - matrix_LR.shape[0]) > 1:
        print("Matrices have more than 1 row difference, cannot merge")
        return None
    elif abs(matrix_ctcf.shape[0] - matrix_LR.shape[0]) == 1:
        print("Matrices have 1 row difference, adding row and column of zeros")
        # Add a row of zeros
        if matrix_ctcf.shape[0] < matrix_LR.shape[0]:
            zeros_row = np.zeros((1, matrix_ctcf.shape[1]))
            matrix_ctcf = np.vstack([matrix_ctcf, zeros_row])
            
            zeros_column = np.zeros((matrix_ctcf.shape[0], 1))
            matrix_ctcf = np.hstack([matrix_ctcf, zeros_column])
            
    return matrix_ctcf, matrix_LR

def calculate_distance(x,y,resolution):
    return abs(x-y) * resolution

def power_law_ctcf(mat_val, distance):
    if distance == 0:
        return mat_val 
    return mat_val * (distance**(-1))

def merge_ctcf_and_LR(matrix_ctcf,matrix_LR,resolution,option="add"):
    '''Resulution must be the LR resolution'''

    if matrix_ctcf.shape != matrix_LR.shape:
        matrix_ctcf,matrix_LR = check_shapes(matrix_ctcf,matrix_LR)  
    else:
        print("Matrices have the same shape")

    print("Merging information")
    matrix_LR = np.nan_to_num(matrix_LR)
    matrix_ctcf = np.nan_to_num(matrix_ctcf)

    print("Scaling CTCF matrix...")

    LR_max = np.max(matrix_LR)

    matrix_ctcf[matrix_ctcf > 0] = LR_max
    print(f"LR matrix shape: {matrix_LR.shape}")
    print(f"CTCF matrix shape: {matrix_ctcf.shape}")

    scaled_ctcf = np.zeros(matrix_ctcf.shape)

    print(scaled_ctcf.shape)

    for i in range(matrix_ctcf.shape[0]):
        for j in range(matrix_ctcf.shape[1]):
            distance = calculate_distance(i,j,resolution)
            scaled_ctcf[i][j] = power_law_ctcf(matrix_ctcf[i][j],distance)
    factor = LR_max * 10000
    scaled_ctcf = scaled_ctcf * factor

    if option == "add":
        matrix = np.add(matrix_LR,scaled_ctcf)
    else:
        matrix = np.add(matrix_LR,scaled_ctcf)
        overlap_indices = (matrix_LR != 0) & (scaled_ctcf != 0)
        matrix[overlap_indices] = scaled_ctcf[overlap_indices]

    matrix = np.nan_to_num(matrix)

    return matrix.astype(float),scaled_ctcf.astype(float)       

def load_goal_matrix(chr,start,end,resolution,type="RCMC"):
    if resolution != 3200 and resolution != 800:
        print("Resolution not available")
        return None
    
    else:
        resolution = str("::resolutions/"+str(resolution))

    if type == "RCMC" or type == "rcmc":
        print("Loading RCMC matrix")
        print("Caution: RCMC matrix is in mm39 coordinates")
        rcmc_file = rcmc_path
        clr = cooler.Cooler(str(rcmc_file)+resolution)
        region = f"{chr}:{start:,}-{end:,}"
        M = clr.matrix().fetch(region)
        M = np.nan_to_num(M)

        return M
    else:
        print("Our goal is to use RCMC matrices")
        return None

def load_interactions(type,chr,start,end,resolution):
    if resolution == 3200:
        resolution = str("::/12")
    if resolution == 800:
        resolution = str("::/14")
    if resolution == 1600:
        resolution = str("::/13")
    
    else:
        print("Resolution not available")
        return None

    if type == "MicroC" or type == "microC" or type == "microc":
        hic_file = microC_path
        
        clr = cooler.Cooler(str(hic_file)+resolution)
        region = f"{chr}:{start:,}-{end:,}"
        P = clr.matrix(balance=True).fetch(region)
        P = np.nan_to_num(P)

        return P
    
    else:
        print("Currently only microC is available")
        return None


def load_ctcf_sites(file,chr):
    print("only resolutions 3200,1600and 800 are available")

    ctcf_files_dir = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/data/ChIP_seq/"
    ctcf_sites = pd.read_csv(ctcf_files_dir + file,
                            sep="\t", header=0)
    ctcf_info = ctcf_sites[(ctcf_sites['seqnames1'] == chr) & (ctcf_sites['seqnames2'] == chr)]
    ctcf_info = ctcf_info[['start1','end2','ctcfSignal']]
    
    
    return ctcf_info


'''ctcf_df : dataframe with ctcf sites result from load_ctcf_sites
    start : start of the region of interest
    end : end of the region of interest
    resolution : resolution of the LR matrix'''

def ctcf_matrix_range(ctcf_df,start,end,resolution):
     
    ctcf_range = ctcf_df[(ctcf_df['start1'] <= end) & (ctcf_df['start1'] >= start) & (ctcf_df['end2'] >= start) & (ctcf_df['end2'] <= end)]
    ctcf_range['start1'] = ctcf_range['start1']-1
    ctcf_range['ctcfSignal'] = 1
    matrix = reconstruct.coordinates2matrix(ctcf_range ,resolution,start,end)

    return matrix


def interaction_file(interact_matrix,hic_file,hic_res,size,d_factor=0,l_res=2):
    filename = os.path.splitext(os.path.basename(hic_file))[0]
    pwd = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning"
    LR_F = int((hic_res * l_res))
    if(clamp_n != True):
        output_dir = os.path.join(pwd, 
                    "results", 
                    "input_matrices", 
                    "{}_data_{}_resolution_{}_Mb_{}_downSample_{}_kb_LR".format(
                        filename, 
                        hic_res, 
                        size, 
                        d_factor, 
                        LR_F ))
    if(clamp_n == True):
        norm_type = "clamp"
        output_dir = os.path.join(pwd, 
                    "results_RCMC",norm_type, 
                    "input_matrices", 
                    "{}_data_{}_resolution_{}_Mb_{}_downSample_{}_kb_LR".format(
                        filename, 
                        hic_res, 
                        size, 
                        d_factor, 
                        LR_F ))
        
    os.makedirs(output_dir, exist_ok=True)

    rm = (int(hic_res)/100)/2 #radius in amstrongs

    pn = "interaction_matrix_"+filename+".txt"
    p = os.path.join(output_dir,pn)
    p2n = "int_mat_seg_len_"+filename+".txt"
    p2 = os.path.join(output_dir,p2n)
    p3n = "whole_matrix_"+filename+".txt" #whole matrix is normalized
    p3 = os.path.join(output_dir,p3n)

    Path(p).touch(exist_ok=True)

    #l = [str(i*l_res)+" "+str((i+1)*l_res-1)+" "+str(j*l_res)+" "+str((j+1)*l_res-1)+" "+str((interact_matrix[4])[i][j])+"\n" for j in range(1,(interact_matrix[4]).shape[0]) for i in range(j)]
    #Giving only interactions that do not have 0 as value
    l = []

    for j in range(1, interact_matrix[1].shape[0]):
        for i in range(j):
            if interact_matrix[1][i][j] != 0:
                s = str(i*l_res) + " " + str((i+1)*l_res-1) + " " + str(j*l_res) + " " + str((j+1)*l_res-1) + " " + str(interact_matrix[1][i][j]) + "\n"
                l.append(s)
    
    m = ((interact_matrix[0]).shape[0]+1)*[str(2*rm)+"\n"]

    n = np.asarray(interact_matrix[0])



    with open(p, "w+") as f:
        f.writelines(l)

    with open(p2, "w+") as g:
        g.writelines(m)

    np.savetxt(p3, n, delimiter=" ") #HR matrix / RCMC matrix
    print(len(m), len(l), n.shape[0])
