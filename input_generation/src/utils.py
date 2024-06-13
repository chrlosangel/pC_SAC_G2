import numpy as np
import scipy as scp
import pandas as pd
import hicstraw
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.pyplot as plt 
import seaborn as sns
import sys
from math import trunc
import os

from scipy.sparse import coo_matrix
from scipy.stats import pearsonr
import subprocess

np.set_printoptions(threshold=sys.maxsize)

clamp_n = True
diag_n = False
norm_q = False
SAVE = True

def run_pcSAC(downsampling_factor,lrf,run_option,number_chains,region_size,hic_res,split_chains,cell,chromosome):
    cmd = [
        "sbatch", "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/02.restrictions.sh",
        "-a", str(chromosome),
        "-d", str(downsampling_factor),
        "-r", str(region_size),
        "-h", str(hic_res),
        "-l", str(lrf),
        "-o", str(run_option),
        "-n", str(number_chains),
        "-s", str(split_chains),
        "-c", str(cell)
    ]
    subprocess.run(cmd)

def run_reconstruction(downsampling_factor,lrf,run_option,number_chains,region_size,hic_res,chromosome,start):
    #Each submatrix is ran with that interaction_distance, therefore main matrix should be ran with the same interaction_distance
    collision=int(hic_res) / 100 #bp to amstrongs
    interaction_distance=int((collision * 2) + 5)
    cmd = [
        "sbatch", "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/03.compare_matrices.sh",
        "-d", str(downsampling_factor),
        "-z", str(region_size),
        "-h", str(hic_res),
        "-l", str(lrf),
        "-o", str(run_option),
        "-n", str(number_chains),
        "-c", str(chromosome),
        "-s", str(start),
        "-i", str(interaction_distance)
    ]
    subprocess.run(cmd)


def run_genData(chr,start,end,downsample_factor,hic_res,lrbin,split_chains,file):
    default_file="/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/data/hic_datasets/GM12878/GM12878-HRC.hic"
    bool_res=0
    cmd = [
        "sbatch", "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/01.generate_data.sh",
        str(file),
        str(chr),
        str(start),
        str(end),
        str(downsample_factor),
        str(hic_res),
        str(lrbin),
        str(bool_res),
        str(split_chains)
    ]
    subprocess.run(cmd)

def interaction_file(interact_matrix,hic_file,hic_res,size,d_factor,l_res,chr="chr20"):
    filename = os.path.splitext(os.path.basename(hic_file))[0]
    pwd = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning"
    LR_F = int((hic_res * l_res) / 1000)
    if(clamp_n != True):
        output_dir = os.path.join(pwd, 
                    "results", 
                    "input_matrices",
                    chr, 
                    "{}_data_{}_resolution_{}_Mb_{}_downSample_{}_kb_LR_{}".format(
                        filename, 
                        hic_res, 
                        size, 
                        d_factor, 
                        LR_F,
                        chr ))
    if(clamp_n == True):
        norm_type = "clamp"
        output_dir = os.path.join(pwd, 
                    "results",norm_type, 
                    "input_matrices", 
                    chr,
                    "{}_data_{}_resolution_{}_Mb_{}_downSample_{}_kb_LR_{}".format(
                        filename, 
                        hic_res, 
                        size, 
                        d_factor, 
                        LR_F,
                        chr ))
        
    os.makedirs(output_dir, exist_ok=True)

    rm = (int(hic_res)/100)/2 #radius in amstrongs

    pn = "interaction_matrix_"+filename+chr+"_.txt"
    p = os.path.join(output_dir,pn)
    p2n = "int_mat_seg_len_"+filename+chr+"_.txt"
    p2 = os.path.join(output_dir,p2n)
    p3n = "whole_matrix_"+filename+chr+"_.txt" #whole matrix is normalized
    p3 = os.path.join(output_dir,p3n)
    p4n = "downsample_matrix_"+filename+chr+"_.txt" #P_N
    p4 = os.path.join(output_dir,p4n)
    p5n = "lr_downsample_matrix_"+filename+chr+"_.txt" #P_N
    p5 = os.path.join(output_dir,p5n)

    Path(p).touch(exist_ok=True)

    #l = [str(i*l_res)+" "+str((i+1)*l_res-1)+" "+str(j*l_res)+" "+str((j+1)*l_res-1)+" "+str((interact_matrix[4])[i][j])+"\n" for j in range(1,(interact_matrix[4]).shape[0]) for i in range(j)]
    #Giving only interactions that do not have 0 as value
    l = []

    for j in range(1, interact_matrix[4].shape[0]):
        for i in range(j):
            if interact_matrix[4][i][j] != 0:
                s = str(i*l_res) + " " + str((i+1)*l_res-1) + " " + str(j*l_res) + " " + str((j+1)*l_res-1) + " " + str(interact_matrix[4][i][j]) + "\n"
                l.append(s)

    m = ((interact_matrix[2]).shape[0]-1)*[str(2*rm)+"\n"]
    n = np.asarray(interact_matrix[2])
    p_arr = np.asarray(interact_matrix[3])
    lr = np.asarray(interact_matrix[4])

    with open(p, "w+") as f:
        f.writelines(l)
    with open(p2, "w+") as g:
        g.writelines(m)
    np.savetxt(p3, n, delimiter=" ")
    np.savetxt(p4, p_arr, delimiter=" ")
    np.savetxt(p5, lr, delimiter=" ")
    print(len(m), len(l), n.shape[0])

def interaction_file_split(interact_matrix,hic_file,hic_res,size,d_factor,l_res):
    filename = os.path.splitext(os.path.basename(hic_file))[0]
    pwd = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning"
    LR_F = int((hic_res * l_res) / 1000)
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
                    "results",norm_type, 
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

    print("Creating inversed input matrix...")
    pninv = "inv_interaction_matrix_"+filename+".txt"
    pinv = os.path.join(output_dir,pninv)
    print(pinv)

    p2n = "int_mat_seg_len_"+filename+".txt"
    p2 = os.path.join(output_dir,p2n)
    p3n = "whole_matrix_"+filename+".txt" #whole matrix is normalized
    p3 = os.path.join(output_dir,p3n)
    p4n = "downsample_matrix_"+filename+".txt" #P_N
    p4 = os.path.join(output_dir,p4n)
    p5n = "lr_downsample_matrix_"+filename+".txt" #P_N
    p5 = os.path.join(output_dir,p5n)

    Path(p).touch(exist_ok=True)

    #l = [str(i*l_res)+" "+str((i+1)*l_res-1)+" "+str(j*l_res)+" "+str((j+1)*l_res-1)+" "+str((interact_matrix[4])[i][j])+"\n" for j in range(1,(interact_matrix[4]).shape[0]) for i in range(j)]
    #Giving only interactions that do not have 0 as value
    l = []

    for j in range(1, interact_matrix[4].shape[0]):
        for i in range(j):
            if interact_matrix[4][i][j] != 0:
                s = str(i*l_res) + " " + str((i+1)*l_res-1) + " " + str(j*l_res) + " " + str((j+1)*l_res-1) + " " + str(interact_matrix[4][i][j]) + "\n"
                l.append(s)
    
    M_I = np.rot90(interact_matrix[4])
    M_I_2 = np.rot90(M_I) # completely inversed matrix

    l_inv = []

    for j in range(1, M_I_2.shape[0]):
        for i in range(j):
            if  M_I_2[i][j] != 0:
                s = str(i*l_res) + " " + str((i+1)*l_res-1) + " " + str(j*l_res) + " " + str((j+1)*l_res-1) + " " + str( M_I_2[i][j]) + "\n"
                l_inv.append(s)

    m = ((interact_matrix[2]).shape[0]-1)*[str(2*rm)+"\n"]
    n = np.asarray(interact_matrix[2])
    p_arr = np.asarray(interact_matrix[3])
    lr = np.asarray(interact_matrix[4])

    with open(p, "w+") as f:
        f.writelines(l)

    with open(pinv, "w+") as inv:
        inv.writelines(l_inv)

    with open(p2, "w+") as g:#int_mat_seg_len, represents the length of each segment in the HR matrix
        g.writelines(m)

    np.savetxt(p3, n, delimiter=" ") #whole matrix
    np.savetxt(p4, p_arr, delimiter=" ") 
    np.savetxt(p5, lr, delimiter=" ")
    print(len(m), len(l), n.shape[0])


def plot_hic_map(dense_matrix, maxcolor, suffix,hic_file,hic_res,size,d_factor,l_res,chr="chr20"):
    LR_F = int((hic_res * l_res) / 1000)
    plt.matshow(dense_matrix, cmap="rocket_r", vmin=0, vmax=maxcolor)
    cbar = plt.colorbar(shrink=0.8)
    cbar.ax.set_ylabel(ylabel='',rotation=270, labelpad=11)
    filename = os.path.splitext(os.path.basename(hic_file))[0]
    pwd = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning"

    #output_dir = pwd + "/" + "results" + "/" + "plots" + "/" + filename + "_data_" + str(hic_res) +"_resolution_"+ str(size) + "_Mb_" + str(d_factor) + "downSample_" + str(LR_F) + "_kb_LR"
    if(clamp_n != True):
        output_dir = os.path.join(pwd, 
                    "results", 
                    "plots",
                    chr, 
                    "{}_data_{}_resolution_{}_Mb_{}_downSample_{}_kb_LR_{}".format(
                        filename, 
                        hic_res, 
                        size, 
                        d_factor, 
                        LR_F,
                         chr ))
        
    if(clamp_n == True):
        norm_type = "clamp"
        output_dir = os.path.join(pwd, 
                    "results",norm_type, 
                    "plots", 
                    chr,
                    "{}_data_{}_resolution_{}_Mb_{}_downSample_{}_kb_LR_{}".format(
                        filename, 
                        hic_res, 
                        size, 
                        d_factor, 
                        LR_F,
                         chr ))
        
    os.makedirs(output_dir, exist_ok=True)
    pln = suffix + '_HiCMatrix.png'
    plp = os.path.join(output_dir, pln)
    plt.savefig(plp)
      

def dense2tag(matrix):
    """converting a square matrix (dense) to coo-based tag matrix"""
    matrix = np.round(matrix).astype(int)
    matrix = np.triu(matrix)
    tag_len = int(np.sum(matrix)) #all non zeros in matrix
    tag_mat = np.zeros((tag_len, 2), dtype=int)# 2D array of zeros with dimensions (tag_len, 2)
    coo_mat = coo_matrix(matrix)# represents the matrix in a compressed format 
    row, col, data = coo_mat.row, coo_mat.col, coo_mat.data #These attributes represent the row indices, column indices, and data values of the non-zero elements in the COO format, respectively. Only data 1
    start_idx = 0 
    for i in range(len(row)):
        end_idx = start_idx + data[i] #getting coordinates
        tag_mat[start_idx:end_idx, :] = (row[i], col[i]) #coordinates of non zero in original matrix
        start_idx = end_idx
    return tag_mat, tag_len

def tag2dense(tag, nsize):
    """coverting a coo-based tag matrix to densed square matrix."""
    coo_data, data = np.unique(tag, axis=0, return_counts=True) #data number of times it was selected
    row, col = coo_data[:, 0], coo_data[:, 1] #putting coordinates back into a matrix
    dense_mat = coo_matrix((data, (row, col)), shape=(nsize, nsize)).toarray()
    dense_mat = dense_mat + np.triu(dense_mat, k=1).T
    return dense_mat

def downsampling(matrix, down_ratio, verbose=False):
    """downsampling method"""
    if verbose: print(f"[Downsampling] Matrix shape is {matrix.shape}")
    tag_mat, tag_len = dense2tag(matrix)
    sample_idx = np.random.choice(tag_len, tag_len//down_ratio)#selecting tag_len//down_ratio from tag_len 
    sample_tag = tag_mat[sample_idx]# getting coordinates from the coo matrix of the selected items
    if verbose: print(f'[Downsampling] Sampling 1/{down_ratio} of {tag_len} reads')
    down_mat = tag2dense(sample_tag, matrix.shape[0])
    return down_mat

#def interactions_map_LR_5(HRmat, lres=4):
    s = HRmat.shape[0]
    if(s%lres>0):
        print("Impossible with LR5 method, trying with LR4.")
        return interactions_map_LR_4(HRmat)
    nbin = s//lres
    return HRmat.reshape(nbin,lres,nbin,lres).mean(3).mean(1)

def interactions_map_LR(HRmat, lres):
    nbin = int(HRmat.shape[0]/lres)
    M = np.triu([[((i)*(2*nbin-i+1)/2)+j-i+1 for j in range(nbin)] for i in range(nbin)])
    nlres=HRmat.shape[0]/nbin
    M_bin = np.array([[M[int(i/nlres), int(j/nlres)] for j in range(HRmat.shape[0])] for i in range(HRmat.shape[0])])
    
    LLR = [HRmat[M_bin == i].mean() for i in range(1, 1+int(nbin*(nbin+1)/2))]
    MLR = np.zeros((nbin,nbin)) 
    inds = np.triu_indices(nbin) 
    MLR[inds] = LLR
    MLR[(inds[1], inds[0])] = LLR
    return MLR

def quantile_norm(matrix):
    matrix_norm = np.copy(matrix).astype(float)
    upper_quartile = np.ceil(np.percentile(matrix_norm, 90)).astype(float)
    matrix_norm /= upper_quartile
    matrix_norm = np.where(matrix_norm > 1, 1, matrix_norm)
    return matrix_norm

def diag_norm(matrix):
    diagonal_elements = np.diag(matrix, k=0).astype(float)
    elements_after_0 = diagonal_elements[diagonal_elements > 0]
    res = np.min(elements_after_0)
    matrix_norm = np.copy(matrix).astype(float)
    matrix_norm = np.where(matrix_norm >= res, res, matrix_norm)
    matrix_norm /= res
    return matrix_norm

def clamping_norm(matrix):
    cutoff = np.ceil(np.percentile(matrix, 99)).astype(float)
    print(cutoff)
    matrix_norm = np.minimum(cutoff, matrix) #Replacing the values above the cutoff with the cutoff value
    matrix_norm = np.maximum(matrix_norm, 0) #Replacing the values below 0 with 0

    #Rescaling
    matrix_norm = matrix_norm / np.max(matrix_norm) #Dividing all values by the maximum value
    return matrix_norm

#At the moment only works with same chromosome and same coordinates plots
def interaction_maps_read_and_average(hic_file,hic_res,chromosome,start,end,bool_res,d_factor,size,lres):
    hic = hicstraw.HiCFile(hic_file)
    matrix_object = hic.getMatrixZoomData(chromosome, chromosome , "observed", "KR", "BP", hic_res)
    numpy_matrix = matrix_object.getRecordsAsMatrix(start, end, start, end)
    M = numpy_matrix
    if(bool_res == 0):
        P = downsampling(M,d_factor)
    else:
        P = M 
    if(norm_q):
        M_N = quantile_norm(M)
        P_N = quantile_norm(P)
    if(diag_n):
        M_N = diag_norm(M)
        P_N = diag_norm(P)
    if(clamp_n):
        M_N = clamping_norm(M)
        P_N = clamping_norm(P)

    P_LR = interactions_map_LR(P_N,lres)
    return M,P,M_N,P_N,P_LR

def gen_mat(hic_file,hic_res,chromosome,start,end,d_factor,bool_res,size,lres):
    hic_file=hic_file
    hic_res=hic_res
    chromosome=chromosome
    start=start
    end=end
    d_factor=d_factor
    bool_res=bool_res
    size=size
    lres=lres
    M,P,M_N,P_N,P_LR = interaction_maps_read_and_average(hic_file,hic_res,chromosome,start,end,bool_res,d_factor,size,lres)
    if(SAVE):
        hic_file=hic_file
        hic_res=hic_res
        d_factor=d_factor
        size=size
        lres=lres
        chr="chr"+str(chromosome)
        interaction_file((M,P,M_N,P_N,P_LR),hic_file,hic_res,size,d_factor,lres,chr)
        plot_hic_map(M_N,1,"M_N",hic_file,hic_res,size,d_factor,lres,chr) # plotting normalized Main matrix
        plot_hic_map(P_N,1,"P_N",hic_file,hic_res,size,d_factor,lres,chr) # plotting normalized subsampled matrix
        plot_hic_map(P_LR,1,"P_LR",hic_file,hic_res,size,d_factor,lres,chr) #plotting normalized subsampled LR matrix


def gen_mat_split_chains(hic_file,hic_res,chromosome,start,end,d_factor,bool_res,size,lres):
    hic_file=hic_file
    hic_res=hic_res
    chromosome=chromosome
    start=start
    end=end
    d_factor=d_factor
    bool_res=bool_res
    size=size
    lres=lres
    M,P,M_N,P_N,P_LR = interaction_maps_read_and_average(hic_file,hic_res,chromosome,start,end,bool_res,d_factor,size,lres)
    if(SAVE):
        print("Splitting Chains for file {}...".format(hic_file))
        hic_file=hic_file
        hic_res=hic_res
        d_factor=d_factor
        size=size
        lres=lres
        interaction_file_split((M,P,M_N,P_N,P_LR),hic_file,hic_res,size,d_factor,lres)
        plot_hic_map(M_N,1,"M_N",hic_file,hic_res,size,d_factor,lres) # plotting normalized Main matrix
        plot_hic_map(P_N,1,"P_N",hic_file,hic_res,size,d_factor,lres) # plotting normalized subsampled matrix
        plot_hic_map(P_LR,1,"P_LR",hic_file,hic_res,size,d_factor,lres) #plotting normalized subsampled LR matrix
        M_In = np.rot90(P_LR)
        M_In_2 = np.rot90(M_In)
        plot_hic_map(M_In_2,1,"P_LR_Inv",hic_file,hic_res,size,d_factor,lres)
