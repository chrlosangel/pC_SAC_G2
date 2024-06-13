#Imports main modules to use 
import argparse
from src import utils,globals,submatrices

# Clamp option is inside utils.py

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description="HiC experimental data processing script")
parser.add_argument("--hic-file", required=True, help="Path of HiC file")
parser.add_argument("--hic-res", required=True, type=int, help="Resolution of HiC data")
parser.add_argument("--chromosome", required=True, help="Target chromosome")
parser.add_argument("--start", required=True, type=int, help="Start coordinate")
parser.add_argument("--end", required=True, type=int, help="End coordinate")
parser.add_argument("--d-factor", required=True, type=int, help="Downsampling factor")
parser.add_argument("--bool-res", required=True, type=int, help="Boolean resolution 0 or 1")
parser.add_argument("--l_res", required=True, type=int, help="LR map factor")
parser.add_argument("--submatrix", required=False, type=str2bool, help="Do you wish to generate submatrices? True or False")
parser.add_argument("--split_chains", required=False, type=str2bool, help="Do you wish to create files to split chains? True/yes or False/no")
# Parse arguments
args = parser.parse_args()

hic_file = args.hic_file
hic_res = args.hic_res
chromosome =  args.chromosome
start = args.start
end = args.end
d_factor = args.d_factor
bool_res = args.bool_res
size = end - start
lres = args.l_res
submatrix = args.submatrix
split_chains = args.split_chains

if(submatrix == True):
    print("Generating submatrices...")
    M,P,M_N,P_N,P_LR = utils.interaction_maps_read_and_average(hic_file,hic_res,chromosome,start,end,bool_res,d_factor,size,lres)
    coords = [(0,150,0,150),(130,310,130,310),(220,400,220,400)]
    submatrices.bin_matrix(M_N,"M_N",coords)
    submatrices.bin_matrix(P_N,"P_N",coords)
    submatrices.bin_matrix(P_LR,"P_LR",coords)
    
else:
    if(split_chains == True):
        print("Generating data, split chains = True...")
        utils.gen_mat_split_chains(hic_file,hic_res,chromosome,start,end,d_factor,bool_res,size,lres)
    else:      
        print("Generating data...")
        utils.gen_mat(hic_file,hic_res,chromosome,start,end,d_factor,bool_res,size,lres)

