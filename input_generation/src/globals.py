'''
    Cutoff percentage used for normalizing the hic matrices
'''
CUTOFF_PERCENTILE = 99.95


'''
    Pre-Computed cutoff values for datasets
'''

dataset_cutoff_values = {
    'GM12878_rao_et_al': 255,
    'GM12878_rao_et_al_replicate': 255,
    'IMR90_rao_et_al': 255,
    'K562_rao_et_al': 255,
    'downsampled_8': 140,
    'downsampled_16': 100,
    'downsampled_25': 80,
    'downsampled_50': 50,
    'downsampled_100': 25,
}
main_dir = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/results/clamp/"
hic_file = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/data/hic_datasets/GM12878/GM12878-HRC.hic"
hic_K562 = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/data/hic_datasets/K562/K562-HRC.hic"
hic_IMR90 = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/data/hic_datasets/IMR90/IMR90-HRC.hic"
chromosome = 20
hic_res = 10000
start = 100000
end = 2100000
size = int(end - start)
sizes="/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/expData_scripts/hg19_resources/hg19.chrom.sizes"
all="/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/results/clamp/plots/"
all_nc="/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/results/plots/"
input_matrices="/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/results/clamp/input_matrices/"
input_matrices_nc="/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/results/input_matrices/"
restrictions_path="/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/results/clamp/restrictions/"
restrictions_path_nc="/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/results/restrictions/"

'''Submatrices'''
submat_main = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/results/clamp/submatrices/"
submat_main_nc = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/results/submatrices/"
sub_reconstruction = "/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/results/clamp/submatrices/restrictions/"
coords =  ["0:180","130:310","220:400"]
reconstructed_matrix="/gpfs/commons/home/cangel/g2lab/projects/01_12_23_highResolutionHiC/scripts/ParameterTuning/results/clamp/restrictions/"
