#pcSAC modules
from . import reconstruct 
from . import utils


import numpy as np

from scipy.stats import pearsonr
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics.pairwise import cosine_similarity

def create_mask(matrix, n_off_diagonals=3):
    mask = np.ones_like(matrix, dtype=bool)
    for k in range(-n_off_diagonals, n_off_diagonals + 1):
        np.fill_diagonal(mask, 0, wrap=True)
        mask = np.logical_and(mask, np.eye(*matrix.shape, k=k, dtype=bool) == 0)
        upper = np.triu(matrix) * mask
    return upper

def get_metrics(upper_reconstructed, upper_whole_M):
    metrics = {
        'SSIM': reconstruct.ssim(upper_reconstructed, upper_whole_M),
        'PCC': reconstruct.pearsons(upper_reconstructed, upper_whole_M),
        'SCC': reconstruct.spearmans(upper_reconstructed, upper_whole_M),
        'MSE': reconstruct.mse(upper_reconstructed, upper_whole_M),
        'MAE': reconstruct.mae(upper_reconstructed, upper_whole_M),
        'PSNR': reconstruct.psnr(upper_reconstructed, upper_whole_M),
        'COSINE': reconstruct.cosine(upper_reconstructed, upper_whole_M)
    }

    return metrics

def mse_gd(dict1, dict2):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    result = {}
    
    for key in shared_keys:
        vector1 = dict1[key]
        vector2 = dict2[key]
        r = mean_squared_error(vector1, vector2)
        result[key] = r

    return result

def mae_gd(dict1, dict2):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    result = {}
    
    for key in shared_keys:
        vector1 = dict1[key]
        vector2 = dict2[key]
        r = mean_absolute_error(vector1, vector2)
        result[key] = r

    return result

def pearson_gd(dict1, dict2):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    result = {}
    
    for key in shared_keys:
        vector1 = dict1[key]
        vector2 = dict2[key]
        if np.std(vector1) == 0 and np.std(vector2) == 0:
            result[key] = 1
        else:
            r,_ = pearsonr(vector1, vector2)
            result[key] = r

    return result
        

def spearman_gd(dict1, dict2):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    result = {}
    
    for key in shared_keys:
        vector1 = dict1[key]
        vector2 = dict2[key]
        if np.std(vector1) == 0 and np.std(vector2) == 0:
            result[key] = 1
        else:
            r,_ = spearmanr(vector1, vector2)
            result[key] = r

    return result

def cosine_gd(dict1, dict2):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    result = {}
    
    for key in shared_keys:
        vector1 = dict1[key]
        vector2 = dict2[key]
        r = cosine_similarity(vector1, vector2)
        result[key] = r 
        
    return result

def GenomicDistance_correlation(reconstructed_matrix, whole_M,resolution,genomic_metric):
    df_rec = reconstruct.matrix2df(reconstructed_matrix)
    df_whole = reconstruct.matrix2df(whole_M)

    gdist_rec = reconstruct.genomicDistance(df_rec, resolution)
    gdist_org = reconstruct.genomicDistance(df_whole, resolution)


    cs_gd = genomic_metric(gdist_rec, gdist_org)
    min_val = resolution * 2
    max_val = resolution * 150
    cs_gd_sel = {k: v for k, v in cs_gd.items() if min_val <= k <= max_val}

    return cs_gd_sel


def calculate_mean(metric_data):
    from collections import defaultdict
    sums = defaultdict(float)
    counts = defaultdict(int)

    for metric_dict in metric_data:
        for key, value in metric_dict.items():
            sums[key] += value
            counts[key] += 1

    return {key: sums[key] / counts[key] for key in sums}
