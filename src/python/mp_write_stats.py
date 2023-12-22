#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 13:34:00 2023

@author: valentin
"""

import multiprocessing
import argparse
from functools import partial
import glob
import numpy as np
from collections import Counter
from Bio import AlignIO
import tqdm
import time
import os

def to_d3(mat):
    ## mat must be(nb_conf, nb_res*3)
    ## output is (nb_conf, nb_res, 3)
    if isinstance(mat, np.ndarray):
        return np.reshape(mat, (len(mat), len(mat.T) // 3, 3))
    else:
        raise ValueError("Unrecognized matrix type")

def to_d2(mat):
    ## mat must be (nb_conf, nb_res, 3)
    ## output is (nb_conf, nb_res*3)
    if isinstance(mat, np.ndarray):
        shape = np.shape(mat)
        return np.reshape(mat, (shape[0], shape[1] * 3))
    else:
        raise ValueError("Unrecognized matrix type")

def eigss_svd(X):
    #input shape (nb_conf, nb_res*3)
    #Vt.T dimensions (nb_res*3, modes)
    _, S, Vt = np.linalg.svd(X.T, full_matrices=False)
    explained_variance_ = (S**2) / (len(X) - 1)
    total_var = explained_variance_.sum()
    explained_variance_ratio = explained_variance_ / total_var
    return explained_variance_ratio, Vt.T

def apply_query_coord_on_missing_data(coord_mat, K_mat, query_id):
    #input coord_mat (nb_res*3, nb_conf)
    add_mat = np.logical_not(K_mat).astype(float) * coord_mat[:, query_id, None]
    coord_mat = coord_mat * K_mat
    coord_mat += add_mat
    return coord_mat

def apply_mean_coord_on_missing_data(coord_mat,K_mat):
    coord_mat = coord_mat+(np.ones(np.shape(coord_mat))*(np.sum(coord_mat,axis=1)/np.sum(K_mat,axis=1))[:,np.newaxis]*(np.ones(np.shape(K_mat))-K_mat))
    return(coord_mat)

def get_query(coordinates, gaps, query_index, coverage = False, normalize=False):
    #Input shape coords : (nb_res*3, nb_conf)
    #Input shape gaps : (nb_res*3, nb_conf)
    #In use for the query centered part
    # Extract the coordinates and gaps for the query index
    query_gaps = gaps[:, query_index:query_index + 1]
    indices = np.where(query_gaps[:, 0] == 1)[0]
    # Filter coordinates and gaps with these indices
    filtered_coordinates = coordinates[indices, :]
    filtered_gaps = gaps[indices, :]

    # Apply the query coordinates on missing data
    coordinates = apply_query_coord_on_missing_data(filtered_coordinates, filtered_gaps, query_index)

    if normalize:
        # Get the query coordinates
        query_coords = coordinates[:, query_index:query_index + 1]
        # Compute the deviations from the query
        deviations = coordinates - query_coords
        # Compute the squared deviations from the query coordinates
        squared_deviations = deviations ** 2
        # Manually compute the mean squared deviation with ddof=1
        # Subtract 1 from the number of observations to adjust for sample standard deviation
        mean_squared_deviation = np.sum(squared_deviations, axis=1, keepdims=True) / (squared_deviations.shape[1] - 1)
        # Compute the standard deviation from the query coordinates
        std_from_query = np.sqrt(mean_squared_deviation)
        # Avoid division by zero in case of zero standard deviation
        std_from_query[std_from_query == 0] = 1
        # Standardize the coordinates
        coordinates = deviations / std_from_query
    else:
        coordinates = coordinates - coordinates[:, query_index:query_index + 1]

    if coverage:
        # Filter gaps based on the query subset
        filtered_gaps_3d = to_d3(filtered_gaps.T)
        coverage_subset = np.sum(filtered_gaps_3d[..., 0], axis=0) / filtered_gaps_3d.shape[0]
        coordinates_3d = to_d3(coordinates.T)
        coordinates_3d = coordinates_3d*coverage_subset[None,:,None]
        coordinates = to_d2(coordinates_3d).T

    return coordinates

def get_mean(coordinates, gaps, coverage=False, normalize=False):
    # Input shape coords: (nb_res*3, nb_conf)
    # Input shape gaps: (nb_res*3, nb_conf)
    # Apply the query coordinates on missing data
    coordinates = apply_mean_coord_on_missing_data(coordinates, gaps)
    

    if normalize:
        mean = np.mean(coordinates, axis=1)[:, None]
        # Compute the deviations from the query
        deviations = coordinates - mean
        # Compute the squared deviations from the mean coordinates
        squared_deviations = deviations ** 2
        mean_squared_deviation = np.sum(squared_deviations, axis=1, keepdims=True) / (squared_deviations.shape[1] - 1)
        std_from_mean = np.sqrt(mean_squared_deviation)
        std_from_mean[std_from_mean == 0] = 1
        # Standardize the coordinates
        coordinates = deviations / std_from_mean
    
    # If not normalizing, just adjust the coordinates by subtracting the mean
    else:
        coordinates = coordinates - np.mean(coordinates, axis=1)[:, None]
    
    # Calculate coverage if required
    if coverage:
        gaps_3d = to_d3(gaps.T)
        coverage = np.sum(gaps_3d[...,0], axis=0) / len(gaps_3d)
        coordinates_3d = to_d3(coordinates.T)
        coordinates_3d = coordinates_3d * coverage[None, :, None]
        coordinates = to_d2(coordinates_3d).T
    
    return coordinates

def load_tensor(filename):
    with open(filename, 'rb') as f:
        numModels = np.frombuffer(f.read(8), dtype=np.int64)[0]
        numSeqs = np.frombuffer(f.read(8), dtype=np.int64)[0]
        numCoords = np.frombuffer(f.read(8), dtype=np.int64)[0]
        data_shape = (numModels, numSeqs * numCoords)
        data = np.frombuffer(f.read(), dtype=np.float64)
        return data.reshape(data_shape)
    
def load_mask(filename):
    with open(filename, 'rb') as f:
        numModels = np.frombuffer(f.read(8), dtype=np.int64)[0]
        numSeqs = np.frombuffer(f.read(8), dtype=np.int64)[0]
        tensor = np.frombuffer(f.read(numModels * numSeqs), dtype=np.uint8).astype(bool).reshape((numModels, numSeqs))
        tensor_3d = np.repeat(tensor[:, :, np.newaxis], 3, axis=2)
        tensor = np.reshape(tensor_3d, (numModels, numSeqs*3))
        return tensor
   
def getCoverage(aln,minCov=0.8):
    cov = 1-np.count_nonzero((aln=='-')|(aln=='X'),axis=0)/len(aln)
    return(np.count_nonzero(cov>minCov)/len(aln.T))

def getMaxScores(aln):
    nb_aas = max(np.count_nonzero((aln!='-') & (aln!='X'),axis=1))
    nb_seqs = len(aln)
    return(nb_seqs*(nb_seqs-1)*nb_aas/2)

def combi2(x):
    return(x*(x-1)/2)

# sum-of-pairs score computed for a MSA, with fixed match, mismatch and gap values
# the score is normalized by the score expected for an ideal alignment 
# (i.e. a MSA with only matches and of the lenght equal to the max number of aa, computed over the sequences in the MSA)
# usage cases:
# 1) % of identity can be retrieved by setting the mismatch penalty to zero
# and disregarding the gaps and Xs:  computeSumOfPairs(readAli("6CCHB"), mismatch=0,  gapChar=c("-","x"))
# 2) global quality of the MSA can be assessed by setting the gap penalty equal to the mismatch 
# 3) quality of the gaps can be assessed and compared by setting the gap penalty to some value

def computeSumOfPairs(aln,match=1,mismatch=-0.5,gap=-0.5):
    nb_seqs = len(aln)
    nb_pos = len(aln.T)
    Ctot = nb_seqs*(nb_seqs-1)/2
    score = 0
    for j in range(nb_pos):
        column = aln[:,j]
        table = Counter(column)
        tgap = 0
        try:
            tgap+=table.pop('-')
        except:
            pass
        try:
            tgap+=table.pop('X')
        except:
            pass
        if len(table)>0:
            Cmatch = np.sum(combi2(np.array([i for i in table.values()])))
        else: 
            Cmatch = 0
        if tgap>0:    
            Cgaps = combi2(tgap) + tgap*(nb_seqs-tgap)
        else:
            Cgaps = 0
        score = score + Cmatch * match + Cgaps * gap + (Ctot - Cgaps - Cmatch) * mismatch 
    return(score/getMaxScores(aln))

def compute_col(comp):
    comp +=1e-20
    compr = np.reshape(comp,(len(comp),int(len(comp.T)/3),3))
    compr = np.sum(compr**2,axis=2)
    col = 1/len(compr.T)*np.exp(-np.sum(compr*np.log(compr),axis=1))
    return(col)
   
def threshold_value(cumsum, threshold):
    return str(np.argmax(cumsum >= threshold) + 1)

def compute_stats(idclu, directory, use_weights=False):
    try:
        thresholds = [0.5, 0.8, 0.85, 0.9, 0.95, 0.99]

        ### msa operations
        al = AlignIO.read(directory + idclu+'_aln.fa', 'fasta')
        aln = np.array(al)
        cov = f'{getCoverage(aln):.3f}'
        globalq = f'{computeSumOfPairs(aln):.3f}'
        pid = f'{computeSumOfPairs(aln,mismatch=0,gap=0):.3f}'
        num_mdl, aln_len = len(aln), len(aln.T)
        mask_res = (aln[0] != '-') & (aln[0] != 'X')
        ref_len = len(aln[0][mask_res])
        ref_name = al[0].id.split('_')[0]
        ref_aln = aln[:, mask_res]
        ref_cov = f'{getCoverage(ref_aln):.3f}'
        ref_globalq = f'{computeSumOfPairs(ref_aln):.3f}'
        ref_pid = f'{computeSumOfPairs(ref_aln,mismatch=0,gap=0):.3f}'

        ### rmsd operations
        mat = np.loadtxt(directory + idclu+'_rmsd.txt')
        ref_rmsd_mat = mat[0][1:]
        ref_rmsd_max, ref_rmsd_mean, ref_rmsd_std = f'{np.nanmax(ref_rmsd_mat):.3f}', f'{np.nanmean(ref_rmsd_mat):.3f}', f'{np.nanstd(ref_rmsd_mat):.3f}'
        mat= mat[np.triu_indices(len(mat),k=1)]
        rmsd_max, rmsd_mean, rmsd_std = f'{np.nanmax(mat):.3f}', f'{np.nanmean(mat):.3f}', f'{np.nanstd(mat):.3f}'

        ### coordinates loading
        coords = load_tensor(directory + idclu+'_raw_coords_ca.bin').T
        mask = load_mask(directory + idclu+'_raw_coords_ca_mask.bin').T
        col_to_keep = np.where(np.sum(mask,axis=1) > 1)[0] #toutes les colonnes où il y a plus de 2 résidus
        mask = mask[col_to_keep]
        coords=coords[col_to_keep]
        

        if use_weights:
            
            eigval_weights_query, eigvec_weights_query = eigss_svd(get_query(coords,mask,0,coverage=True))
            eigval_weights_mean, eigvec_weights_mean = eigss_svd(get_mean(coords,mask,coverage=True))
            cumsum_weights_query = np.cumsum(eigval_weights_query)
            cumsum_weights_mean = np.cumsum(eigval_weights_mean)
            nb_modes_weights_query = np.argmax(cumsum_weights_query>=0.95)+1
            nb_modes_weights_mean = np.argmax(cumsum_weights_mean>=0.95)+1
            comp_pca_weights_query = eigvec_weights_query[:,0:nb_modes_weights_query]
            comp_pca_weights_mean = eigvec_weights_mean[:,0:nb_modes_weights_mean]
            col_pca_weights_query = compute_col(comp_pca_weights_query.T)
            col_pca_weights_mean = compute_col(comp_pca_weights_mean.T)
            
            eigval_weights_query_norm, eigvec_weights_query_norm = eigss_svd(get_query(coords, mask, 0, coverage=True, normalize=True))
            eigval_weights_mean_norm, eigvec_weights_mean_norm = eigss_svd(get_mean(coords, mask, coverage=True, normalize=True))
            cumsum_weights_query_norm = np.cumsum(eigval_weights_query_norm)
            cumsum_weights_mean_norm = np.cumsum(eigval_weights_mean_norm)
            nb_modes_weights_query_norm = np.argmax(cumsum_weights_query_norm >= 0.95) + 1
            nb_modes_weights_mean_norm = np.argmax(cumsum_weights_mean_norm >= 0.95) + 1
            comp_pca_weights_query_norm = eigvec_weights_query_norm[:, 0:nb_modes_weights_query_norm]
            comp_pca_weights_mean_norm = eigvec_weights_mean_norm[:, 0:nb_modes_weights_mean_norm]
            col_pca_weights_query_norm = compute_col(comp_pca_weights_query_norm.T)
            col_pca_weights_mean_norm = compute_col(comp_pca_weights_mean_norm.T)
            
            with open(directory + idclu+'_var_col_weights.csv','w') as f:
                f.write('%var_mean,col_pca_mean,%var_query,col_pca_query\n')
                for valm,colvalm,valq,colvalq in zip(eigval_weights_mean,col_pca_weights_mean,eigval_weights_query,col_pca_weights_query):
                    f.write(f'{valm:.5f}'+','+f'{colvalm:.3f}'+','+f'{valq:.5f}'+','+f'{colvalq:.3f}'+'\n')


            with open(directory + idclu+'_var_col_weights_norm.csv','w') as f:
                f.write('%var_mean,col_pca_mean,%var_query,col_pca_query\n')
                for valm,colvalm,valq,colvalq in zip(eigval_weights_mean_norm,col_pca_weights_mean_norm,eigval_weights_query_norm,col_pca_weights_query_norm):
                    f.write(f'{valm:.5f}'+','+f'{colvalm:.3f}'+','+f'{valq:.5f}'+','+f'{colvalq:.3f}'+'\n')

            line = (
            #centered on mean
            f"{idclu},{ref_name},{num_mdl},{aln_len},{ref_len},{cov},{pid},{globalq},"
            f"{rmsd_max},{rmsd_mean},{rmsd_std},"
            f"{','.join([threshold_value(cumsum_weights_mean, t) for t in thresholds])},"
            f"{eigval_weights_mean[0]:.3f},{col_pca_weights_mean[0]:.3f},"
            
            #centered on query
            f"{ref_cov},{ref_pid},{ref_globalq},{ref_rmsd_max},{ref_rmsd_mean},{ref_rmsd_std},"
            f"{','.join([threshold_value(cumsum_weights_query, t) for t in thresholds])},"
            f"{eigval_weights_query[0]:.3f},{col_pca_weights_query[0]:.3f}"

            #centered on mean normalized
            f",{','.join([threshold_value(cumsum_weights_mean_norm, t) for t in thresholds])},"
            f"{eigval_weights_mean_norm[0]:.3f},{col_pca_weights_mean_norm[0]:.3f},"

            #centered on query normalized
            f"{','.join([threshold_value(cumsum_weights_query_norm, t) for t in thresholds])},"
            f"{eigval_weights_query_norm[0]:.3f},{col_pca_weights_query_norm[0]:.3f}"

            )

        if not use_weights:
            eigval_query, eigvec_query = eigss_svd(get_query(coords, mask, 0))
            eigval_mean, eigvec_mean = eigss_svd(get_mean(coords, mask))
            cumsum_query = np.cumsum(eigval_query)
            cumsum_mean = np.cumsum(eigval_mean)
            nb_modes_query = np.argmax(cumsum_query >= 0.95) + 1
            nb_modes_mean = np.argmax(cumsum_mean >= 0.95) + 1
            comp_pca_query = eigvec_query[:, 0:nb_modes_query]
            comp_pca_mean = eigvec_mean[:, 0:nb_modes_mean]
            col_pca_query = compute_col(comp_pca_query.T)
            col_pca_mean = compute_col(comp_pca_mean.T)

            eigval_query_norm, eigvec_query_norm = eigss_svd(get_query(coords, mask, 0, normalize=True))
            eigval_mean_norm, eigvec_mean_norm = eigss_svd(get_mean(coords, mask, normalize=True))
            cumsum_query_norm = np.cumsum(eigval_query_norm)
            cumsum_mean_norm = np.cumsum(eigval_mean_norm)
            nb_modes_query_norm = np.argmax(cumsum_query_norm >= 0.95) + 1
            nb_modes_mean_norm = np.argmax(cumsum_mean_norm >= 0.95) + 1
            comp_pca_query_norm = eigvec_query_norm[:, 0:nb_modes_query_norm]
            comp_pca_mean_norm = eigvec_mean_norm[:, 0:nb_modes_mean_norm]
            col_pca_query_norm = compute_col(comp_pca_query_norm.T)
            col_pca_mean_norm = compute_col(comp_pca_mean_norm.T)
            with open(directory + idclu+'_var_col.csv', 'w') as f:
                f.write('%var_mean,col_pca_mean,%var_query,col_pca_query\n')
                for valm, colvalm, valq, colvalq in zip(eigval_mean, col_pca_mean, eigval_query, col_pca_query):
                    f.write(f'{valm:.5f}'+','+f'{colvalm:.3f}'+','+f'{valq:.5f}'+','+f'{colvalq:.3f}'+'\n')
            with open(directory + idclu+'_var_col_norm.csv', 'w') as f:
                f.write('%var_mean,col_pca_mean,%var_query,col_pca_query\n')
                for valm, colvalm, valq, colvalq in zip(eigval_mean_norm, col_pca_mean_norm, eigval_query_norm, col_pca_query_norm):
                    f.write(f'{valm:.5f}'+','+f'{colvalm:.3f}'+','+f'{valq:.5f}'+','+f'{colvalq:.3f}'+'\n')
            line = (
                #centered on mean
                f"{idclu},{ref_name},{num_mdl},{aln_len},{ref_len},{cov},{pid},{globalq},"
                f"{rmsd_max},{rmsd_mean},{rmsd_std},"
                f"{','.join([threshold_value(cumsum_mean, t) for t in thresholds])},"
                f"{eigval_mean[0]:.3f},{col_pca_mean[0]:.3f},"
                
                #centered on query
                f"{ref_cov},{ref_pid},{ref_globalq},{ref_rmsd_max},{ref_rmsd_mean},{ref_rmsd_std},"
                f"{','.join([threshold_value(cumsum_query, t) for t in thresholds])},"
                f"{eigval_query[0]:.3f},{col_pca_query[0]:.3f}"

                #centered on mean normalized
                f",{','.join([threshold_value(cumsum_mean_norm, t) for t in thresholds])},"
                f"{eigval_mean_norm[0]:.3f},{col_pca_mean_norm[0]:.3f},"

                #centered on query normalized
                f"{','.join([threshold_value(cumsum_query_norm, t) for t in thresholds])},"
                f"{eigval_query_norm[0]:.3f},{col_pca_query_norm[0]:.3f}"

                )
            
        return line
    
    except Exception as e:
        print(f"Error processing {idclu}")
        return f"ERROR_{idclu}"
    
def main(use_weights, directory):
    start_time = time.time()

    compute_stats_configured = partial(compute_stats, directory=directory, use_weights=use_weights)

    input_path = os.path.join(directory, '*_aln.fa')
    liste_id = glob.glob(input_path)
    liste_id = [os.path.basename(i).rsplit('_aln.fa', 1)[0] for i in liste_id]

    print('Computing statistics:')
    with multiprocessing.Pool() as pool:
        results = list(tqdm.tqdm(pool.imap_unordered(compute_stats_configured, liste_id), total=len(liste_id)))

    output_filename = 'stats_weighted.csv' if use_weights else 'stats.csv'
    output_path = os.path.join(directory, output_filename)

    with open(output_path, 'w') as f:
        header = 'name_file,name_ref,nb_members,aln_len,ref_len,coverage,percent_id,global_quality,rmsd_max,rmsd_mean,rmsd_std,50%,80%,85%,90%,95%,99%,%var_1st,col_1st,ref_coverage,ref_percent_id,ref_global_quality,ref_rmsd_max,ref_rmsd_mean,ref_rmsd_std,ref_50%,ref_80%,ref_85%,ref_90%,ref_95%,ref_99%,ref_%var_1st,ref_col_1st,50%_norm,80%_norm,85%_norm,90%_norm,95%_norm,99%_norm,%var_1st_norm,col_1st_norm,ref_50%_norm,ref_80%_norm,ref_85%_norm,ref_90%_norm,ref_95%_norm,ref_99%_norm,ref_%var_1st_norm,ref_col_1st_norm'
        f.write(header + '\n')
        for result_line in results:
            if not result_line.startswith("ERROR_"):
                f.write(result_line + '\n')
            else:
                error_idclu = result_line.split("_")[1]
                print(f"Error processing {error_idclu}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    #print(f"Finished in {elapsed_time:.2f} seconds")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute statistics.")
    parser.add_argument("--use_weights", action="store_true", help="Use this option to enable weights.")
    parser.add_argument("--directory", type=str, default=".", help="Directory for input and output files.")
    args = parser.parse_args()

    main(args.use_weights, args.directory)