# This source code is licensed under the license found in the
# LICENSE file in the root directory of this source tree.
# --------------------------------------------------------
# a script to help convert .npz files to .cool files.
# --------------------------------------------------------

import argparse
import os
import numpy as np
from tqdm import tqdm
from Arg_Parser import *

import cooler
import pandas as pd

def read_npz(
        file_name, 
        hic_caption='hic', 
        bound = None,
        multiple = 1,
        include_additional_channels = True,
        compact_idx_caption = 'compact',
        norm_caption = 'norm',
        ):

    data = np.load(file_name)
    if include_additional_channels:
        hic_matrix = data[hic_caption]
        hic_matrix[0] *= multiple
        
    else:
        hic_matrix = data[hic_caption]
        if len(hic_matrix.shape) >= 3:
            hic_matrix = hic_matrix[0]
        
        hic_matrix *= multiple
    
    if bound is not None:
        mask = np.zeros((hic_matrix.shape[-2], hic_matrix.shape[-1]))
        for i in range(mask.shape[0]):
            for j in range(max(i-bound+1, 0), min(i+bound, mask.shape[1])):         
                mask[i][j] = 1
        
        hic_matrix = hic_matrix * mask

    if compact_idx_caption in data:
        compact_idx = data[compact_idx_caption]
    else:
        compact_idx = [i for i in range(hic_matrix.shape[-2])]

    if norm_caption in data:
        norm = data[norm_caption]
    else:
        norm = np.ones(hic_matrix.shape[-2])

    return hic_matrix, compact_idx, norm

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument('--data-dir', type=str, required = True)
    parser.add_argument('--hic-caption', type=str, default = 'hic')
    parser.add_argument('--external-norm-file', type=str, 
                        default = 'NONE')
    parser.add_argument('--resolution', type=str, default='10kb')

    parser.add_argument('--bound', type=int, default=200)
    parser.add_argument('--multiple', type=int, default=255)

    parser.add_argument('-c', '--cell-line', default='GM12878')
    parser.add_argument('-s', '--dataset', default='test' )
    args = parser.parse_args()

    # data_dir = args.data_dir
    hic_caption = args.hic_caption
    external_norm_file = args.external_norm_file
    res = args.resolution
    cell_line = args.cell_line
    dataset = args.dataset
    bound = args.bound
    multiple = args.multiple
    resolution = res_map[res.split('_')[0]]
    data_dir = os.path.join(root_dir, 'predict', cell_line)

    chr_list = set_dict[dataset]
    abandon_chromosome = abandon_chromosome_dict.get(cell_line, [])

    bins = {"chrom": [], "start": [], "end": []}
    reads = {"bin1_id": [], "bin2_id": [], "count": []}
    for n in tqdm(chr_list):
        if n in abandon_chromosome:
            continue

        in_file = os.path.join(data_dir, f'predict_chr{n}_{res}.npz')

        matrix, compact_idx, norm = read_npz(in_file, hic_caption = hic_caption, bound = bound, multiple=multiple, include_additional_channels=False)

        if external_norm_file != 'NONE':
            if '#(CHR)' in external_norm_file:
                norm = read_singlechromosome_norm(external_norm_file, n, cell_line)
            else:
                raise NotImplementedError
        
        s = len(bins['chrom'])

        for i in range(matrix.shape[0]):
            bins['chrom'].append(f'chr{n}')
            bins['start'].append(i*resolution)
            bins['end'].append((i+1)*resolution)

        for i in tqdm(range(matrix.shape[0]), leave=False):
            for j in range(i, min(i+bound, matrix.shape[1])):
                RAW = float(matrix[i][j])*norm[i]*norm[j]
                if RAW > 0:
                    reads['bin1_id'].append(s+i)
                    reads["bin2_id"].append(s+j)
                    reads["count"].append(matrix[i][j])

    cool_file = os.path.join(data_dir, f"bound{bound}_{res}.cool")
    cooler.create_cooler(cool_file, bins=pd.DataFrame.from_dict(bins), pixels=pd.DataFrame.from_dict(reads), dtypes={'count': float})

