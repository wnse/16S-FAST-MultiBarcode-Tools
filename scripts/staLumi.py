# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.10.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import logging
import pandas as pd
import os
from Bio.Seq import Seq
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
import numpy as np
from cutLadapter import cutLadapter


# %%
def rc_seq(seq):
    if seq == seq:
        return str(Seq(seq).reverse_complement())
    else:
        return None

def get_umi_info(f, col='left'):
    col_name = ['Read name',
            'Number of errors',
            'start',
            'end',
            'left',
            'matched',
            'right',
            'Name',
            'left_Q',
            'matched_Q',
            'right_Q',
           ]
    df_tmp = pd.read_csv(f, sep='\t',  names=col_name)
    df_tmp = df_tmp.set_index('Read name')
    df_tmp.index = df_tmp.index.str.split('\s').str[0]
    return df_tmp[['Number of errors', col, 'Name']].rename(columns={col:'umi'})

def get_df_umi(umi1_fq, umi1_fq_1, umi2_fq, umi2_fq_1, suffix = '.cutadapt.info.file'):    
    df_10 = get_umi_info(umi1_fq+suffix)
    df_11 = get_umi_info(umi1_fq_1+suffix, col='right')
    df_20 = get_umi_info(umi2_fq+suffix)
    df_21 = get_umi_info(umi2_fq_1+suffix, col='right')
    
    return df_10, df_11, df_20, df_21

def get_df_merge(umi1_fq, umi1_fq_1, umi2_fq, umi2_fq_1, umi_lmin=10, umi_lmax=20):
    df_10, df_11, df_20, df_21 = get_df_umi(umi1_fq, umi1_fq_1, umi2_fq, umi2_fq_1)
    df_10 = df_10.rename(columns={'Number of errors':'errors_10', 'umi':'umi_10', 'Name':'Name_10'})
    df_11 = df_11.rename(columns={'Number of errors':'errors_11', 'umi':'umi_11', 'Name':'Name_11'})
    df_20 = df_20.rename(columns={'Number of errors':'errors_20', 'umi':'umi_20','Name':'Name_20'})
    df_21 = df_21.rename(columns={'Number of errors':'errors_21', 'umi':'umi_21', 'Name':'Name_21'})
    df_merge = pd.concat([df_10, df_11, df_20, df_21], axis=1)
    total_reads = df_merge.shape[0]
    logging.info('L reads:\t{}'.format(total_reads))
    df_merge_f = df_merge[(df_merge[df_merge.columns[df_merge.columns.str.contains('errors')]]==0).all(axis=1)]
    logging.info('L reads after cut adapter/linkers:\t{}'.format(df_merge_f.shape[0]))

    df_merge_f = df_merge_f[(df_merge_f['umi_10'] == df_merge_f['umi_21'].apply(rc_seq)) & 
                            (df_merge_f['umi_11'].apply(rc_seq) == df_merge_f['umi_20'])
                           ]
    logging.info('L reads after fiter R1_F==R2_R and R1_R==R2_F:\t{}'.format(df_merge_f.shape[0]))
    df_merge_f = df_merge_f[(df_merge_f['Name_21'].str[0]=='F')&(df_merge_f['Name_11'].str[0]=='R')]
    logging.info('L reads after fiter FBarcode!="F" or RBarcode!="R":\t{}'.format(df_merge_f.shape[0]))

    df_merge_f['Barcode'] = df_merge_f['Name_21'] + '|' + df_merge_f['Name_11']
    df_merge_f = df_merge_f[['Barcode', 'umi_21', 'umi_11', 'Name_21', 'Name_11']]
    

    df_merge_f = df_merge_f[(df_merge_f['umi_11'].apply(len).between(umi_lmin, umi_lmax)) &
                            (df_merge_f['umi_21'].apply(len).between(umi_lmin, umi_lmax))
                           ]
    total_umi_reads = df_merge_f.shape[0]
    logging.info('L reads after filter length(umi) between {}-{}bp:\t{}'.format(umi_lmin, umi_lmax, total_umi_reads))
    df_merge_f = df_merge_f.reset_index().groupby(['Barcode', 'umi_21', 'umi_11'])['Read name'].count()
    df_merge_f = df_merge_f.sort_values(ascending=False).reset_index()
    
    return df_merge_f, total_reads, total_umi_reads


# %%
def get_overlap_umi(df_merge_f, umi='umi_21'):
    tmp = df_merge_f.groupby(['Barcode', umi])['Read name'].sum().reset_index()
    out_dict={}
    for i in combinations(tmp['Barcode'].unique(), 2):
        barcode1 = i[0]
        barcode2 = i[1]
        key = "{} vs {}".format(barcode1, barcode2)
        value = (tmp[tmp['Barcode']==barcode2][umi].isin(tmp[tmp['Barcode']==barcode1][umi])).sum()
        if value > 0:
            out_dict[key] = value
    return tmp['Barcode'].value_counts().to_dict(), out_dict

def get_all_sample_umi(df_Lumi):
    out_list = {}
    out_list['sample_reads'] = df_Lumi.groupby('Barcode')['Read name'].sum().to_dict()
    
    sample_umi, sample_overlap_umi = get_overlap_umi(df_Lumi)
    out_list['sample_umi_left'] = sample_umi
    out_list['sample_overlap_umi_left'] = sample_overlap_umi
    sample_umi, sample_overlap_umi = get_overlap_umi(df_Lumi, umi='umi_11')
    out_list['sample_umi_right'] = sample_umi
    out_list['sample_overlap_umi_right'] = sample_overlap_umi
    
    df_out = pd.DataFrame.from_dict(out_list)
    return df_out

# %%
# logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')
# out_dir = '../test/out_dir'
# File_Tag='test'
# L_umi_info_file = os.path.join(out_dir, File_Tag+'_L_umi.csv')
# L_umi_sta = os.path.join(out_dir, File_Tag+'_L_umi_sta.csv')

# umi1_fq = os.path.join(out_dir,File_Tag+'_cut_primer_10.fq')
# umi1_fq_1 = os.path.join(out_dir,File_Tag+'_cut_primer_11.fq')
# umi2_fq = os.path.join(out_dir,File_Tag+'_cut_primer_20.fq')
# umi2_fq_1 = os.path.join(out_dir,File_Tag+'_cut_primer_21.fq')
# df_Lumi, total_reads, total_umi_reads = get_df_merge(umi1_fq, umi1_fq_1, umi2_fq, umi2_fq_1)

# df_out = get_all_sample_umi(df_Lumi)

# L_umi_info_file = os.path.join(out_dir, File_Tag+'_L_umi.csv')
# L_umi_sta = os.path.join(out_dir, File_Tag+'_L_umi_sta.csv')
# df_Lumi.to_csv(L_umi_info_file)
# df_out.to_csv(L_umi_sta)
