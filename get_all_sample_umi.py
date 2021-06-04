import pandas as pd
import os
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
import numpy as np

def rc_seq(seq):
    if seq == seq:
        return str(Seq(seq,generic_dna).reverse_complement())
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

def get_df_umi(ana_dir, File_Tag):    
    umi1_fq = os.path.join(ana_dir,File_Tag+'_cut_primer_10.fq')
    umi1_fq_log = os.path.join(ana_dir,File_Tag+'_cut_primer_10.log')
    umi1_fq_1 = os.path.join(ana_dir,File_Tag+'_cut_primer_11.fq')
    umi1_fq_log_1 = os.path.join(ana_dir,File_Tag+'_cut_primer_11.log')
    
    umi2_fq = os.path.join(ana_dir,File_Tag+'_cut_primer_20.fq')
    umi2_fq_log = os.path.join(ana_dir,File_Tag+'_cut_primer_20.log')
    umi2_fq_1 = os.path.join(ana_dir,File_Tag+'_cut_primer_21.fq')
    umi2_fq_log_1 = os.path.join(ana_dir,File_Tag+'_cut_primer_21.log')
    
    
    suffix = '.cutadapt.info.file'
    
    ul = [10,11,20,21]
    df_10 = get_umi_info(umi1_fq+suffix)
    df_11 = get_umi_info(umi1_fq_1+suffix, col='right')
    df_20 = get_umi_info(umi2_fq+suffix)
    df_21 = get_umi_info(umi2_fq_1+suffix, col='right')
    
    return df_10, df_11, df_20, df_21



def get_df_merge(ana_dir, File_Tag):
    df_10, df_11, df_20, df_21 = get_df_umi(ana_dir, File_Tag)
    df_11['umi_11'] = df_11['umi'].apply(rc_seq)
    df_21['umi_21'] = df_21['umi'].apply(rc_seq)
    df_merge = df_10.rename(columns={'Number of errors':'errors_10', 'umi':'umi_10', 'Name':'Name_10'})
    df_tmp = df_11.rename(columns={'Number of errors':'errors_11','Name':'Name_11'}).drop('umi', axis=1)
    df_merge = pd.merge(df_merge, df_tmp, left_index=True, right_index=True)
    df_tmp = df_20.rename(columns={'Number of errors':'errors_20', 'umi':'umi_20','Name':'Name_20'})
    df_merge = pd.merge(df_merge, df_tmp, left_index=True, right_index=True)
    df_tmp = df_21.rename(columns={'Number of errors':'errors_21','Name':'Name_21'}).drop('umi', axis=1)
    df_merge = pd.merge(df_merge, df_tmp, left_index=True, right_index=True)
    total_reads = df_merge.shape[0]
    print('total reads:{}'.format(total_reads))
    df_merge_f = df_merge[(df_merge[df_merge.columns[df_merge.columns.str.contains('errors')]]==0).all(axis=1)]
    print('cut adapter or linkers:{}'.format(df_merge_f.shape[0]))

    df_merge_f = df_merge_f[((df_merge_f['umi_10'] == df_merge_f['umi_21']) & (df_merge_f['umi_11'] == df_merge_f['umi_20']))]
    print('R1_F equal R2_R and reverse:{}'.format(df_merge_f.shape[0]))

    df_merge_f = df_merge_f[((df_merge_f['umi_10'].apply(len) == 14) & 
                             (df_merge_f['umi_11'].apply(len) == 14) &
                             (df_merge_f['umi_20'].apply(len) == 14) & 
                             (df_merge_f['umi_21'].apply(len) == 14))
                           ]
    umis_reads = df_merge_f.shape[0]
    print('14bp left:{}'.format(umis_reads))
    return df_merge_f, total_reads, umis_reads

def get_overlap_umi(df_merge_f, umi='umi_10'):
    tmp = df_merge_f.reset_index().groupby(['R1_Barcode', umi])['Read name'].count().reset_index()
    print(tmp['R1_Barcode'].value_counts())
    out_dict={}
    for i in combinations(tmp['R1_Barcode'].unique(), 2):
        barcode1 = i[0]
        barcode2 = i[1]
        key = "{} vs {}".format(barcode1, barcode2)
        value = (tmp[tmp['R1_Barcode']==barcode2][umi].isin(tmp[tmp['R1_Barcode']==barcode1][umi])).sum()
        out_dict[key] = value
    return tmp['R1_Barcode'].value_counts().to_dict(), out_dict

def get_all_sample_umi(ana_dir, File_Tag):
#     ana_dir = 'cut_adapter_out/'
#     File_Tag = '15'
    out_list = {}
    # df_10, df_11, df_20, df_21 = get_df_umi(os.path.join(ana_dir,File_Tag), File_Tag)
    out_file = File_Tag+'.umi_info.txt'
    df_merge_f, total_reads, umis_reads = get_df_merge(os.path.join(ana_dir,File_Tag), File_Tag)
    out_list['total_reads'] = total_reads
    out_list['umi_reads'] = umis_reads
    
    df_merge_f['R1_Barcode'] = df_merge_f['Name_10'] + '|' + df_merge_f['Name_11']
    df_merge_f['R2_Barcode'] = df_merge_f['Name_21'] + '|' + df_merge_f['Name_20']
    df_merge_f.to_csv(out_file)
    out_list['sample_reads'] = df_merge_f['R1_Barcode'].value_counts().to_dict()
    
    sample_umi, sample_overlap_umi = get_overlap_umi(df_merge_f)
    out_list['sample_umi_left'] = sample_umi
    out_list['sample_overlap_umi_left'] = sample_overlap_umi
    sample_umi, sample_overlap_umi = get_overlap_umi(df_merge_f, umi='umi_11')
    out_list['sample_umi_right'] = sample_umi
    out_list['sample_overlap_umi_right'] = sample_overlap_umi
    
    df_out = pd.DataFrame.from_dict(out_list)
    df_out.columns = pd.MultiIndex.from_product([[File_Tag], df_out.columns])
    return df_out

ana_dir = 'cut_adapter_out/'
df_total = pd.DataFrame()
for i in os.listdir(ana_dir):
    File_Tag = i
    print(File_Tag)
    try:
        df_tmp = get_all_sample_umi(ana_dir, File_Tag)
        print(df_tmp)
        if df_total.empty:
            df_total = df_tmp
        else:
            df_total = pd.concat([df_total, df_tmp], axis=1)
    except Exception as e:
        print(e)
df_total.to_csv('cut_adapter_out_'+'umi_info.txt')
