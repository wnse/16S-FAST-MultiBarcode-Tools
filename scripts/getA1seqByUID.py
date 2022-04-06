# -*- coding: utf-8 -*-
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
import pandas as pd
import logging
import os
import matplotlib.pyplot as plt
from pathos.multiprocessing import ProcessingPool as Pool
from mkdir import mkdir
from Bio import SeqIO
import gzip


# %%
def mkSeqIDFile(seqID, df_Aread2uID, out_dir='./'):
    seqIDfile = os.path.join(out_dir, seqID+".id")
    df_Aread2uID[df_Aread2uID['umiID']==seqID]['AreadID'].to_csv(seqIDfile, index=False, header=False)
    return seqIDfile


def getSeqByID(seqID, rawdata, df_Aread2uID, out_dir='./', seqkit_path='../venv/bin/seqkit'):
    mkdir(out_dir)
    seqIDfile = mkSeqIDFile(seqID, df_Aread2uID, out_dir)
    cmd = "{seqkit} grep --pattern-file {seqIDfile} {rawdata} >{seqID}.fq".format(seqkit=seqkit_path, 
                                                                                  seqIDfile=seqIDfile, 
                                                                                  rawdata=rawdata,
                                                                                  seqID=os.path.join(out_dir, seqID),
                                                                                 )
    try:
        os.system(cmd)
        os.remove(seqIDfile)
        seq_count = len(list(SeqIO.parse(os.path.join(out_dir, seqID+'.fq'), 'fastq')))
        return seq_count
    except Exception as e:
        logging.info(e)
        return 0
    
    
    
def getSeqById_multi(seqID_list, rawdata, df_Aread2uID, out_dir='./', seqkit_path='../venv/bin/seqkit', 
                     func=getSeqByID, dir_files=1000):
    def pool_sub(func, seqID_list_tmp, rawdata_list, df_Aread2uID_list, out_dir_tmp_list, seqkit_path_list):
        seq_count = 0
        pool_write = Pool()
        result = pool_write.map(func, seqID_list_tmp, rawdata_list, df_Aread2uID_list, out_dir_tmp_list, seqkit_path_list)
        pool_write.close()
        pool_write.join()
        pool_write.clear()
        for get in result:
            seq_count += get
        return seq_count
    
    total_seq_count = 0
    n = len(seqID_list)
    range_list = range(int(len(seqID_list)/dir_files))
    for i, v in enumerate(range_list):
        out_dir_tmp = os.path.join(out_dir, 'tmp_'+str(v))
        if i==len(range_list)-1:
            seqID_list_tmp = seqID_list[v*dir_files:]
        else:
            seqID_list_tmp = seqID_list[v*dir_files:range_list[i+1]*dir_files]
        total_seq_count += pool_sub(func, seqID_list_tmp, [rawdata]*n, [df_Aread2uID]*n, [out_dir_tmp]*n, [seqkit_path]*n)
    return total_seq_count


# %%
def write_to_file(file_name, seq_list):
    count = 0
#     file_name = tmp_list[0]
#     file_name_seq_dict = tmp_list[1]
    with open(file_name, 'a') as handle:
        count += SeqIO.write(seq_list, handle, 'fastq')
    return count

def queu_group_umi_seq(seq, out_dir, aUMI_dict):
    '''基于umi将序列分组（多进程）
    参数：
        seq: 输入序列文件（fastq）
        out_dir: 分组序列输出目录
        aUMI_dict: 序列与uID的对应关系（type=dict)
        threads: 进程数
        paired: 默认参数
    返回：
        umi_counts_sub: 每个umi的序列数（type=dict）
        umi_id_counts_sub: 每个umi ID的序列数（type=dict）
        umi2seq_tmp: 分组过程检测的序列数据（type=list）
    '''
#     func = 'group_into_umis'
#     if paired == 0:
#         func = 'group_into_umis_unpaired'

    umi_counts_sub = {}
    umi_id_counts_sub = {}
    umi_file_list_sub = {}
    umi_seq_dir_tmp_sub = {}
    total_check_seq = 0

    file_open = 'open'
    if os.path.splitext(seq)[1] == '.gz':
        file_open = 'gzip.open'

    with eval(file_open)(seq, 'rt') as seq_handle:
        seqs = SeqIO.parse(seq_handle, 'fastq')
        total_umi_seq_dict = {}
        queue = []
        check_umiID = {}
        count = 0
        count_seq = 0
        umi2seq_tmp = []
        for seq_tmp in seqs:
            if seq_tmp.id in aUMI_dict.keys():
                umiID = aUMI_dict[seq_tmp.id]
                count_seq += 1
                total_umi_seq_dict[umiID] = total_umi_seq_dict.get(umiID, [])
                total_umi_seq_dict[umiID].append(seq_tmp)
                total_check_seq += 1
                if (len(total_umi_seq_dict.keys()) > 5000 or total_check_seq >= 500000):
                    file_name_seq_dict = {}
                    queue_write_files = []
                    queue_write_seqs = []
                    for uid, record in total_umi_seq_dict.items():
#                         umiID = umi2ID_dict[u]
                        check_umiID[uid] = check_umiID.get(uid, 0) + len(record)
                        dir_tmp = 'tmp_' + str(int(len(check_umiID.keys()) / 1000))
                        dir_tmp_path = os.path.join(out_dir, dir_tmp)
                        mkdir(dir_tmp_path)
                        umi_seq_dir_tmp_sub[uid] = umi_seq_dir_tmp_sub.get(uid, dir_tmp)
                        file_name = os.path.join(out_dir, umi_seq_dir_tmp_sub[uid], uid + '.fastq')
                        queue_write_files.append(file_name)
                        queue_write_seqs.append(record)
#                         file_name_seq_dict[file_name] = file_name_seq_dict.get(file_name, {})
#                         file_name_seq_dict[file_name] = record
#                     for file_name in file_name_seq_dict.keys():
#                         queue_write.append([file_name, file_name_seq_dict[file_name]])
                    total_umi_seq_dict = {}
#                     file_name_seq_dict = {}
                    total_check_seq = 0
#                     pool_write = multiprocessing.Pool(processes=threads)
#                     result_write = pool_write.map(write_to_file, queue_write_files, queue_write_seqs)
#                     pool_write.close()
#                     pool_write.join()
                    pool_write = Pool()
                    result = pool_write.map(write_to_file, queue_write_files, queue_write_seqs)
                    pool_write.close()
                    pool_write.join()
                    pool_write.clear()
                    for get in result:
                        count += get
                    logging.info(' write {} / {} seqs of {} umi_IDs to files'. \
                                 format(count, count_seq, len(check_umiID.keys())))
                    umi2seq_tmp.append([count, count_seq, len(check_umiID.keys())])
    if total_umi_seq_dict:
        file_name_seq_dict = {}
        queue_write_files = []
        queue_write_seqs = []
        for uid, record in total_umi_seq_dict.items():
            check_umiID[uid] = check_umiID.get(uid, 0) + len(record)
            dir_tmp = 'tmp_' + str(int(len(check_umiID.keys()) / 1000))
            dir_tmp_path = os.path.join(out_dir, dir_tmp)
            mkdir(dir_tmp_path)
            umi_seq_dir_tmp_sub[uid] = umi_seq_dir_tmp_sub.get(uid, dir_tmp)
            file_name = os.path.join(out_dir, umi_seq_dir_tmp_sub[uid], uid + '.fastq')
            queue_write_files.append(file_name)
            queue_write_seqs.append(record)
        total_umi_seq_dict = {}
        total_check_seq = 0
        pool_write = Pool()
        result = pool_write.map(write_to_file, queue_write_files, queue_write_seqs)
        pool_write.close()
        pool_write.join()
        pool_write.clear()
        for get in result:
            count += get
        logging.info(' write {} / {} seqs of {} umi_IDs to files'. \
                     format(count, count_seq, len(check_umiID.keys())))
        umi2seq_tmp.append([count, count_seq, len(check_umiID.keys())])

    return umi2seq_tmp

# %% tags=[]
# logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')

# out_dir = '../test/out_dir_2/test/F1_R3'
# uIDseq_dir = '../test/out_dir_2/test/F1_R3'
# File_Tag = 'F1_R3'
# seqkit_path = '../venv/bin/seqkit'
# # A1_file = '../test/mnt/data/16S_FAST_data/210125_A00262_0590_AHVYVTDSXY/1-P_L2_Y0000085Y0000712.R1.fastq.gz'
# A1_file = '../test/rawdata/1-P_R1_test.fastq.gz'

# A_uID_file = os.path.join(out_dir, File_Tag+'_A_reads2uID.csv')
# df_Aread2uID = pd.read_csv(A_uID_file, index_col=0)[['AreadID', 'umiID']]

# uID_list = df_Aread2uID.head(1000)['umiID'].unique()

# %% tags=[]
# aUMI_dict = df_Aread2uID[df_Aread2uID['umiID'].isin(uID_list)].set_index('AreadID').to_dict()['umiID']
# Seq2uID_sta = queu_group_umi_seq(A1_file, uIDseq_dir, aUMI_dict)

# %% tags=[]
# logging.info("uID count:\t{}".format(len(uID_list)))
# total_seq_count = getSeqById_multi(uID_list, A1_file, df_Aread2uID, uIDseq_dir, seqkit_path, dir_files=10)
# logging.info("total seqs write to uID fq file:\t{}".format(total_seq_count))
