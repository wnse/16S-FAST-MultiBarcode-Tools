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
import argparse
import logging
import os
import pandas as pd
from Bio import SeqIO
import re
from submit_cutadapt import submit_cutadapt
from get_consensus_seq_from_cdhit import get_consensus_seq_from_cdhit
from remove_chimeras_uchime import remove_chimeras_uchime
from mkdir import mkdir

def merge_contigs(contig_list, merge_fa_name):
    '''合并组装后contigs
    参数：
        contig_list: contig文件路径列表（type=list）
        merge_fa_name: 输出文件
    返回：
        df: 每个umiID的contig数目（type=pandas）
    '''
    seq_no_dict = {}
    with open(merge_fa_name, 'w') as w:
        m = 0
        for record in contig_list:
            name = record[0]
            fa = record[1]
            tmp = []
            n = 0
            for seq in SeqIO.parse(fa, 'fasta'):
                seq.description = str(name) + "|" + str(seq.description)
                seq.id = str(m) + '_' + str(n)
                tmp.append(seq)
                n += 1
            seq_no = SeqIO.write(tmp, w, 'fasta')
            seq_no_dict[name] = [m, seq_no]
            m += 1
    return pd.DataFrame.from_dict(seq_no_dict, orient='index',
                                  columns=(['umi_id', 'NO_of_Contigs']))

def cut_fa_by_len(input_fa, output_fa, min_len=0, max_len=0):
    '''筛选fasta文件符合长度要求的序列
    参数：
        input_fa: 输入fasta文件
        output_fa: 输出fasta文件
        min_len: 序列最小长度
        max_len: 序列最大长度
    返回：
        c: 过滤后序列数量
    '''
    tmp = []
    if max_len == 0:
        for seq in SeqIO.parse(input_fa, 'fasta'):
            l = len(seq.seq)
            if l >= min_len:
                seq.description = 'len=' + str(l)
                tmp.append(seq)
    else:
        for seq in SeqIO.parse(input_fa, 'fasta'):
            l = len(seq.seq)
            if l >= min_len and l < max_len:
                seq.description = 'len=' + str(l)
                tmp.append(seq)
    c = SeqIO.write(tmp, output_fa, 'fasta')
    return c

def submit_cdhit(infa, outfa, identity, path, threads=1):
    '''序列聚类
    参数：
        infa: 输入fasta序列
        outfa: 输出fasta代表序列
        identity: 相似度
        path: cd-hit路径
        threads: 线程数
    返回：
        无
    '''
    cdhit = path
    str_join = ' '
    cmd = str_join.join([cdhit,
                         '-i', infa,
                         '-o', outfa,
                         '-c', str(identity),
                         '-T', str(threads),
                         '-M', str(0)]
                        )
    logging.info(' %s' % cmd)
    status = os.system(cmd)
    if status == 0:
        logging.info(' done')
    else:
        logging.info(' exit')
        
def update_seqID(infa, outfa, tag):
    '''重新命名序列
    参数：
        infa: 输入fasta文件
        outfa: 输出fasta文件
        tag: 序列标签
    返回：
        id_track: 序列ID对应关系（type=dict）
    '''
    n = 0
    id_track = {}
    seqs = []
    for rec in SeqIO.parse(infa, 'fasta'):
        new_id = tag + '_' + str(n)
        id_track[new_id] = rec.id
        rec.id = new_id
        rec.name = new_id
        rec.description = str(rec.description).split()[-1]
        n += 1
        seqs.append(rec)
    count = SeqIO.write(seqs, outfa, 'fasta')
    return id_track

def get_final_contigs_from_id(rawfa, idlist, outfile, prefix='contig'):
    raw=[]
    if os.path.splitext(rawfa)[1] in ['.fa','.fasta','.FASTA', '.Fasta', '.FA']:
        raw = list(SeqIO.parse(rawfa, 'fasta'))
    else:
        try:
            with gzip.open(rawfa, 'rt') as h:
                raw = list(SeqIO.parse(h, 'fasta'))
        except Exception as e:
            print('error: {}'.format(e))
    contigs = [i for i in raw if i.id in idlist]
    for i in contigs:
        i.id = str(prefix) + '_' + str(i.id)
        i.description = ''
    count = SeqIO.write(contigs, outfile, 'fasta')
    return count


# %%
def filterContigs(contigs_dir, output_dir, 
                  File_Tag='test',
                  minlength=1200, maxlength=1700, 
                  cdhit='cd-hit',
                  cutadapt='cutadapt',
                  usearch='usearch11',
                  file_primer_rc='adapter_fa/Lrc_1-8.fa',
                  file_primer='adapter_fa/L_1-8.fa',
                  threads=2,
                 ):
    spades_path = contigs_dir
    mkdir(output_dir)
    result_dir = output_dir
    ana_dir = output_dir
    merge_fa = os.path.join(result_dir, File_Tag+'.merged.fasta')
    merge_trim_fa1 = os.path.join(ana_dir, File_Tag+'.merged.trim1.fasta')
    merge_trim_log1 = os.path.join(ana_dir, File_Tag+'.merged.trim1.log')
    merge_trim_fa2 = os.path.join(result_dir, File_Tag+'.merged.trim.fasta')
    merge_trim_log2 = os.path.join(ana_dir, File_Tag+'.merged.trim2.log')
    merge_filter_fa = os.path.join(ana_dir, File_Tag+'.merged.filter.fasta')
    merge_filter_cluster_fa = os.path.join(ana_dir, File_Tag+'.merged.filter.rep.fa')
    clust_fa = os.path.join(ana_dir, File_Tag+'.clust.fasta')
    clust_fa_clstr_table = os.path.join(result_dir, File_Tag+'.clust.clstr.info')
    clust_uchimeout = os.path.join(ana_dir, File_Tag+'.clust.uchimeout.out.txt')
    clust_ch_fa = os.path.join(ana_dir, File_Tag+'.clust.ch.fasta')
    clust_nonch_fa = os.path.join(ana_dir, File_Tag+'.clust.nonch.fasta')
    final_fa = os.path.join(result_dir, File_Tag+'.final.fasta')
    final_tab = os.path.join(result_dir, File_Tag+'.final.size.tab.txt')
    ID_info = os.path.join(result_dir, File_Tag+'.final.id.info')
    ID_contig_info = os.path.join(result_dir, File_Tag+'.final.contig.id.info')
    final_contig_fa = os.path.join(result_dir, File_Tag+'.final.contig.fasta')

    #合并组装后的contigs
    assemble_list = []
    for tmpdir in os.listdir(spades_path):
        for d in os.listdir(os.path.join(spades_path, tmpdir)):
            contig_file = os.path.join(spades_path, tmpdir, d, 'contigs.fasta')
            if os.path.exists(contig_file):
                tmp = [d,contig_file]
                assemble_list.append(tmp)
    df_merge_fa = merge_contigs(assemble_list,merge_fa)
    
    #去除contig中引物序列
    logging.info('remove primer in contigs')
    submit_cutadapt(merge_fa, 
                    merge_trim_fa1, 
                    merge_trim_log1,
                    "file:"+file_primer_rc,
                    'a',
                    cutadapt,
                    threads=threads
                   )
    submit_cutadapt(merge_trim_fa1,
                    merge_trim_fa2,
                    merge_trim_log2,
                    "file:"+file_primer,
                    'g',
                    cutadapt,
                    threads=threads
                   )
    #长度过滤 大于minlength小于maxlength
    tmp = cut_fa_by_len(merge_trim_fa2,
                        merge_filter_fa,
                        minlength,
                        maxlength)
    logging.info('Contigs after length filtered({}-{}):\t{}'.format(minlength, maxlength, tmp))
    
    #聚类 100%相似度 cd-hit
    logging.info('contig cluster start')
    submit_cdhit(merge_filter_fa, 
                 merge_filter_cluster_fa, 
                 1, 
                 cdhit, 
                 threads)
    
    #去除同一umiID中未聚类在同一组的contig
    logging.info('get consensus seq')
    consensus_seq_count, rep_seq_tab = get_consensus_seq_from_cdhit(
        merge_filter_cluster_fa+'.clstr', 
        merge_filter_fa, 
        clust_fa)
    
    pd.DataFrame(rep_seq_tab,
                 columns=(['clust_rep_id','size','seq_id'])).to_csv(clust_fa_clstr_table, 
                                                                  sep='\t',
                                                                  index=False)
    #去除嵌合
    logging.info('remove cluster chimeras start')
    remove_chimeras_uchime(clust_fa,
                           clust_uchimeout,
                           clust_ch_fa,
                           clust_nonch_fa,
                           usearch)
    dict_update_id = update_seqID(clust_nonch_fa, final_fa, File_Tag)
    
    # 生成cluster size表格
    tab_list = []
    total_size = 0
    total_cluster = 0
    for i in SeqIO.parse(final_fa, 'fasta'):
        size = re.search('size=(\d+)', i.description).group(1)
        tab_list.append([i.id, size])
        total_size += int(size)
        total_cluster += 1
    pd.DataFrame(tab_list, columns=(['id', File_Tag])).set_index('id').to_csv(final_tab, sep='\t')
    logging.info('Cluster after remove Chimerias:\t{}'.format(total_cluster))
    logging.info('Contigs after remove Chimerias:\t{}'.format(total_size))
    
    #合并并输出代表序列对应的ID信息
    df_merge_fa = df_merge_fa.reset_index().rename(columns={'index': 'umiID'})
    df_update_id = pd.DataFrame.from_dict(dict_update_id,
                                          orient='index',
                                          columns=(['contig_id']))
    df_update_id['umi_id'] = pd.to_numeric(df_update_id['contig_id'].str.split('_').str[0])
    df_ID_info = pd.merge(df_update_id.reset_index().rename(columns={'index': 'ID'}),
                          df_merge_fa,
                          left_on='umi_id',
                          right_on='umi_id',
                          how='left')
    df_ID_info.to_csv(ID_info, sep='\t')
    #提取并输出所有cluster的contigs对应的umiID信息
    df_tmp = df_merge_fa[['umi_id','umiID']].set_index('umi_id')
    # df_tmp.to_csv(os.path.join(result_dir, 'df_merge_fa.csv'))
    total_id = df_tmp.to_dict()['umiID']
    final_clust = df_update_id['contig_id'].to_list()
    final_total_contig = {}
    for c in rep_seq_tab:
        if c[0] in final_clust:
            for i in c[2].split('|'):
                umiid = int(i.split('_')[0])
                final_total_contig[i] = total_id[umiid]
    c = get_final_contigs_from_id(merge_trim_fa2,
                                  final_total_contig.keys(),
                                  final_contig_fa,
                                  prefix=File_Tag)
    return final_tab, final_fa, c

# %%
# logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')

# contigs_dir = '../test/out_dir_2/test/F1_R3_spades'
# output_dir = '../test/out_dir_2/test/F1_R3/'
# File_Tag = 'F1_R3'
# filterContigs(contigs_dir, output_dir, 
#               File_Tag=File_Tag, minlength=1200, maxlength=1700, 
#               cdhit='../venv/bin/cd-hit',
#               cutadapt='../venv/bin/cutadapt',
#               file_primer_rc='../adapter_fa/Lrc_1-8.fa',
#               file_primer='../adapter_fa/L_1-8.fa',
#               threads=2,
#              )
