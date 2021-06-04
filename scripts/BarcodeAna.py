# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.9.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import logging
import os
import pandas as pd
import re
from submit_Trimmomatic_SE import submit_Trimmomatic_SE
from submit_Trimmomatic_SE import submit_Trimmomatic_PE
from getAseqByBarcode import reads2uID
from getA1seqByUID import queu_group_umi_seq
from assembleUIDseq import assemble
from filterContigs import filterContigs
from submit_mothur import submit_mothur
from mkdir import mkdir
import shutil
import sys


# %%
def AnnoTax(fa, final_tab,
            mothur_db_fa, mothur_db_tax, mothur, 
            result_dir, threads=4):
    logging.info(' taxonomy start')
    submit_mothur(fa, mothur_db_fa, mothur_db_tax, mothur, threads)
    tax_file = (re.search(r'(.*)\.fasta',fa).group(1) + 
                re.search(r'(\..*)\.tax',mothur_db_tax).group(1) + 
                '.wang.taxonomy')
    tax_ratio_file = os.path.join(result_dir, 'tax.ratio.xlsx')
    tax_reads_file = os.path.join(result_dir, 'tax.reads.xlsx')
    get_tax_from_mothur_with_level(tax_file,
                                   final_tab,
                                   tax_ratio_file,
                                   tax_reads_file)


# %%
def get_tax_from_mothur_with_level(tax_file, asv_file, tax_ratio_file, tax_reads_file):
    '''统计物种分类结果
    参数：
        tax_file: mothur注释结果（.taxonomy文件）
        asv_file: asv或者otu表格（每个asv或者otu的数量）
        tax_ratio_file: 输出文件（物种序列比例）
        tax_reads_file: 输出文件（物种序列数量）
    返回：
        无
    '''
    df_count = pd.DataFrame()
    df_asv_reads = pd.read_csv(asv_file, sep='\t', index_col=0)
    df_asv_ratio = df_asv_reads / df_asv_reads.sum()
    df_count['Contigs'] = df_asv_reads.sum()
    df_count['Asv'] = (df_asv_ratio > 0).sum()
    dict_tag = {'Kingdom': 0,
                'Phylum': 1,
                'Class': 2,
                'Order': 3,
                'Family': 4,
                'Genus': 5,
                'Species': 6}
    dict_tag_r = dict(list((v, k) for k, v in dict_tag.items()))
    df_asv_taxon = pd.read_csv(tax_file,
                               sep='\t',
                               header=None,
                               index_col=0)
    df_asv_taxon.columns = (['Taxon'])
    for pos, tax in dict_tag_r.items():
        tmp_list = list(tax[0].lower() for i in range(df_asv_taxon.shape[0]))
        df_asv_taxon[tax] = tmp_list
        df_asv_taxon[tax] = df_asv_taxon[tax].str.cat(
            df_asv_taxon['Taxon'].str.replace('\(\d+\)', '').str.split(';').str[pos], sep='__')
        if pos > 0:
            df_asv_taxon[tax] = df_asv_taxon[dict_tag_r[pos - 1]].str.cat(df_asv_taxon[tax], sep=';')
    df_asv_taxon.drop('Taxon', axis=1, inplace=True)
    # df_asv_taxon['ID'] = pd.to_numeric(df_asv_taxon['Rep_Seq'].str.split('_').str[1])
    # df_asv_taxon = df_asv_taxon.set_index('ID').drop('Rep_Seq',axis=1)
    df_taxon_ratio = pd.merge(df_asv_ratio,
                              df_asv_taxon,
                              left_index=True,
                              right_index=True)
    df_taxon_reads = pd.merge(df_asv_reads,
                              df_asv_taxon,
                              left_index=True,
                              right_index=True)
    writer_ratio = pd.ExcelWriter(tax_ratio_file)
    writer_reads = pd.ExcelWriter(tax_reads_file)
    sample_name = list(df_asv_ratio.columns)
    for tax in dict_tag.keys():
        sample_name_tmp = list(df_asv_ratio.columns)
        sample_name_tmp.append(tax)
        df_tmp_ratio = df_taxon_ratio[sample_name_tmp].groupby([tax]).sum().reset_index().rename(
            columns={tax: 'Taxonomy'})
        df_tmp_ratio.insert(0, tax, df_tmp_ratio['Taxonomy'].str.split('__').str[-1])
        df_tmp_reads = df_taxon_reads[sample_name_tmp].groupby([tax]).sum().reset_index().rename(
            columns={tax: 'Taxonomy'})
        df_tmp_reads.insert(0, tax, df_tmp_ratio['Taxonomy'].str.split('__').str[-1])
        df_tmp_ratio.sort_values(by=sample_name, ascending=False).reset_index().drop('index', axis=1).to_excel(
            writer_ratio, tax)
        df_tmp_reads.sort_values(by=sample_name, ascending=False).reset_index().drop('index', axis=1).to_excel(
            writer_reads, tax)
        df_count[tax] = (df_tmp_ratio[sample_name] > 0).sum()
    df_count.to_excel(writer_ratio, 'Sta')
    df_count.to_excel(writer_reads, 'Sta')
    writer_ratio.save()
    writer_reads.save()


# %%
def BarcodeAna(A1, A_uID_file, A2=None, Fbarcode='F1', Rbarcode='R3', barcode_name='', out_dir='./', **kargs):
    FB = Fbarcode
    RB = Rbarcode
    if barcode_name=='':
        barcode_name = FB+"_"+RB
    File_Tag_Barcode = out_dir #os.path.join(out_dir, barcode_name)
    ana_dir = os.path.join(File_Tag_Barcode, 'analysis')
    mkdir(File_Tag_Barcode)
    mkdir(ana_dir)
    
    A1_file = os.path.join(ana_dir, 'A1.fastq.gz')
    logging.info('trim A low quality start')
    if A2:
        A2_file = os.path.join(ana_dir, 'A2.fastq.gz')
        submit_Trimmomatic_PE(A1, A2, A1_file, A2_file, kargs['trimmomatic'], kargs['threads'])
    else:
        submit_Trimmomatic_SE(A1, A1_file, kargs['trimmomatic'], kargs['threads'])


    
    uIDseq_dir = os.path.join(ana_dir, barcode_name+"_uIDseq")
    spades_path = os.path.join(ana_dir, barcode_name+'_spades')
    spades_log_path = os.path.join(ana_dir, barcode_name+'_spades_log')

    df_Aread2uID = pd.read_csv(A_uID_file, index_col=0)
    df_Aread2uID = df_Aread2uID[df_Aread2uID['Sample']==str(FB)+'|'+str(RB)][['AreadID', 'umiID']]
    
#     uID_list = df_Aread2uID['umiID'].unique()[0:110]    ############  修改uID数量
#     df_Aread2uID = df_Aread2uID[df_Aread2uID['umiID'].isin(uID_list)]
    logging.info("uID count:\t{}".format(len(df_Aread2uID['umiID'].unique())))
    
    aUMI_dict = df_Aread2uID.set_index('AreadID').to_dict()['umiID']
    Seq2uID_sta = queu_group_umi_seq(A1_file, uIDseq_dir, aUMI_dict)
    if A2:
        Seq2uID_sta = queu_group_umi_seq(A2_file, uIDseq_dir, aUMI_dict)
    
    mkdir(spades_path)
    mkdir(spades_log_path)
    tmp_list=[]
    for tmp_dir in os.listdir(uIDseq_dir):
        spades_path_tmp = os.path.join(spades_path, tmp_dir)
        spades_log_tmp = os.path.join(spades_log_path, tmp_dir)
        mkdir(spades_path_tmp)
        mkdir(spades_log_tmp)
        file_path_tmp = os.path.join(uIDseq_dir, tmp_dir)
        files = list(map(lambda x:os.path.join(file_path_tmp,x),os.listdir(file_path_tmp)))
        logging.info('assemble {} files from {} using {}\t'
                     .format(len(files), file_path_tmp, kargs['spades_py_path']))
        tmp_list.append(assemble(files, spades_path_tmp, spades_log_tmp, kargs['spades_py_path']))
    assmble_suc = [i[0] for i in tmp_list]
    logging.info('Contigs assembled:\t{}'.format(sum(assmble_suc)))
    
    
    if sum(assmble_suc) == 0:
        try:
            sys.exit()
        finally:
            print('Exit with no contigs')
    else:
        final_tab, final_fa, c = filterContigs(spades_path, File_Tag_Barcode,
                                                File_Tag=barcode_name, 
                                                minlength=1200, maxlength=1700,
                                                cdhit=kargs['cdhit_path'],
                                                cutadapt=kargs['cutadapt_path'],
                                                usearch=kargs['usearch_path'],
                                                file_primer_rc=kargs['file_primer_rc'],
                                                file_primer=kargs['file_primer'],
                                                threads=kargs['threads'],
                                               )
        
        if c > 0:
            AnnoTax(final_fa, final_tab, kargs['mothur_db_fa'], 
                    kargs['mothur_db_tax'], kargs['mothur_path'], 
                    File_Tag_Barcode, threads=kargs['threads'])
    shutil.rmtree(ana_dir)


# %%

# %%

# %%
