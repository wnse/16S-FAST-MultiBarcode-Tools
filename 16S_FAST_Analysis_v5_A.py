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

# %% tags=[]
import sys
import os
sys.path.append(os.path.join(sys.path[0], './scripts/'))
import logging

import psutil
import shutil
import argparse
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import time

date = time.strftime("%Y%m%d%H%M%S", time.localtime())
def del_dir(debug, ana_dir):
    if debug:
        logging.info(' remove dir {}'.format(ana_dir))
        shutil.rmtree(ana_dir)

script_path = sys.path[0]
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Barcode Analysis')
    parser.add_argument('-c', '--config_file', help='config file including program path', 
                        default=os.path.join(script_path,'config_path_linux.txt'))
    parser.add_argument('-tag','--tag', default='test', help='tag for data')
    parser.add_argument('-ar1','--Assemble_Read_1', required=True, help='Assemble Read 1 fastq data')
    parser.add_argument('-ar2','--Assemble_Read_2', required=True, help='Assemble Read 2 fastq data')
    parser.add_argument('-lr1','--Linker_Read_1', required=True, help='Linker Read 1 fastq data')
    parser.add_argument('-lr2','--Linker_Read_2', required=True, help='Linker Read 2 fastq data')
    parser.add_argument('-FB', '--Fbarcode', help='F barcode list, eg. F1 F2 F3...split by space', required=True, nargs="+")
    parser.add_argument('-RB', '--Rbarcode', help='R barcode list, eg. R1 R2 R3...split by space', required=True, nargs="+")
    parser.add_argument('-N', '--Name', help='Sample names for each barcode,split by space', required=True, nargs="+")
    parser.add_argument('-o', '--out_dir', help='result out dir', required=True)
    parser.add_argument('-d', '--debug', help='for debug', action='store_false')
    parser.add_argument('-remoteA1', '--remoteA1', help='assemble read 1 fastq data in bucket', default=None)
    parser.add_argument('-remoteA2', '--remoteA2', help='assemble read 2 fastq data in bucket', default=None)
    parser.add_argument('-remark', '--remark', help='remark time for date', default=date)
    parser.add_argument('-L_uID_file', '--L_uID_file', help='L_uID_dropRepeat.csv path for uID', required=True)
    
    args = parser.parse_args()
    
    from mkdir import mkdir
    out_dir = args.out_dir
    ana_dir = os.path.join(out_dir, 'analysis')
    mkdir(ana_dir)
    
    logging.basicConfig(filename=os.path.join(out_dir, 'log'), format='%(asctime)s %(levelname)s:%(message)s',level='INFO')
#     logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')

    from submit_Trimmomatic_SE import submit_Trimmomatic_SE
    from plotCover import plotCover
    from submit_cutadapt import submit_cutadapt
    from cutLadapter import cutLadapter
    from staLumi import get_df_merge
    from staLumi import get_all_sample_umi
    from pairLumi import dropRepeatUid
    from pairLumi import pairLumi
    from pairLumi import plot_Lreads
    from cutAadapter import cutAadapter
    from staAumi import staAumi
    from getAseqByBarcode import reads2uID_total
    from plotAumiReads import plotAreads
    from BarcodeAna import BarcodeAna
    from qsctl_sync import qsctl_sync
    from submit_FAST_Barcode import submit_FAST_Barcode

    logging.info(' {}'.format(args.__dict__))
    contig_path = pd.read_csv(args.config_file, sep=':', header=None).set_index(0).to_dict()[1]
    contig_path['threads'] = psutil.cpu_count()
    logging.info(contig_path)    
#     sys.exit()
    File_Tag = args.tag
    L1 = args.Linker_Read_1
    L2 = args.Linker_Read_2
    A2 = args.Assemble_Read_2
    A1 = args.Assemble_Read_1
    
    A1_file = os.path.join(ana_dir, 'A1.fastq.gz')
    logging.info('trim A low quality start')
    mkdir(ana_dir)
    
#    submit_Trimmomatic_SE(A1, A1_file, contig_path['trimmomatic'], contig_path['threads'])

#     L_umi_file = os.path.join(ana_dir, File_Tag+'_L_umi.csv')
#     L_umi_sta_file = os.path.join(out_dir, File_Tag+'_L_umi_sta.csv')
#     L_uID_file = os.path.join(out_dir, File_Tag+'_L_uID_dropRepeat.csv')
#     L_uID_repeat_file = os.path.join(out_dir, File_Tag+'_L_uID_Repeat.csv')
    A_umi_file = os.path.join(out_dir, File_Tag+'_A_umi_dropRepeat.csv')
    A_uID_file = os.path.join(out_dir, File_Tag+'_A_uID_dropRepeat.csv')
    A_uID_repeat_file = os.path.join(out_dir, File_Tag+'_A_uID_Repeat.csv')

#     logging.info('plotCover')
#     plotCover(A1_file, 
#               bowtie_db=contig_path['bowtie_db'], 
#               bowtie_path=contig_path['bowtie_path'],
#               out_dir=ana_dir,
#               out_png=os.path.join(out_dir, File_Tag+'_A1_100000.cover'),
#               threads=contig_path['threads']
#              )
    
    logging.info('cut L adapter&linker')

#     umi1_fq = os.path.join(L_dir,File_Tag+'_cut_primer_10.fq')
#     umi1_fq_1 = os.path.join(L_dir,File_Tag+'_cut_primer_11.fq')
#     umi2_fq = os.path.join(L_dir,File_Tag+'_cut_primer_20.fq')
#     umi2_fq_1 = os.path.join(L_dir,File_Tag+'_cut_primer_21.fq')

#     umi1_fq, umi1_fq_1, umi2_fq, umi2_fq_1 = cutLadapter(L1, L2, out_dir=ana_dir,
#                                                          Linker=contig_path['file_primer'],
#                                                          LinkerRC=contig_path['file_primer_rc'],
#                                                          Adapter=contig_path['file_adapter'],
#                                                          AdapterRC=contig_path['file_adapter_rc'],
#                                                          cutadapt=contig_path['cutadapt_path'],
#                                                         )
    
#     logging.info('get L umi')
#     df_Lumi, total_reads, total_umi_reads = get_df_merge(umi1_fq, umi1_fq_1, umi2_fq, umi2_fq_1)
#     L_umi_sta = get_all_sample_umi(df_Lumi)
#     df_Lumi.to_csv(L_umi_file)
#     L_umi_sta.to_csv(L_umi_sta_file)

    logging.info('get L uID')
#     df_uID, df_uID_repeat = pairLumi(df_Lumi)
#     df_uID.to_csv(L_uID_file)
#     df_uID_repeat.to_csv(L_uID_repeat_file)
#     logging.info('plot L umi&uID')
#     try:
#         plot_Lreads([df_uID, df_uID_repeat], 
#                     [os.path.join(out_dir, File_Tag+'_'+i) for i in ['dropRepeat', 'Repeat']])
#     except Exception as e:
#         logging.info(e)
    df_uID = pd.read_csv(args.L_uID_file, index_col=0)
    
    logging.info('cut A2 adapter&linker')
    A2_cut_adapter_info, A2_umi_file_info, A2_umi_file = cutAadapter(A2, out_dir=ana_dir,
                                                                     Linker=contig_path['file_primer'],
    #                                                                  Linker=os.path.join('adapter_fa/', 'L_1-8.fa'),
                                                                     Adapter=contig_path['file_adapter'],
    #                                                                  Adapter=os.path.join('adapter_fa/', 'A.fa'),
                                                                     cutadapt=contig_path['cutadapt_path'],
    #                                                                  cutadapt=os.path.join('venv/bin', 'cutadapt')
                                                                    )
    
#     A2_cut_adapter = os.path.join(out_dir,File_Tag+'_cut_adapter_P.fq.gz')
#     A2_cut_adapter_info = A2_cut_adapter + '.cutadapt.info.file'
#     A2_umi_file = os.path.join(out_dir,File_Tag+'_cut_linker_P.fq')
#     A2_umi_file_info = A2_umi_file + '.cutadapt.info.file'
    
    logging.info('get A2 umi')
    umi_set = set(df_uID['umi_21'].to_list() + df_uID['umi_11'].to_list())
    df_umi_paired_dropRepeat, df_umi_unpaired_dropRepeat = staAumi(A2_cut_adapter_info, A2_umi_file_info, umi_set)
    df_umi_paired_dropRepeat.to_csv(A_umi_file)
    
    logging.info('get A2 reads2uID')
    df_Aread2uID, df_Aread2uID_repeat = reads2uID_total(df_umi_paired_dropRepeat, df_uID)
    df_Aread2uID.to_csv(A_uID_file)
    
    logging.info('plot A2 umi&uID')
    try:
        plotAreads(df_umi_paired_dropRepeat, 
                   df_umi_unpaired_dropRepeat, 
                   df_Aread2uID, 
                   os.path.join(out_dir, File_Tag)
                  )
    except Exception as e:
        logging.info(e)
        
    del_dir(args.debug, ana_dir)
    remark = args.remark
    remote_path = qsctl_sync(out_dir, remark)
    logging.info('remote path {}'.format(remote_path))
    A_uID_file = os.path.join(remote_path, File_Tag+'_A_uID_dropRepeat.csv')
    A1_file_remote = args.remoteA1
    A2_file_remote = args.remoteA2
    logging.info(A_uID_file)
    logging.info(A1_file_remote)
    if A1_file_remote:
        for i in zip(args.Fbarcode, args.Rbarcode, args.Name):
            Fbarcode = i[0]
            Rbarcode = i[1]
            barcode_name = i[2]
            logging.info('Barcode {} {} analysis'.format(barcode_name, Fbarcode+"_"+Rbarcode))
            result_path = '/mnt/data/yangkai/16S_FAST/V5/'
            submit_FAST_Barcode(A1_file_remote, A_uID_file, Fbarcode, Rbarcode, remark, File_Tag, barcode_name, result_path)
