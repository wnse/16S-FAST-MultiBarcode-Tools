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
import multiprocessing
from pathos.multiprocessing import ProcessingPool as Pool
import os
import re
from mkdir import mkdir


def spades_sub(file, outdir, logfile, spades_path):
    '''提交spades组装任务
    参数：
        list: [file, outdir, logfile, spades_path]
            file: 用于组装的fastq序列文件
            outdir: 组装结果输出目录
            logfile: 组装日志输出文件
            spades_path: spades.py程序路径
    返回：
        程序执行的返回值（0,1）
    '''
#     file = file_list[0]
#     outdir = file_list[1]
#     logfile = file_list[2]
#     spades_path = file_list[3]
    cmd = spades_path + ' --s1 ' + file + ' -o ' + outdir + ' -t 1 ' \
          + ' 1>' + logfile + ' 2>' + logfile
    tmp_sys = os.system(cmd)
    rm_cmd = ('/bin/rm -rf outdir/assembly_graph* '
              'outdir/before_rr.fasta '
              'outdir/contigs.paths '
              'outdir/corrected '
              'outdir/dataset.info '
              'outdir/input_dataset.yaml '
              'outdir/misc '
              'outdir/K* '
              'outdir/params.txt '
              'outdir/tmp '
              'outdir/scaffolds.* '
              'outdir/spades.log '
              'outdir/warnings.log '
              'outdir/configs '
              )
    try:
        tmp_sys_rm = os.system(re.sub('outdir', outdir, rm_cmd))
        
    except OSError as e:
        logging.info('rm error {} - {}'.format(e.filename, e.strerror))
    return tmp_sys

def assemble(file_list, out_dir, log_dir, spades_py_path):
    '''多进程提交序列组装任务
    参数：
        file_list: 用于组装的序列文件列表
        out_dir: 组装结果输出目录
        log_dir: 组装日志输出目录
        spades_path: spades.py程序路径
        threads: 进程数
    返回：
        ass_success: 组装成功(结果中有contigs.fasta)数量
        ass_failed: 组装失败数量
    '''
    # tmp_file_list=list(file_list)[:]
    ass_success = 0
    ass_failed = 0
    mkdir(out_dir)
    mkdir(log_dir)
#     queu_file = []
    assemble_outdir = []
    logfile = []
    bin_path = []
    
    for f in file_list:
        f_name = os.path.splitext(os.path.basename(f))[0]
        assemble_outdir.append(os.path.join(out_dir, f_name))
        logfile.append(os.path.join(log_dir, f_name))
        bin_path.append(spades_py_path)
    pool_spade = Pool()
    spades_return = pool_spade.map(spades_sub, file_list, assemble_outdir, logfile, bin_path)
#     spades_return.wait()

    pool_spade.close()
    pool_spade.join()
    pool_spade.clear()
    
    for get in spades_return:
        if get == 0:
            ass_success += 1
        else:
            ass_failed += 1
            
    return ass_success, ass_failed

# %%
# logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')

# out_dir = '../test/out_dir/'
# uIDseq_dir = '../test/out_dir/uIDseq'
# File_Tag = 'test'
# spades_py_path = '../venv/lib/SPAdes-3.12.0-Darwin/bin/spades.py'

# spades_path = os.path.join(out_dir,'spades')
# spades_log_path = os.path.join(out_dir,'spades_log')
# mkdir(spades_path)
# mkdir(spades_log_path)
# tmp_list=[]
# for tmp_dir in os.listdir(uIDseq_dir):
#     spades_path_tmp = os.path.join(spades_path, tmp_dir)
#     spades_log_tmp = os.path.join(spades_log_path, tmp_dir)
#     mkdir(spades_path_tmp)
#     mkdir(spades_log_tmp)
#     file_path_tmp = os.path.join(uIDseq_dir, tmp_dir)
#     files = list(map(lambda x:os.path.join(file_path_tmp,x),os.listdir(file_path_tmp)))
#     logging.info('assemble {} files from {} using {}\t'
#                  .format(len(files), file_path_tmp, spades_py_path))
#     tmp_list.append(assemble(files, spades_path_tmp, spades_log_tmp, spades_py_path))
# assmble_suc = [i[0] for i in tmp_list]
# logging.info(' No.of Contigs assembled:\t{}'.format(sum(assmble_suc)))

# if sum(assmble_suc) == 0:
#     sta_file = os.path.join(result_dir, 'sta.txt')
#     pd.DataFrame(sta_list,columns=(['description', File_Tag])).to_csv(sta_file,
#                                                                       sep='\t',
#                                                                       index=False)
#     del_dir(args.debug)
#     try:
#         sys.exit()
#     finally:
#         print('Exit with no contigs')
