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
import os
import logging
from submit_cutadapt import submit_cutadapt
from mkdir import mkdir
from Bio.Seq import Seq


# %%
def cutAadapter(A2, File_Tag='test', out_dir='./',
                Linker=os.path.join('../adapter_fa/', 'L_1-8.fa'),
                Adapter=os.path.join('../adapter_fa/', 'A.fa'),
                cutadapt=os.path.join('../venv/bin', 'cutadapt'),
                error=0,
               ):
    
    mkdir(out_dir)

    A2_cut_adapter = os.path.join(out_dir,File_Tag+'_cut_adapter_P.fq.gz')
    A2_cut_adapter_log = os.path.join(out_dir,File_Tag+'_cut_adapter_P.log')
    A2_cut_adapter_info = A2_cut_adapter + '.cutadapt.info.file'
    A2_umi_file = os.path.join(out_dir,File_Tag+'_cut_linker_P.fq')
    A2_umi_file_log = os.path.join(out_dir,File_Tag+'_cut_linker_P.log')
    A2_umi_file_info = A2_umi_file + '.cutadapt.info.file'
    
    submit_cutadapt(A2,A2_cut_adapter,A2_cut_adapter_log,
                    'file:'+Adapter,'a',cutadapt, error=error, info=1)
    submit_cutadapt(A2_cut_adapter,A2_umi_file,A2_umi_file_log,
                    'file:'+Linker,'g',cutadapt, error=error, info=1)
    
    return A2_cut_adapter_info, A2_umi_file_info, A2_umi_file

# %%
# logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')
# A2 = '../test/rawdata/1-P_R2_test.fastq.gz'
# out_dir = '../test/out_dir'
# cutAadapter(A2, out_dir=out_dir)
