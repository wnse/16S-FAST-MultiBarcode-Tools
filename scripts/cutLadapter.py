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
import os
import logging

from submit_cutadapt import submit_cutadapt
from mkdir import mkdir
from Bio.Seq import Seq


# %%

# %%
def cutLadapter(L1, L2, File_Tag='test', out_dir='./',
                Linker=os.path.join('../adapter_fa/', 'L_1-8.fa'),
                LinkerRC=os.path.join('../adapter_fa/', 'Lrc_1-8.fa'),
                Adapter=os.path.join('../adapter_fa/', 'A.fa'),
                AdapterRC=os.path.join('../adapter_fa/', 'Arc.fa'),
                cutadapt=os.path.join('../venv/bin', 'cutadapt'),
                error=0,
               ):
    
    mkdir(out_dir)

    L1_cut_adapter = os.path.join(out_dir,File_Tag+'_cut_adapter_10.fq.gz')
    L1_cut_adapter_log = os.path.join(out_dir,File_Tag+'_cut_adapter_10.log')
    L1_cut_adapter_1 = os.path.join(out_dir,File_Tag+'_cut_adapter_11.fq.gz')
    L1_cut_adapter_log_1 = os.path.join(out_dir,File_Tag+'_cut_adapter_11.log')

    umi1_fq = os.path.join(out_dir,File_Tag+'_cut_primer_10.fq')
    umi1_fq_log = os.path.join(out_dir,File_Tag+'_cut_primer_10.log')
    umi1_fq_1 = os.path.join(out_dir,File_Tag+'_cut_primer_11.fq')
    umi1_fq_log_1 = os.path.join(out_dir,File_Tag+'_cut_primer_11.log')

    L2_cut_adapter = os.path.join(out_dir,File_Tag+'_cut_adapter_20.fq.gz')
    L2_cut_adapter_log = os.path.join(out_dir,File_Tag+'_cut_adapter_20.log')
    L2_cut_adapter_1 = os.path.join(out_dir,File_Tag+'_cut_adapter_21.fq.gz')
    L2_cut_adapter_log_1 = os.path.join(out_dir,File_Tag+'_cut_adapter_21.log')

    umi2_fq = os.path.join(out_dir,File_Tag+'_cut_primer_20.fq')
    umi2_fq_log = os.path.join(out_dir,File_Tag+'_cut_primer_20.log')
    umi2_fq_1 = os.path.join(out_dir,File_Tag+'_cut_primer_21.fq')
    umi2_fq_log_1 = os.path.join(out_dir,File_Tag+'_cut_primer_21.log')


    submit_cutadapt(L1,L1_cut_adapter,L1_cut_adapter_log,
                    'file:'+AdapterRC ,'g',cutadapt, error=error, info=1)
    submit_cutadapt(L1_cut_adapter,umi1_fq,umi1_fq_log,
                    'file:'+LinkerRC,'a',cutadapt, error=error, info=1)
    submit_cutadapt(L1, L1_cut_adapter_1, L1_cut_adapter_log_1,
                    'file:'+Adapter, 'a', cutadapt, error=error, info=1)
    submit_cutadapt(L1_cut_adapter_1, umi1_fq_1, umi1_fq_log_1,
                    'file:'+Linker, 'g', cutadapt, error=error, info=1)

    submit_cutadapt(L2,L2_cut_adapter,L2_cut_adapter_log,
                    'file:'+AdapterRC ,'g',cutadapt, error=error, info=1)
    submit_cutadapt(L2_cut_adapter,umi2_fq,umi2_fq_log,
                    'file:'+LinkerRC,'a',cutadapt, error=error, info=1)
    submit_cutadapt(L2, L2_cut_adapter_1, L2_cut_adapter_log_1,
                    'file:'+Adapter, 'a', cutadapt, error=error, info=1)
    submit_cutadapt(L2_cut_adapter_1, umi2_fq_1, umi2_fq_log_1,
                    'file:'+Linker,'g', cutadapt, error=error, info=1)
    return umi1_fq, umi1_fq_1, umi2_fq, umi2_fq_1
