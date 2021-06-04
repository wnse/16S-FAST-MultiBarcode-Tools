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
import re
import time

def qsctl_sync(local_dir, remark='', bucket='qs://bc-input-prod', bucket_path='/mnt/data/16S_FAST_V5_UMI/'):
    if not remark:
        remark = time.strftime("%Y%m%d%H%M%S", time.localtime()) 
    remote_path = os.path.join(bucket_path, remark, os.path.split(re.sub('/$', '', local_dir))[1])
    remote = bucket+remote_path
    cmd = ' '.join(['qsctl sync', local_dir, remote])
    logging.info(cmd)
    status = os.system(cmd)
    return remote_path
    #qsctl sync local_dir 'qs://bc-input-prod/mnt/data/16S_FAST_V5_UMI/test1/

