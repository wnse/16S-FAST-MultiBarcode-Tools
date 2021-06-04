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
import argparse
import logging
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# %%
def reads2uID(df_umi_paired_dropRepeat, df_uID_dropRepeat, Fbarcode='F1', Rbarcode='R3'):
    df_uID = df_uID_dropRepeat[df_uID_dropRepeat['Barcode']==Fbarcode+'|'+Rbarcode]
    umi_list = set(df_uID['umi_21'].to_list() + df_uID['umi_11'].to_list())
    df_umi = df_umi_paired_dropRepeat[df_umi_paired_dropRepeat['umi'].isin(umi_list)]
    
    df_Aread2uID = pd.merge(df_umi, df_uID.set_index('umi_21')['umiID'], left_on='umi', right_index=True, how='left')
    df_Aread2uID = pd.merge(df_Aread2uID, df_uID.set_index('umi_11')['umiID'], left_on='umi', right_index=True, how='left')

    check1 = df_Aread2uID[df_Aread2uID[['umiID_x', 'umiID_y']].isna().all(axis=1)].index
    check2 = df_Aread2uID[(~df_Aread2uID[['umiID_x', 'umiID_y']].isna()).all(axis=1)].index
    logging.info('umis not in any uID of barcode {}:\t{}'.format(Fbarcode+"|"+Rbarcode, len(check2)))
    logging.info('umis in multiple uID of barcode {}:\t{}'.format(Fbarcode+"|"+Rbarcode, len(check2)))

    if len(check1) != 0 or  len(check2) != 0:
        df_Aread2uID = df_Aread2uID.drop(np.concatenate([check1, check2]))

    df_Aread2uID.loc[df_Aread2uID['umiID_x'].isna(), 'umiID'] = df_Aread2uID.loc[df_Aread2uID['umiID_x'].isna(), 'umiID_y']
    df_Aread2uID.loc[df_Aread2uID['umiID_y'].isna(), 'umiID'] = df_Aread2uID.loc[df_Aread2uID['umiID_y'].isna(), 'umiID_x']
    df_Aread2uID = df_Aread2uID.drop(['umiID_x', 'umiID_y'], axis=1)
    df_Aread2uID = df_Aread2uID.drop_duplicates()
    
    return df_Aread2uID


# %%
def reads2uID_total(df_umi_paired_dropRepeat, df_uID):
    def check_BinS(df_Aread2uID, s='Sample_x'):
        df_tmp = df_Aread2uID[~(df_Aread2uID[s].isna())]
        idx = df_tmp.apply(lambda df: df['Barcode'] in df[s].split('|'), axis=1)
        if len(idx)>0:
            logging.info('\n{}'.format(df_tmp[~idx]))
        return df_tmp[~idx].index
        
#     df_uID = df_uID_dropRepeat[df_uID_dropRepeat['Barcode']==Fbarcode+'|'+Rbarcode]
    umi_list = set(df_uID['umi_21'].to_list() + df_uID['umi_11'].to_list())
    df_umi = df_umi_paired_dropRepeat[df_umi_paired_dropRepeat['umi'].isin(umi_list)]
    df_uID_tmp = df_uID.set_index('umi_21')[['umiID','Barcode']].rename(columns={'Barcode':'Sample'})
    df_Aread2uID = pd.merge(df_umi, df_uID_tmp, left_on='umi', right_index=True, how='left')
    df_uID_tmp = df_uID.set_index('umi_11')[['umiID','Barcode']].rename(columns={'Barcode':'Sample'})
    df_Aread2uID = pd.merge(df_Aread2uID, df_uID_tmp, left_on='umi', right_index=True, how='left')

    check1 = df_Aread2uID[df_Aread2uID[['umiID_x', 'umiID_y']].isna().all(axis=1)].index
    check2 = df_Aread2uID[(~df_Aread2uID[['umiID_x', 'umiID_y']].isna()).all(axis=1)].index
    logging.info('\n{}'.format(df_Aread2uID.loc[check1]))
    logging.info('\n{}'.format(df_Aread2uID.loc[check2]))
    logging.info('umis not in any uID:\t{}'.format(len(check1)))
    logging.info('umis in multiple uID:\t{}'.format(len(check2)))

    df_Aread2uID.loc[df_Aread2uID['umiID_x'].isna(), 'umiID'] = df_Aread2uID.loc[df_Aread2uID['umiID_x'].isna(), 'umiID_y']
    df_Aread2uID.loc[df_Aread2uID['umiID_y'].isna(), 'umiID'] = df_Aread2uID.loc[df_Aread2uID['umiID_y'].isna(), 'umiID_x']
        
    check3 = check_BinS(df_Aread2uID, s='Sample_x')
    check4 = check_BinS(df_Aread2uID, s='Sample_y')
    logging.info('A2 Barcode inconsistent with L F barcode:\t{}'.format(len(check3)))
    logging.info('A2 Barcode inconsistent with L R barcode:\t{}'.format(len(check4)))
    
    df_Aread2uID.loc[df_Aread2uID['Sample_x'].isna(), 'Sample'] = df_Aread2uID.loc[df_Aread2uID['Sample_x'].isna(), 'Sample_y']
    df_Aread2uID.loc[df_Aread2uID['Sample_y'].isna(), 'Sample'] = df_Aread2uID.loc[df_Aread2uID['Sample_y'].isna(), 'Sample_x']


    if len(check1) != 0 or  len(check2) != 0 or len(check3) != 0 or len(check4) != 0:
        idx = np.concatenate([check1, check2, check3, check4])
        df_Aread2uID_repeat = df_Aread2uID.loc[idx]
        df_Aread2uID = df_Aread2uID.drop(idx)
    else:
        df_Aread2uID_repeat = pd.DataFrame()
    
    df_Aread2uID = df_Aread2uID.drop(['umiID_x', 'umiID_y'], axis=1)
    df_Aread2uID = df_Aread2uID.drop(['Sample_x', 'Sample_y'], axis=1)
    df_Aread2uID = df_Aread2uID.drop_duplicates()
    return df_Aread2uID, df_Aread2uID_repeat


# df_Aread2uID, df_Aread2uID_repeat = reads2uID_total(df_umi_paired_dropRepeat, df_uID)

# %%
# out_dir = '../test/out_dir_2'
# File_Tag = 'test'
# FB = 'F1'
# RB = 'R3'
# # A_uID_file = os.path.join(out_dir, File_Tag+'_A_reads2uID.csv')


# A_umi_info_file = os.path.join(out_dir, File_Tag+'_A_umi_dropRepeat.csv')
# df_umi_paired_dropRepeat = pd.read_csv(A_umi_info_file, index_col=0)
# L_uID_info_file = os.path.join(out_dir, File_Tag+'_L_uID_dropRepeat.csv')
# df_uID_dropRepeat = pd.read_csv(L_uID_info_file, index_col=0)
# df_Aread2uID, df_Aread2uID_repeat = reads2uID_total(df_umi_paired_dropRepeat, df_uID_dropRepeat)
# # df_Aread2uID.to_csv(A_uID_file)
