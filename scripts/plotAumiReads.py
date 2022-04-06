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
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
import logging

def plot_barcode_AumiReads(ax, df_list, Barcode_col='Barcode', Barcode='F1', label_list=[]):
    def get_df_label(df):
        tmp_count = df[df[Barcode_col]==Barcode]
        if tmp_count.empty:
            logging.info('{} not in here'.format(Barcode))
            return tmp_count['AreadID']
        else:
            tmp_count = tmp_count['AreadID'].sort_values(ascending=False)
            return tmp_count
    
    n = len(df_list)
    color_list = ['blue','red','green']
    describe_data = []
    title = ' : '.join(['sum','count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max'])
    ax = ax or plt.gca()
    for i, df in enumerate(df_list):
        if label_list:
            label_tmp = label_list[i]
        else:
            label_tmp = ''
        df_tmp = get_df_label(df)
        if not df_tmp.empty:
            df_tmp.replace(np.nan, 0, inplace=True)
            df_tmp.replace(np.inf, 0, inplace=True)
            ax.plot(df_tmp.to_list(), color=color_list[i], label=label_tmp)
            sta_tmp = df_tmp.describe().round(1).astype(str).to_list() #astype(int)
            sta_tmp = [str(np.sum(df_tmp))] + sta_tmp
            describe_data.append(sta_tmp)
    text_x = ax.get_xlim()[1]*0.05
    text_y = ax.get_ylim()[1]*0.95
    ax.annotate(title, (text_x, text_y), (0, 2), textcoords='offset points')
    ax.set_title(Barcode)
    for i, s in enumerate(describe_data):
        logging.info("{} describe {}".format(Barcode, " : ".join(s)))
        ax.annotate(" : ".join(s), (text_x, text_y), (0, -15*(i+1)), textcoords='offset points', c=color_list[i])

def plot_barcode_AumiReads_total(df_list, out_png, Barcode_list=[], Barcode_col='Barcode'):
    M = 4
    N = math.ceil(len(Barcode_list)/M)
    if M>len(Barcode_list):
        M = len(Barcode_list)
        N = 1
    fig, ax = plt.subplots(N, M, sharex=True, sharey=True, figsize=(40,22.5))
    m = 0
    n = 0
    barcode_no = 0
    for i in range(N):
        for j in range(M):
            if M==len(Barcode_list):
                ax_tmp = ax[j]
            else:
                ax_tmp = ax[i,j]
            plot_barcode_AumiReads(ax_tmp, 
                                   df_list,
                                   Barcode_col=Barcode_col,
                                   Barcode=Barcode_list[barcode_no], 
                                   label_list=[]
                                  )
            m+=1
            if m==M:
                n+=1
                m=0
            barcode_no += 1
            if barcode_no >= len(Barcode_list):
                fig.savefig(out_png)
                return 0

def plot_barcode_AumiReads_box(df_list, out_png, Barcode_col='Barcode', label_col='label', Barcode_list=[]):
    tmp_data = pd.concat(df_list)
    fig, ax = plt.subplots(figsize=(40,22.5))
    if Barcode_list==[]:
        Barcode_list = tmp_data[Barcode_col].unique()
    if label_col in tmp_data.columns:
        sns.boxplot(x=Barcode_col, y='AreadID', hue=label_col, data=tmp_data, color="1", order=Barcode_list, dodge=True, ax=ax)
        sns.stripplot(x=Barcode_col, y='AreadID', hue=label_col, data=tmp_data, order=Barcode_list, dodge=True, ax=ax)
    else:
        sns.boxplot(x=Barcode_col, y='AreadID', data=tmp_data, color="1", order=Barcode_list, ax=ax)
        sns.stripplot(x=Barcode_col, y='AreadID', data=tmp_data, order=Barcode_list, ax=ax)
#     plt.show()
    fig.savefig(out_png)

def plotAreads(df_umi_paired_dropRepeat, df_umi_unpaired_dropRepeat, df_Aread2uID, outfile):
    def get_count(df, col1='Barcode', col2='umi'):
        df_count = df.groupby([col1, col2])['AreadID'].count()
        logging.info('A {}&{} describe:\t{}'.format(col1, col2, ":".join(df_count.describe().round(2).astype(str).to_list())))
        return df_count.reset_index()
    df_Aumi_count_unpaired = get_count(df_umi_unpaired_dropRepeat)
    df_Aumi_count = get_count(df_umi_paired_dropRepeat)
    
#     df_AuID_count = df_Aread2uID.groupby(['Sample','umiID'])['AreadID'].count().reset_index()
    Barcode_list = df_umi_paired_dropRepeat['Barcode'].value_counts().index.to_list()
    plot_barcode_AumiReads_total([df_Aumi_count, df_Aumi_count_unpaired], 
                                 outfile+'_Aumi_counts.png', 
                                 Barcode_list=Barcode_list, 
                                 Barcode_col='Barcode'
                                )
    
    tmp_data1 = df_Aumi_count.copy()
    tmp_data1['label'] = 'InL'
    tmp_data2 = df_Aumi_count_unpaired.copy()
    tmp_data2['label'] = 'NotInL'
    plot_barcode_AumiReads_box([tmp_data1, tmp_data2], 
                               outfile+'_A_umi_counts_box.png', 
                               Barcode_col='Barcode', 
                               label_col='label',
                               Barcode_list=Barcode_list,
                              )
    
    samples = df_Aread2uID['Sample'].value_counts().head(20).index.to_list()
    df_AuID_count = get_count(df_Aread2uID, col1='Sample', col2='umiID')
    tmp_data = df_AuID_count[df_AuID_count['Sample'].isin(samples)]
    plot_barcode_AumiReads_total([tmp_data],
                                 outfile+'_A_uID_counts.png', 
                                 Barcode_list=samples, 
                                 Barcode_col='Sample')
    
    plot_barcode_AumiReads_box([tmp_data],
                               outfile+'_A_uID_counts_box.png',
                               Barcode_col='Sample',
                               Barcode_list=samples, 
                              )
    return df_Aumi_count, df_Aumi_count_unpaired, df_AuID_count

# %%
# import os
# logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')

# out_dir = '/Users/yk/work/Microbiology/16S-FAST/16S_FAST_result/Project/YF/20210128_multisample/20210128/NP_1_DS/DS0/'
# File_Tag='DS0'
# A_umi_file = os.path.join(out_dir, File_Tag+'_A_umi_dropRepeat.csv')
# A_uID_file = os.path.join(out_dir, File_Tag+'_A_uID_dropRepeat.csv')
# A_umi_repeat_file = os.path.join(out_dir, File_Tag+'_A_umi_unpaired_dropRepeat.csv')

# df_umi_paired_dropRepeat = pd.read_csv(A_umi_file, index_col=0)
# df_umi_unpaired_dropRepeat = pd.read_csv(A_umi_repeat_file, index_col=0)
# df_Aread2uID = pd.read_csv(A_uID_file, index_col=0)
# df_Aumi_count, df_Aumi_count_unpaired, df_AuID_count = plotAreads(df_umi_paired_dropRepeat, 
#                                                                   df_umi_unpaired_dropRepeat, 
#                                                                   df_Aread2uID, 
#                                                                   out_dir+'test_DS0_fig',
#                                                                  )

# df_Aumi_count_unpaired = df_umi_unpaired_dropRepeat.groupby(['Barcode','umi'])['AreadID'].count().reset_index()
# df_Aumi_count = df_umi_paired_dropRepeat.groupby(['Barcode','umi'])['AreadID'].count().reset_index()
# Barcode_list = df_umi_paired_dropRepeat['Barcode'].value_counts().index
# plot_barcode_AumiReads_total([df_Aumi_count, df_Aumi_count_unpaired], 
#                              'Aumi_counts.png', 
#                              Barcode_list=Barcode_list, 
#                              Barcode_col='Barcode'
#                             )
