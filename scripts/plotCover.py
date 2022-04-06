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
import pandas as pd
import re
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from submit_bowtie import submit_bowtie

def get_map_info(map_sam_list):
    map_info = []
    for i in map_sam_list:
        if re.match('@', i):
            next
        else:
            info = i.split('\t')
            if len(info)<5:
                next
            else:
                if re.search('(\d+)M', info[5]):
                    map_len = re.search('(\d+)M', info[5]).group(1)
                    map_pos = info[3]
                    map_info.append([info[0], map_pos, map_len])
    return map_info

def cal_ref_cover(map_info, norm_len=1500):
    ref_cover_info = np.zeros(norm_len)
    for i in map_info:
        norm_s = int(i[1])
        norm_e = int(i[2])
#           norm_s = int(np.floor(int(i[1])/ref_fa_len*norm_len))
#           norm_e = int(np.floor(int(i[2])/ref_fa_len*norm_len))
        ref_cover_info[norm_s:norm_s+norm_e]+=1
    return ref_cover_info

def get_map_info_by_ref(map_sam_out):
    ref_dict = {}
    for i in open(map_sam_out).readlines():
        if re.match('@', i):
            next
        else:
            info = i.split('\t')
            if len(info)<5:
                next
            else:
                ref_dict[info[2]] = ref_dict.get(info[2], [])
                ref_dict[info[2]].append(i)
    return ref_dict            


def plot_contig_cover_box(refs_cover_info, ax=None):
    ax = ax or plt.gca()
    ax = sns.boxplot(x='variable', y='value', data=pd.DataFrame(refs_cover_info).melt(), color='w', width=0.5, linewidth=0.3, fliersize=0.2)
    # for patch in ax.artists:
    #     r, g, b, a = patch.get_facecolor()
    #     patch.set_facecolor((0, 0, 0, .3))
    # ax.plot(np.array(refs_cover_info).mean(axis=0), color='black')
    x_pos = range(0, pd.DataFrame(refs_cover_info).shape[1]+1, 100)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_pos)
    ax.set_xlabel('16S DNA Position')
    ax.set_ylabel('Cover Depth')
    return ax

def plot_contig_cover_line(refs_cover_info, ax=None):
    ax = ax or plt.gca()
    ax = sns.lineplot(x='variable', y='value', data=pd.DataFrame(refs_cover_info).melt())
    x_pos = range(0, pd.DataFrame(refs_cover_info).shape[1]+1, 100)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_pos)
    ax.set_xlabel('16S DNA Position')
    ax.set_ylabel('Cover Depth')
#     ax.set_title('Contigs Coverage')
    return ax


def plotCover(A1_file, bowtie_db, bowtie_path, out_dir='./', out_png='', threads=2):
    out_file = os.path.join(out_dir, A1_file+'.100000.bowtie.sam')
    log_file = out_png + '_bowtie.log'
#     out_png = os.path.join(out_dir, A1_file+'.100000.cover')
    submit_bowtie(A1_file, out_file, bowtie_db, log_file, threads, bowtie_path)
    ref_dict = get_map_info_by_ref(out_file)
    ref_cover_info = []
    for i in ref_dict.keys():
        map_info_tmp = get_map_info(ref_dict[i])
        ref_cover_info_tmp = cal_ref_cover(map_info_tmp)
        ref_cover_info.append(ref_cover_info_tmp)
    
    fig, ax = plt.subplots()
    plot_contig_cover_line(ref_cover_info)
    plt.xticks(rotation=60)
    fig.savefig(out_png+'_line.png')
    fig, ax = plt.subplots()
    plot_contig_cover_box(ref_cover_info)
    plt.xticks(rotation=60)
    fig.savefig(out_png+'_box.png')
