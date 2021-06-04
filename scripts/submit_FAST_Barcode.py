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
import json
import os
import argparse
import time
import logging
                
def submit_FAST_Barcode(ar1, A_uID_file,
                        FB, RB,
                        remark, lib_name, sam_name,
                        result_path,
                        *args
                       ):
    result_path = os.path.join(result_path,str(remark),str(lib_name),str(sam_name))
    cmd_list = ['/root/anaconda3/bin/python3.7', 
                '/Bioinfo/16S-FAST-Tools-2021/16S_FAST_Analysis_v5_Barcode.py',
                '-f ${ar1}',
                '-u ${A_uID_file}',
                '-FB', FB,
                '-RB', RB,
                '-o', '${outdir}',
                '-n', sam_name,
               ]
    if args:
        cmd_list.extend(list(args))
    cmd=' '.join(cmd_list)
    # print(cmd)
    config={
                'application':
                {
                    'cmd':cmd,
                    'id':2,
                    'name':'16S-FAST-Tools-2021',
                },
                'baseData':
                {
                    'remark':remark
                },
                'datas':
                [
                    {
                        'objectKey':ar1,'paramKey':'ar1'
                    },
                    {
                        'objectKey':A_uID_file,'paramKey':'A_uID_file'
                    },
                ],
                'name':sam_name,
                'result':
                {
                    'path':result_path,
                    'paramKey':'outdir'
                },
                'type':'aln'    
    }
    config_json=(json.dumps(config))
    commond="curl -X POST https://bc.dev.germountx.com/api/tasks -H 'Content-Type: application/json' -d " + \
            "'" + config_json + "'"
    out=os.popen(commond)
    # print(ar1, ar2, lr1, lr2, name, result_path)
    logging.info('{}\t{}\t{}'.format(out.read(), sam_name, result_path))
    # logging.info('{} {}'.format(name,out.read()))


# %%
# A_uID_file = '/mnt/data/16S_FAST_V5_UMI/test/LibaryName/test_A_uID_dropRepeat.csv'
# A1_file = '/mnt/data/16S_FAST_V5_UMI/test/1-P_R1_test.fastq.gz'
# FB = 'F1'
# RB = 'R5'
# result_path = '/mnt/data/yangkai/16S_FAST/V5/'
# remark = '202103080653'
# library_name = 'Library'
# sample_name = 'sample'
# submit_FAST_Barcode(A1_file, A_uID_file, FB, RB, remark, library_name, sample_name, result_path)

# %%
# if __name__ == '__main__':
#     logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level=logging.INFO) 
#     date = time.strftime("%Y%m%d%H%M%S", time.localtime()) 
#     __version__ = "0.1_20190918"
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--version', action='version',version='%(prog)s {version}'.format(version=__version__))
#     parser.add_argument('-ar1','--Assemble_Read_1',help='Assemble Read 1 fastq data absolute path in bc-input',required=True)
#     parser.add_argument('-uID','--A_uID_file',help='AuIDfile in bc-input',required=True)
#     parser.add_argument('-n','--name',help='sample name.',required=True)
#     parser.add_argument('-remark','--remark',help='batch name. default:Current Time',default=date)
#     parser.add_argument('-o','--outdir',help='result dir in bc-output. default:/mnt/data/yangkai/16S_FAST/remark/name',\
#                         default='/mnt/data/yangkai/16S_FAST/V5')
#     parser.add_argument('-add_para', '--add_parameter', default='', help='additional parameter', type=str)
#     args = parser.parse_args() 
#     ar1 = args.Assemble_Read_1
#     A_uID = args.A_uID_file
#     remark = args.remark
#     name = args.name
#     outdir = args.outdir
#     add_para = args.add_parameter
#     submit_16S_FAST(ar1, A_uID, FB, RB, remark, name, outdir, add_para)
