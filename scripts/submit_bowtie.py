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
import os

def submit_bowtie(in_file,out_file,db,log_file,threads,bowtie_path):
    '''
    bowtie alignment
    '''
    str_join=' '
#     bowtie2 -u 10000 -p 8 -x index -1 reads1.fq -2 reads2.fq -S out.sam
    cmd=str_join.join([bowtie_path,
                       '-u', str(100000),
                       '-p', str(threads),
                       '-x',db,
                       '-q',in_file,
                       '-S',out_file,
                       '1>',log_file,
                       '2>',log_file])
    logging.info(' %s'%cmd)
    status=os.system(cmd)
    if status==0:
        logging.info(' done')
    else:
        logging.info(' exit')

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description=\
#                                      'this is for alignment of fasta \
#                                      using bowtie2')
#     parser.add_argument('-i','--in_fasta',help='input fasta file',required=True)
#     parser.add_argument('-o','--out_sam',help='output sam file',required=True)
#     parser.add_argument('-l','--log',help='output log file. default=out.log',\
#                         default='')
#     parser.add_argument('-d', '--db', help='bowtie db', default='/Bioinfo/Database/16S_BSI/GM_BSI_v3')
#     parser.add_argument('-t','--threads',help='threads for bowtie. default=1',default=1,type=int)
#     parser.add_argument('-path','--bowtie2_path',help='bowtie2 absolute path.',
#                         default='/root/anaconda3/bin/bowtie2')
#     args = parser.parse_args()
#     logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level=logging.INFO)
#     log_check = args.log if args.log else args.out_sam + '.log'
#     submit_bowtie(args.in_fasta,args.out_sam,args.db,log_check,args.threads,args.bowtie2_path)  
