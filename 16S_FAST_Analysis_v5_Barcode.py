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
import sys
import os
sys.path.append(os.path.join(sys.path[0], './scripts/'))
import pandas as pd
import logging
import argparse
import psutil

script_path = sys.path[0]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Barcode Analysis')
    parser.add_argument('-c', '--config_file', help='config file including program path', 
                        default=os.path.join(script_path,'config_path_linux.txt'))
    parser.add_argument('-f', '--A1_file', help='Assemble data R1 fastq', required=True)
    parser.add_argument('-f2', '--A2_file', help='Assemble data R2 fastq', default=None)
    parser.add_argument('-u', '--uID_file', help='R1 reads ID to umiID file, with columns [umiID|AreadsID|Sample]', required=True)
    parser.add_argument('-FB', '--Fbarcode', help='F barcode, eg. F1,F2,F3...', required=True)
    parser.add_argument('-RB', '--Rbarcode', help='R barcode, eg. R1,R2,R3...', required=True)
    parser.add_argument('-n', '--Name', help='Sample name in uID file', required=True)
    parser.add_argument('-o', '--out_dir', help='result out dir', required=True)
    
    args = parser.parse_args()
    
    logging.basicConfig(filename=os.path.join(args.out_dir, 'log'), format='%(asctime)s %(levelname)s:%(message)s',level='INFO')
#     logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')
    contig_path = pd.read_csv(args.config_file, sep=':', header=None).set_index(0).to_dict()[1]
    contig_path['threads'] = psutil.cpu_count()
    logging.info(contig_path)
    
    from BarcodeAna import BarcodeAna
    BarcodeAna(args.A1_file, args.uID_file, 
               A2=args.A2_file,
               Fbarcode=args.Fbarcode, 
               Rbarcode=args.Rbarcode, 
               barcode_name=args.Name,
               out_dir=args.out_dir,
               **contig_path
              )
