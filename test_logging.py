import os
import time
import sys
sys.path.append(os.path.join(sys.path[0], './scripts/'))
from mkdir import mkdir
from submit_Trimmomatic_SE import submit_Trimmomatic_SE
from plotCover import plotCover
from submit_cutadapt import submit_cutadapt
from cutLadapter import cutLadapter
from staLumi import get_df_merge
from staLumi import get_all_sample_umi
from pairLumi import dropRepeatUid
from pairLumi import pairLumi
from pairLumi import plot_Lreads
from cutAadapter import cutAadapter
from staAumi import staAumi
from getAseqByBarcode import reads2uID_total
from plotAumiReads import plotAreads
from BarcodeAna import BarcodeAna
from qsctl_sync import qsctl_sync
from submit_FAST_Barcode import submit_FAST_Barcode

import psutil
import shutil
import argparse
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import time

def test_def():
	logging.info(' remove dir ')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Barcode Analysis')
	parser.add_argument('-c', '--config_file', help='config file including program path', default='config_path_linux.txt')
	args = parser.parse_args()
	logging.basicConfig(filename=os.path.join('./', 'log'), format='%(asctime)s %(levelname)s:%(message)s',level='INFO')
	logging.info('test info')
	test_def()
	# time.sleep(10)
	test_def()
	
