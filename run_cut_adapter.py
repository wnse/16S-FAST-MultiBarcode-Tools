import os
import logging

logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')

def submit_cutadapt(in_file, out_file, log_file, adapter, g_a, path, error=0.1, threads=1, info=0):
    '''去除接头序列
    参数：
        in_file: 输入fastq文件
        out_file: 输出fastq文件
        log_file: 输出日志
        adapter: 接头序列
        g_a: 去除方式（g/a）
        path: cutadapt路径
        threads: 进程数
        info: 是否输出每条序列的接头信息
    返回：
        无
    '''
    str_join = ' '
    info_file = out_file + '.cutadapt.info.file'
    cutadapt = path
    if info:
        cmd = str_join.join(
            [cutadapt, '-%s' % g_a, adapter, '-o', out_file, in_file, '-e', str(error), '--info-file', info_file, '>', log_file])
    else:
        cmd = str_join.join(
            [cutadapt, '-%s' % g_a, adapter, '-o', out_file, in_file, '-e', str(error), '-j', str(threads), '>', log_file])
    logging.info(' %s' % cmd)
    # print(cmd)
    status = os.system(cmd)
    if status == 0:
        logging.info(' done')
    else:
        logging.info(' exit')

def mkdir(dir_name):
    '''建立目录'''
    if not os.path.isdir(dir_name):
        os.makedirs(dir_name)
        logging.info(' mkdir ' + dir_name)

cutadapt = '/root/anaconda3/bin/cutadapt'
fa_dir = '/mnt/data/work/Project/YF/20201202_multisample/20201231/adapter_fa/'
F_a_rc = os.path.join(fa_dir, 'F_a_rc.fa')
F_l_rc = os.path.join(fa_dir, 'F_l_rc.fa')
R_a = os.path.join(fa_dir, 'R_a.fa')
R_l = os.path.join(fa_dir, 'R_l.fa')

R_a_rc = os.path.join(fa_dir, 'R_a_rc.fa')
R_l_rc = os.path.join(fa_dir, 'R_l_rc.fa')
F_a = os.path.join(fa_dir, 'F_a.fa')
F_l = os.path.join(fa_dir, 'F_l.fa')


d = './cut_adapter_out_e'
for i in open('rawdata.list', 'rt').readlines():
    temp = i.strip().split('\t')
    File_Tag = temp[0]
    L1 = temp[1]
    L2 = temp[2]
    ana_dir = os.path.join(d, File_Tag)
    threads = 2
    mkdir(ana_dir)

    L1_cut_adapter = os.path.join(ana_dir,File_Tag+'_cut_adapter_10.fq.gz')
    L1_cut_adapter_log = os.path.join(ana_dir,File_Tag+'_cut_adapter_10.log')
    L1_cut_adapter_1 = os.path.join(ana_dir,File_Tag+'_cut_adapter_11.fq.gz')
    L1_cut_adapter_log_1 = os.path.join(ana_dir,File_Tag+'_cut_adapter_11.log')

    umi1_fq = os.path.join(ana_dir,File_Tag+'_cut_primer_10.fq')
    umi1_fq_log = os.path.join(ana_dir,File_Tag+'_cut_primer_10.log')
    umi1_fq_1 = os.path.join(ana_dir,File_Tag+'_cut_primer_11.fq')
    umi1_fq_log_1 = os.path.join(ana_dir,File_Tag+'_cut_primer_11.log')

    L2_cut_adapter = os.path.join(ana_dir,File_Tag+'_cut_adapter_20.fq.gz')
    L2_cut_adapter_log = os.path.join(ana_dir,File_Tag+'_cut_adapter_20.log')
    L2_cut_adapter_1 = os.path.join(ana_dir,File_Tag+'_cut_adapter_21.fq.gz')
    L2_cut_adapter_log_1 = os.path.join(ana_dir,File_Tag+'_cut_adapter_21.log')

    umi2_fq = os.path.join(ana_dir,File_Tag+'_cut_primer_20.fq')
    umi2_fq_log = os.path.join(ana_dir,File_Tag+'_cut_primer_20.log')
    umi2_fq_1 = os.path.join(ana_dir,File_Tag+'_cut_primer_21.fq')
    umi2_fq_log_1 = os.path.join(ana_dir,File_Tag+'_cut_primer_21.log')


    submit_cutadapt(L1,L1_cut_adapter,L1_cut_adapter_log,
                    'file:'+F_a_rc ,'g',cutadapt, threads=threads)
    submit_cutadapt(L1_cut_adapter,umi1_fq,umi1_fq_log,
                    'file:'+F_l_rc,'a',cutadapt, info=1)
    submit_cutadapt(L1, L1_cut_adapter_1, L1_cut_adapter_log_1,
                    'file:'+R_a, 'a', cutadapt, threads=threads)
    submit_cutadapt(L1_cut_adapter_1, umi1_fq_1, umi1_fq_log_1,
                    'file:'+R_l, 'g', cutadapt, info=1)

    submit_cutadapt(L2,L2_cut_adapter,L2_cut_adapter_log,
                    'file:'+R_a_rc ,'g',cutadapt, threads=threads)
    submit_cutadapt(L2_cut_adapter,umi2_fq,umi2_fq_log,
                    'file:'+R_l_rc,'a',cutadapt, info=1)
    submit_cutadapt(L2, L2_cut_adapter_1, L2_cut_adapter_log_1,
                    'file:'+F_a, 'a', cutadapt, threads=threads)
    submit_cutadapt(L2_cut_adapter_1, umi2_fq_1, umi2_fq_log_1,
                    'file:'+F_l,'g', cutadapt, info=1)

