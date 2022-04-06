#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import logging
import os


def submit_cutadapt(in_file, out_file, log_file, 
                    adapter, g_a, 
                    path='/Bio/bin/cutadapt', 
                    error=0, threads=1, info=0):
    '''trimming adapter of fastq using cutadapt
    Parameters
    ----------
        in_file: input fastq.
        out_file: output fastq.
        log_file: output log file.
        adapter: adapter sequences.
        g_a: g(5') or a(3') for trimming.
        path: cutadapt absolute path.
        error: adapter sequences error.
        threads: threads for computing.
        info: output info file named out_fastq.cutadapt.info.file, if set -i,threads not useful.
        
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
    status = os.system(cmd)
    return status
