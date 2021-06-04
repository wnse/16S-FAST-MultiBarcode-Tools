
python /Bioinfo/16S-FAST-Tools-2021/16S_FAST_Analysis_v5.py \
-tag NP_1 \
-ar1 /mnt/data/16S-FAST_data/20210128/rawdata/210125_A00262_0590_AHVYVTDSXY/2-P_L2_Y0000716Y0000712.R1.fastq.gz \
-ar2 /mnt/data/16S-FAST_data/20210128/rawdata/210125_A00262_0590_AHVYVTDSXY/2-P_L2_Y0000716Y0000712.R2.fastq.gz \
-lr1 /mnt/data/16S-FAST_data/20210128/rawdata/210125_A00262_0590_AHVYVTDSXY/2-L_L2_Y0000680Y0000712.R1.fastq.gz \
-lr2 /mnt/data/16S-FAST_data/20210128/rawdata/210125_A00262_0590_AHVYVTDSXY/2-L_L2_Y0000680Y0000712.R2.fastq.gz \
-FB F1 F1 F1 F1 F1 F2 F2 F2 \
-RB R2 R3 R4 R5 R8 R4 R6 R7 \
-N F1R2 F1R3 F1R4 F1R5 F1R8 F2R4 F2R6 F2R7 \
-o /mnt/data/work/Project/YF/20210128_multisample/20210128/NP_2 \
-remoteA1 /mnt/data/16S_FAST_data/210125_A00262_0590_AHVYVTDSXY/1-P_L2_Y0000085Y0000712.R1.fastq.gz \
-remark 20210128 \
> /mnt/data/work/Project/YF/20210128_multisample/20210128/NP_2/log \

shutdown -h 10





python /Bioinfo/16S-FAST-Tools-2021/16S_FAST_Analysis_v5.py \
-tag DS0_1 \
-ar1 /root/data_size/rawdata/subseq/s0.r1.fq.gz \
-ar2 /root/data_size/rawdata/subseq/s0.r2.fq.gz \
-lr1 /root/data_size/rawdata/mnt/data/16S_FAST_data/201128_A00869_0348_AHL3JMDSXY/1-L_L4_Y0000680Y0000712.R1.fastq.gz \
-lr2 /root/data_size/rawdata/mnt/data/16S_FAST_data/201128_A00869_0348_AHL3JMDSXY/1-L_L4_Y0000680Y0000712.R2.fastq.gz \
-FB F1 F1 F2 \
-RB R1 R2 R1 \
-N F1R1 F1R2 F2R1 \
-o /root/data_size/resutl/DS0_1/ \
-remoteA1 /mnt/data/16S_FAST_V5_UMI/test_rawdata/s0.r1.fq.gz \
-remark 20210309 
# > /root/data_size/resutl/DS0_1/log \

shutdown -h 10










