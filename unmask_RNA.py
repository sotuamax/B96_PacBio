#!/usr/bin/env python

"""

This script was developed as a part of gene annotation pipeline. 
For soft/hard-masked genome, its gene spaces can sometimes be masked and affect following gene prediction process. To alleviate over-masking problem, the program takes soft-masked genome and RNA-alignment information to unmask regions which can include transcripts. 

For any questions, please directly report to the developer: 
Meiyuan Ji (meiyuan.ji@utah.edu)
Last change: May/06/2022

"""

import time 
import pysam
from Bio import SeqIO
import argparse
import pandas as pd 
import re 
import os 
import subprocess
from mpi4py import MPI
import sys
import shutil 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description=" \n\
        Usage: mpiexec -n 10 python unmask_RNA.py -fasta <ref> -bam <bam> -O <prefix>")
    parser.add_argument("-fasta", "--fasta", help = "reference fasta file using soft-masked. ")
    parser.add_argument("-bam", "--bam", help = "bam file of RNA-seq alignment. ")
    parser.add_argument("-step", "--step", default = 5000, type = int, help = "step size for searching soft-masked region (default: 5000). ")
    parser.add_argument("-coverage", "--coverage", default = 3, type = int, help = "RNA-seq coverage to count as gene-space (default: > 3). ")
    parser.add_argument("-MQ", "--mapping_quality", default = 60, type = int, help = "mapping quality for read to be count as effective (default: >= 60). ")
    parser.add_argument("-O", "--output", help = "output file prefix name. ")
    args = parser.parse_args()
    return args

def parse_fasta(fa):
    seq = SeqIO.parse(fa, "fasta")
    seq_dict = SeqIO.to_dict(seq)
    return seq_dict

def fasta_mask_report(fasta_dict, seq_ids):
    fasta_mask_len = dict()
    for seq_id in seq_ids:
        sstr = str(fasta_dict[seq_id].seq)
        seq_len = len(sstr)
        mask_len = len(re.findall("[atcg]", sstr))
        fasta_mask_len[seq_id] = (mask_len, seq_len)
    return fasta_mask_len

def bam_fasta_match(args, seq_ids):
    bam = args.bam
    bam_handle = pysam.AlignmentFile(bam, "rb")
    bam_seqs = sorted([i.split("SN:")[-1].split("\t")[0] for i in str(bam_handle.header).split("\n") if i.startswith("@SQ")])
    if bam_seqs == seq_ids:
        return True
    else:
        return False

def pos_pileup(bam_handle, chr, pos, coverage, quality):
    count_read = list()
    for read_column in bam_handle.pileup(chr, pos, pos+1, min_mapping_quality = quality):
        if read_column.reference_pos == pos:
            for bam_read in read_column.pileups:
                if bam_read.query_position != None:
                    read_name = bam_read.alignment.query_name
                    count_read.append(read_name)
                    read_dep = len(count_read)
                    if read_dep >= coverage:
                        break
    return len(count_read)

def parse_bam(args, seq_dict, seq_ids, out):
    bam = args.bam
    step = args.step
    coverage = args.coverage
    quality = args.mapping_quality
    bam_handle = pysam.AlignmentFile(bam, "rb")
    fw = open(out + ".tmp", "w") 
    fw.write("chromosome\tposition\n")
    ###
    for seq_id in seq_ids:
        seq_str = str(seq_dict[seq_id].seq)
        seq_start = 0
        seq_index = 0
        int_list = list()
        while seq_index != None:
            pattern_search = re.search("[actg]", seq_str[seq_start:])
            if pattern_search == None:
                seq_index = None
            else:
                seq_index = pattern_search.start()
                bam_search_pos = seq_index + seq_start # 0-based
                count_num = pos_pileup(bam_handle, seq_id, bam_search_pos, coverage, quality)
                if count_num >= coverage:
                    fw.write(f"{seq_id}\t{bam_search_pos}\n")
                bam_search_pos += step
                bam_region_int = bam_search_pos//1000000
                if not bam_region_int in int_list:
                    int_list.append(bam_region_int)
                    print(f"### Processing {seq_id} on {bam_region_int} million bp ......")
                seq_start = bam_search_pos
    fw.close()
    ###
    rank_df = pd.read_table(out + ".tmp", sep = "\t", header = 0)
    return rank_df

def fasta_unmask(args, seq_dict, seq_ids, out):
    step = args.step
    unmask_df = pd.read_table(out + ".tmp", sep = "\t", header = 0)
    unmask_chr = sorted(set(unmask_df["chromosome"]))
    fa_fw = open(out + ".tmp.fasta", "w")
    unchange_seqs = set(seq_ids).difference(set(unmask_chr))
    for chr in unmask_chr:
        i = 0
        seq_str = str(seq_dict[chr].seq)
        unmask_pos_list = list(unmask_df[unmask_df["chromosome"] == chr]["position"])
        for p in unmask_pos_list:
            if i == 0:
                seq_str = seq_str[:p] + seq_str[p:p+step].upper() + seq_str[p+step:]
            else:
                seq_str = seq_str[:p-step] + seq_str[p-step:p+step].upper() + seq_str[p+step:]
            i += 1
        fa_fw.write(f">{chr}\n{seq_str}\n")
    for chr in unchange_seqs:
        seq_str = str(seq_dict[chr].seq)
        fa_fw.write(f">{chr}\n{seq_str}\n")
    fa_fw.close()

def hard_mask(new_seq_dict, out):
    fw = open(out + ".hard.fasta", "w")
    for s in new_seq_dict:
        seq_str = str(new_seq_dict[s].seq)
        seq_str = seq_str.replace("a", "N")
        seq_str = seq_str.replace("t", "N")
        seq_str = seq_str.replace("c", "N")
        seq_str = seq_str.replace("g", "N")
        fw.write(f">{s}\n{seq_str}\n")
    fw.close()

def endtime(start):
    end = time.time()
    t = end-start
    if t < 60:
        print('{:.2f} seconds elapsed'.format(t))
    elif t < 3600:
        print('{:.2f} minutes elapsed'.format(t/60))
    else:
        print('{:.2f} hours elapsed'.format(t/3600))
    return end

def rm_temp(args, size):
    output = args.output
    if output + ".txt" in os.listdir():
        for s in range(size):
            os.remove("rank_" + str(s) + ".tmp")
            os.remove("rank_" + str(s) + ".tmp.fasta")

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    ### 
    args = args_parser()
    out = args.output
    seq_dict = parse_fasta(args.fasta)
    seq_ids = sorted([seq_dict[entry].id for entry in seq_dict]) # seq_ids never change in the process
    if bam_fasta_match(args, seq_ids) == False:
        if rank == 0:
            print("FASTA and BAM file not match! ")
            print("Please provide correct FASTA and BAM file! ")
        sys.exit()
    ###
    if rank == 0:
        if not os.path.isdir("tmp"):
            os.makedirs("tmp")
        print("### Build fasta file directory ......")
        start = time.time()
        worker_tasks = {w:[] for w in range(size)}
        w_idx = 0
        for s in seq_ids:
            worker_tasks[w_idx].append(s)
            w_idx = (w_idx + 1) % size
        print("split work ......")
    else:
        worker_tasks = None
        start = None
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    start = comm.bcast(start, root = 0)
    ###
    s = worker_tasks[rank] # s is seq_ids in a list
    mask_seq = fasta_mask_report(seq_dict, list(s))
    mask_seqs_list = comm.gather(mask_seq, root = 0)
    if rank == 0:
        old_seq_mask = {i:j for d in mask_seqs_list for i,j in d.items()}
    rank_df = parse_bam(args, seq_dict, list(s), os.path.join("tmp", "rank_"+str(rank)))
    endtime(start)
    print(f"### Finish process for bam-searching on rank {rank} ...... ")
    fasta_unmask(args, seq_dict, list(s), os.path.join("tmp", "rank_" + str(rank)))
    endtime(start)
    print(f"### Finish process for unmask-correction on rank {rank} ...... ")
    df = comm.gather(rank_df, root = 0)
    comm.barrier()
    if rank == 0:
        df_all = pd.concat(df, axis = 0)
        df_all.to_csv(out + ".txt", sep = "\t", index = False)
        fa_fs = sorted([os.path.join("tmp", f) for f in os.listdir("tmp") if f.endswith(".tmp.fasta")])
        fa_tmps = " ".join(fa_fs)
        command = f"cat {fa_tmps} > {out}.soft.fasta"
        subprocess.call(command, shell = True)
        new_seq_dict = parse_fasta(out + ".soft.fasta")
        # seq_ids_sorted = sorted([new_seq_dict[entry].id for entry in new_seq_dict])
        print("### Write new soft-masked genome ...... ")
        new_fw_fa = open(out + ".soft.fasta", "w")
        for s in seq_ids:
            new_seq_str = str(new_seq_dict[s].seq)
            new_fw_fa.write(f">{s}\n{new_seq_str}\n")
        new_fw_fa.close()
    comm.barrier()
    new_seq_dict = parse_fasta(out + ".soft.fasta")
    s = worker_tasks[rank]
    mask_seq = fasta_mask_report(new_seq_dict, s)
    mask_seqs_list = comm.gather(mask_seq, root = 0)
    if rank == 0:
        new_seq_mask = {i:j for d in mask_seqs_list for i,j in d.items()}
        print("### Write hard-masked genome ...... ")
        hard_mask(new_seq_dict, out)
        fw_report = open(out + ".log", "w")
        fw_report.write("chromosome\told_mask_length\tnew_mask_length\ttotal\n")
        for sid in seq_ids:
            omask_len = old_seq_mask[sid][0]
            olen = old_seq_mask[sid][1]
            nmask_len = new_seq_mask[sid][0]
            nlen = new_seq_mask[sid][1]
            if olen == nlen:
                fw_report.write(f"{sid}\t{omask_len}\t{nmask_len}\t{olen}\n")
            else:
                print(f"{sid} sequence changes in length. ")
        fw_report.close()
        if out + ".txt" in os.listdir() and out + ".soft.fasta" in os.listdir() and out + ".hard.fasta" in os.listdir():
            shutil.rmtree("tmp")
            endtime(start)
            print("### Successfully finish all! ")


##############
### Run it ###
##############

if __name__ == "__main__":
    main()

