#!/usr/bin/env python

import os

def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))


def gbk2faa(input_path, ouput_prefix):
    from Bio import SeqIO
    
    print("readinmg records")
    records = SeqIO.parse(input_path, "genbank")
    
    o = open(ouput_prefix + '.faa', 'w')
    print(ouput_prefix + '.faa')
    for record in records:
        print(record.id)
        for feature in record.features:
            if feature.type == 'CDS' and 'translation' in feature.qualifiers and 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers["locus_tag"][0]
                seq = str(feature.qualifiers["translation"][0])
                o.write(f">{locus_tag}\n")
                for chunk in chunkstring(seq, 60):
                    o.write(f'{chunk}\n')
    o.close()
    


def locus2gene_and_product(gbk_file):
    
    from Bio import SeqIO

    locus2annotation = {}
    locus2length = {}
    records = SeqIO.parse(gbk_file, "genbank")
    for record in records:
        print(record.id)
        for feature in record.features:
            if feature.type == 'CDS' and 'translation' in feature.qualifiers and 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers["locus_tag"][0]
                try:
                    gene = feature.qualifiers["gene"][0]
                except KeyError:
                    gene = '-'
                try:
                    product = feature.qualifiers["product"][0]
                except KeyError:
                    product = '-'

                locus2annotation[locus_tag] = f'{gene}: {product}'
                locus2length[locus_tag] = len(feature.qualifiers["translation"][0])
    return locus2annotation,locus2length


def locus2best_hit(blast_tab_file):
    import pandas 
    
    header = ["qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore"]

    ssearch_df = pandas.read_csv(blast_tab_file, names=header, sep="\t")
    
    locus2best_hit = {}
    for n,row in ssearch_df.iterrows():
        if row["qseqid"] not in locus2best_hit:
            locus2best_hit[row["qseqid"]] = row        

    return locus2best_hit
    

def RBBH_table(gbk_file_list,
               ssearch_result_list):
    
    from itertools import combinations 
    
       
    locus2annot = {}
    locus2len = {}
    for gbk in gbk_file_list:
        print(gbk)
        tmp_locus2annot, tmp_locus2len = locus2gene_and_product(gbk)
        locus2annot.update(tmp_locus2annot)
        locus2len.update(tmp_locus2len)
    
    prefix_lst = [os.path.splitext(i)[0] for i in gbk_file_list]
    
    comb = combinations(prefix_lst, 2)
    
    for genome1, genome2 in comb:
        locus2best_hit_a_vs_b = locus2best_hit(f"{genome1}_vs_{genome2}.tab")
        locus2best_hit_b_vs_a = locus2best_hit(f"{genome2}_vs_{genome1}.tab")
        print(list(locus2best_hit_a_vs_b.keys())[0:10])
        print(list(locus2best_hit_b_vs_a.keys())[0:10])
        
        o = open(f"{genome1}_{genome2}_RBBH.tab", "w")
        header = ["locus_a",
                  "locus_b",
                  "length_a",
                  "length_b",
                  "coverage_a",
                  "coverage_b",
                  "identity_a_vs_b",
                  "identity_b_vs_a",
                  "annot_a",
                  "annot_b"]
        o.write("\t".join(header) + '\n')
        for locus_a in locus2best_hit_a_vs_b:
            locus_b = locus2best_hit_a_vs_b[locus_a]["sseqid"]
            
            try: 
                reverse_search = locus2best_hit_b_vs_a[locus_b]["sseqid"]
            except KeyError:
                # not match reverse (can happen if one sequence is much shorter for instance => only significant onw way)
                continue
            if reverse_search == locus_a:
                identity_a_vs_b = locus2best_hit_a_vs_b[locus_a]["pident"]
                identity_b_vs_a = locus2best_hit_b_vs_a[locus_b]["pident"]
                len_locus_a = locus2len[locus_a]
                len_locus_b = locus2len[locus_b]
                len_alignment_a = (locus2best_hit_a_vs_b[locus_a]["qend"] - locus2best_hit_a_vs_b[locus_a]["qstart"]) + 1
                len_alignment_b = (locus2best_hit_a_vs_b[locus_a]["send"] - locus2best_hit_a_vs_b[locus_a]["sstart"]) + 1
                coverage_a = round((float(len_locus_a)/len_alignment_a) * 100, 2)
                coverage_b = round((float(len_locus_b)/len_alignment_b) * 100, 2)
                o.write(f'{locus_a}\t{locus_b}\t{len_locus_a}\t{len_locus_b}\t{coverage_a}\t{coverage_b}\t{identity_a_vs_b}\t{identity_b_vs_a}\t{locus2annot[locus_a]}\t{locus2annot[locus_b]}\n')
    
    