#!python3

import pandas as pd
import numpy as np

# ---- Useful Functions ---- #


def obtain_vcf_manifest(chrom, pop):
    """Obtain a VCF file using a chromosome and population notation."""
    vcf_file = manifest[
        (manifest.chrom == chrom) & (manifest.population == pop)
    ].vcf_file.values[0]
    return vcf_file


def obtain_all_chroms(pop):
    """Obtain all chromosomes for a specific population."""
    return np.unique(manifest[manifest.population == pop].chrom.values)


# --- Rules --------------------------- #
rule conv_vcf2smc_1pop:
    """Converts vcf2smc files for multiple focal individuals at once."""
    input:
        vcf=lambda wildcards: obtain_vcf_manifest(wildcards.chrom, wildcards.focal_pop),
        panel=lambda wildcards: config["datasets"][wildcards.focal_pop]["popfile"],
        chrom_lengths=lambda wildcards: config["datasets"][wildcards.focal_pop][
            "contig_lengths"
        ],
    output:
        smc_out=temp(
            "results/smcpp_format/{focal_pop}/{focal_pop}.{chrom}.{focal}.smcpp.gz"
        ),
    params:
        exclusion_mask=lambda wildcards: f'-m {config["datasets"][wildcards.focal_pop]["mask"]}'
        if config["datasets"][wildcards.focal_pop]["mask"] != ""
        else "",
    resources:
        time="2:00:00",
        mem_mb="2G",
    conda:
        "../envs/smcpp.yaml"
    shell:
        """
        contig_len=$(awk \'$1 == \"{wildcards.chrom}\" {{print $2}}\' {input.chrom_lengths})
        pop_str=$(awk \'{{print $1}}\' {input.panel} | paste -s -d, - | sed -e \'s/^/{wildcards.focal_pop}:/\')
        smc++ vcf2smc {params.exclusion_mask} --length $contig_len -d {wildcards.focal} {wildcards.focal} {input.vcf} {output.smc_out} {wildcards.chrom} $pop_str
        """


rule smcpp_estimate_single_pop:
    """Run SMC++ estimation of trajectory."""
    input:
        smcpp_files=lambda wildcards: expand(
            "results/smcpp_format/{{focal_pop}}/{{focal_pop}}.{chrom}.{focal}.smcpp.gz",
            chrom=obtain_all_chroms(wildcards.focal_pop),
            focal=config["datasets"][wildcards.focal_pop]["focal_indiv"],
        ),
    output:
        smc_out="results/smcpp_output_mult_final/{focal_pop}_t1_{t1}_knots_{knots}_filt/model.final.json",
    params:
        t1=lambda wildcards: f"--timepoints {wildcards.t1} 30000"
        if wildcards.t1 != "None"
        else "",
        knots=lambda wildcards: f"--knots {wildcards.knots}"
        if wildcards.knots != "None"
        else "",
        mu=config["mu"],
    conda:
        "../envs/smcpp.yaml"
    threads: 8
    resources:
        time="6:00:00",
        mem_mb="8G",
    shell:
        "smc++ estimate -o results/smcpp_output_mult_final/{wildcards.focal_pop}_t1_{wildcards.t1}_knots_{wildcards.knots}_filt/  {params.t1} {params.knots} --cores {threads} -v {params.mu} {input.smcpp_files}"


rule gen_smcpp_csvs_single:
    """Generate CSVs of SMC++ trajectories."""
    input:
        smc_out=rules.smcpp_estimate_single_pop.output.smc_out,
    output:
        smc_csv="results/smcpp_output_mult_final/{focal_pop}_t1_{t1}_knots_{knots}_filt/{focal_pop}_t1_{t1}_knots_{knots}_filt.csv",
        smc_pdf=temp(
            "results/smcpp_output_mult_final/{focal_pop}_t1_{t1}_knots_{knots}_filt/{focal_pop}_t1_{t1}_knots_{knots}_filt.pdf"
        ),
    resources:
        time="0:30:00",
        mem_mb="1G",
    conda:
        "../envs/smcpp.yaml"
    shell:
        "smc++ plot {output.smc_pdf} {input.smc_out} --csv"


# # ---- Estimating Split Times ---- #
# rule conv_vcf2smc_2pop:
# input:
# vcf1= VCF_PATH_MAIN + '{pop1}.chr{CHROM}.smcpp.primary.vcf.gz',
# vcf2= VCF_PATH_MAIN + '{pop2}.chr{CHROM}.smcpp.primary.vcf.gz',
# panel1 = VCF_PATH_MAIN + 'aux.indlist.sampleid.primary.smcpp.{pop1}.txt',
# panel2 = VCF_PATH_MAIN + 'aux.indlist.sampleid.primary.smcpp.{pop2}.txt',
# inclusion_mask = DATA_PATH + 'mappability_masks/hs37d5_chr{CHROM}.mask.bed.gz',
# contig_length = DATA_PATH + 'contig_lengths.txt'
# output:
# temp_vcf_merged = temp(VCF_PATH_MAIN + 'data/smcpp_input/{pop1}_{pop2}/{pop1}.{pop2}.chr{CHROM, \d+}.{pop1}.split.vcf.gz'),
# temp_vcf_merged_index = temp(VCF_PATH_MAIN + 'data/smcpp_input/{pop1}_{pop2}/{pop1}.{pop2}.chr{CHROM, \d+}.{pop1}.split.vcf.gz.tbi'),
# temp_pop_file1 = temp(VCF_PATH_MAIN + 'data/smcpp_input/{pop1}_{pop2}/{pop1}.chr{CHROM, \d+}.{pop1}.split.list'),
# temp_pop_file2 = temp(VCF_PATH_MAIN + 'data/smcpp_input/{pop1}_{pop2}/{pop2}.chr{CHROM, \d+}.{pop2}.split.list'),
# smc_pop1_out1 = VCF_PATH_MAIN + 'smcpp_format/split/{pop1}_{pop2}/{pop1}.{pop2}.chr{CHROM, \d+}.split.smcpp.gz',
# smc_pop2_out2 = VCF_PATH_MAIN + 'smcpp_format/split/{pop1}_{pop2}/{pop2}.{pop1}.chr{CHROM, \d+}.split.smcpp.gz'
# run:
# shell("bcftools merge -R {input.inclusion_mask} {input.vcf1} {input.vcf2} | bgzip -c > {output.temp_vcf_merged}; tabix {output.temp_vcf_merged}")
# shell("tac {input.panel1} | awk \'{{print $1}}\' > {output.temp_pop_file1}")
# shell("tac {input.panel2} | awk \'{{print $1}}\' > {output.temp_pop_file2}")
# pop_str1 = ""
# with open(output.temp_pop_file1, 'r') as f:
# for line in f:
# pop_str1 += line.rstrip() + ','
# pop_str1 = pop_str1.rstrip(',')
# print(pop_str1)
# pop_str2 = ""
# with open(output.temp_pop_file2, 'r') as f:
# for line in f:
# pop_str2 += line.rstrip() + ','
# pop_str2 = pop_str2.rstrip(',')
# print(pop_str2)
# contig_len = ""
# with open(input.contig_length, 'r') as f:
# for line in f:
# l_splt = line.split()
# print(l_splt)
# if l_splt[0] == wildcards.CHROM:
# contig_len = l_splt[1]; break
# shell("smc++ vcf2smc --length {contig_len} -c 50000 {output.temp_vcf_merged} {output.smc_pop1_out1} {wildcards.CHROM} {wildcards.pop1}:{pop_str1} {wildcards.pop2}:{pop_str2}")
# shell("smc++ vcf2smc --length {contig_len} -c 50000 {output.temp_vcf_merged} {output.smc_pop2_out2} {wildcards.CHROM} {wildcards.pop2}:{pop_str2} {wildcards.pop1}:{pop_str1}")
# '''
# SMC++ for estimating split times between populations
# '''
# rule run_smcpp_split:
# input:
# smc_model1 = "data/smcpp_output_mult_final/{pop1}_t1_{t1}_knots_{knots}_filt/model.final.json",
# smc_model2 = "data/smcpp_output_mult_final/{pop2}_t1_{t1}_knots_{knots}_filt/model.final.json",
# smc_joint_files1 = expand(VCF_PATH_MAIN + 'smcpp_format/split/{{pop1}}_{{pop2}}/{{pop1}}.{{pop2}}.chr{CHROM}.split.smcpp.gz', CHROM=CHROMS),
# smc_joint_files2 = expand(VCF_PATH_MAIN + 'smcpp_format/split/{{pop1}}_{{pop2}}/{{pop2}}.{{pop1}}.chr{CHROM}.split.smcpp.gz', CHROM=CHROMS)
# output:
# smc_split_out = "data/smcpp_output_mult_final/split/{pop1}_{pop2}/t1_{t1}_knots_{knots}/model.final.json"
# run:
# shell("smc++ split -o data/smcpp_output_mult_final/split/{wildcards.pop1}_{wildcards.pop2}/t1_{wildcards.t1}_knots_{wildcards.knots}/ {input.smc_model1} {input.smc_model2}  {VCF_PATH_MAIN}/smcpp_format/split/{wildcards.pop1}_{wildcards.pop2}/*.smcpp.gz")
# '''
# Plotting the split time-estimation
# '''
# rule plot_smcpp_split:
# input:
# smc_split = "data/smcpp_output_mult_final/split/{pop1}_{pop2}/t1_{t1}_knots_{knots}/model.final.json"
# output:
# smc_split_csv = "data/smcpp_output_mult_final/split/{pop1}_{pop2}/t1_{t1}_knots_{knots}/{pop1}_{pop2}_t1_{t1}_knots_{knots}.csv",
# smc_split_pdf = "data/smcpp_output_mult_final/split/{pop1}_{pop2}/t1_{t1}_knots_{knots}/{pop1}_{pop2}_t1_{t1}_knots_{knots}.pdf"
# run:
# shell('smc++ plot {output.smc_split_pdf} {input.smc_split} --csv')
