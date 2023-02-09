#!python3

'''
    Sanity-Checking Analyses 
'''

import gzip as gz
import random


CHROMS = [i for i in range(22,0, -1)]

VCF_PATH_SANITY = "/home/abiddanda/novembre_lab/share/smcpp_input/20170809/"


# ---- Useful Functions ---- #

base = lambda x:os.path.splitext(x)[0]


# ---- Definitions ---- #

CEU_INDIVS = ['NA12878', 'NA12889', 'NA12890', 'NA12891']

TSI_INDIVS = ['NA20502', 'NA20509', 'NA20510', 'NA20511']


# --- Rules --------------------------- # 


'''
    Check the SFS of distinguished individuals vs non-distinguished 
'''
rule check_sfs:
    input:
        vcf = VCF_PATH_SANITY  + '{POP}.chr{CHROM}.smcpp.primary.vcf.gz',
        vcf_idx = VCF_PATH_SANITY + '{POP}.chr{CHROM}.smcpp.primary.vcf.gz.tbi',
        samplist = VCF_PATH_SANITY + 'aux.indlist.sampleid.primary.smcpp.{POP}.txt'
    output:
        focal_inds = temp('data/sanity/{POP}_{CHROM}_focal.ind'), 
        non_focal_inds = temp('data/sanity/{POP}_{CHROM}_nonfocal.ind'),
        focal_miss_sites = temp('data/sanity/{POP}_{CHROM}_focal_miss.sites'),
        non_focal_miss_sites = temp('data/sanity/{POP}_{CHROM}_nonfocal_miss.sites'),
        sfs_nonfocal_miss = 'data/sanity/sfs/{POP}_{CHROM}_focal_miss.sfs',
        sfs_focal_miss = 'data/sanity/sfs/{POP}_{CHROM}_nonfocal_miss.sfs'
    run:
        shell('awk \'$3 == \"DEEPSEQ\" {{print $1}} \' {input.samplist} > {output.focal_inds}')
        shell('awk \'$3 != \"DEEPSEQ\" {{print $1}} \' {input.samplist} > {output.non_focal_inds}')
        shell('vcftools --gzvcf {input.vcf} --min-alleles 2 --max-alleles 2  --keep {output.non_focal_inds} --missing-site --stdout | awk -v OFS=\'\t\' \'NR > 1 && $6 == 1 {{print $1,$2}}\' > {output.non_focal_miss_sites} ')
        shell('vcftools --gzvcf {input.vcf} --min-alleles 2 --max-alleles 2 --keep {output.focal_inds} --missing-site --stdout | awk -v OFS=\'\t\' \'NR > 1 && $6 == 1 {{print $1,$2}}\' > {output.focal_miss_sites} ')
        shell('vcftools --gzvcf {input.vcf} --keep {output.focal_inds} --positions {output.non_focal_miss_sites} --counts2 --stdout > {output.sfs_nonfocal_miss}')
        shell('vcftools --gzvcf {input.vcf} --keep {output.non_focal_inds} --positions {output.focal_miss_sites} --counts2 --stdout > {output.sfs_focal_miss}')


rule check_split_seeds:
    input:
        smcpp_input = expand(VCF_PATH_SANITY + 'smcpp_format/split/{{POP1}}_{{POP2}}/{{POP1}}.{{POP2}}.chr{i}.split.smcpp.gz', i = CHROMS),
        smc_model1 = "data/smcpp_output_mult_final/{POP1}_t1_{t1}_knots_{knots}_filt/model.final.json",
        smc_model2 = "data/smcpp_output_mult_final/{POP2}_t1_{t1}_knots_{knots}_filt/model.final.json"
    output:
        smc_out = "data/smcpp_output_split_check/{POP1}_{POP2}_t1_{t1}_knots_{knots}/run_{n}/model.final.json"
    run:
        out_dir = "data/smcpp_output_split_check/{wildcards.POP1}_{wildcards.POP2}_t1_{wildcards.t1}_knots_{wildcards.knots}/run_{wildcards.n}/"
        data_dir = "{VCF_PATH_SANITY}smcpp_format/split/{wildcards.POP1}_{wildcards.POP2}"
        shell("smc++ split --seed %d --ftol 0.0001  -o %s {input.smc_model1} {input.smc_model2} %s/*.smcpp.gz " % (random.randint(0,1000), out_dir, data_dir))


