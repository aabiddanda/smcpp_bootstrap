#!python3

import gzip as gz
import numpy as np 
import pandas as pd

# ---- Parameter Defs ---- #

CHROMS = [i for i in range(22,0, -1)]

mu=1.25e-8

FOCAL=4

# ---- Useful Functions ---- #

base = lambda x:os.path.splitext(x)[0]
            
# --- Rules --------------------------- # 

'''
    Simulating CEU under the Tennessen et al or Gravel et al demography 
'''
rule simulate_CEU_demography:
    input:
        recomb_map = 'data/genetic_maps/genetic_map_GRCh37_chr{CHROM}.txt'
    output:
       tmp_vcf = temp('data/sims_vcf/CEU_{m, \d+}/ceu_chr{CHROM, \d+}_n{n,\d+}_filt.vcf'),
       vcf = 'data/sims_vcf/CEU_{m,\d+}/ceu_chr{CHROM, \d+}_n{n,\d+}_filt.vcf.gz',
       vcf_idx = 'data/sims_vcf/CEU_{m,\d+}/ceu_chr{CHROM, \d+}_n{n,\d+}_filt.vcf.gz.tbi'
    run:
        shell('python src/sim_demo.py -n {wildcards.n} -m {wildcards.m} -o {output.tmp_vcf} -r {input.recomb_map}')
        shell('bgzip < {output.tmp_vcf} > {output.vcf};  tabix {output.vcf}')


'''
    Converts VCF to SMC++ format for a simulated data
'''
rule conv_vcf2smc_1pop_sim:
    input:
        vcf = 'data/sims_vcf/CEU_{m}/ceu_chr{CHROM}_n{n}_filt.vcf.gz',
        contig_length = DATA_PATH + 'contig_lengths.txt',
    output:
        temp_pop_file = temp('data/smcpp_input_sim/CEU_{m,\d+}/focal_{f,\d+}/ceu_chr{CHROM, \d+}_n{n,\d+}.filt.{i,\d+}.list'),
        smc_out = 'data/smcpp_input_sim/CEU_{m, \d+}/focal_{f, \d+}/ceu_chr{CHROM, \d+}_n{n, \d+}.filt.{i, \d+}.smcpp.gz'
    run:
        shell("bcftools query -l {input.vcf} | shuf > {output.temp_pop_file}")
        pop_str = ""
        with open(output.temp_pop_file, 'r') as f:
            for line in f:
                    pop_str += line.rstrip() + ','
            pop_str = pop_str.rstrip(',')
            print(pop_str)
	    contig_len = ""
	    with open(input.contig_length, 'r') as f:
                for line in f:
                    l_splt = line.split()
                    print(l_splt)
		    if l_splt[0] == wildcards.CHROM:
                        contig_len = l_splt[1]
		        break
            shell("smc++ vcf2smc --length {contig_len} -c 50000 {input.vcf} {output.smc_out} 1 CEU:{pop_str}")



'''
    SMC++ in the 1-population setting for estimation (1-focal individual)
'''
rule run_smcpp_estimate_sim:
    input:
        smc_files = expand('data/smcpp_input_sim/CEU_{{m}}/focal_{{f}}/ceu_chr{CHROM}_n{{n}}.filt.{i}.smcpp.gz', CHROM=CHROMS, i = [x for x in range(FOCAL)])
    output:
        smc_out = "data/smcpp_output_sim/CEU_n{n}_m{m,\d+}/focal_{f,\d+}/t1_{t1,\d+}_knots_{knots,\d+}/model.final.json",
        smc_pdf = "data/smcpp_output_sim/CEU_n{n}_m{m,\d+}/focal_{f,\d+}/t1_{t1}_knots_{knots}/ceu_t1_{t1}_knots_{knots}_filt.pdf",
        smc_csv = "data/smcpp_output_sim/CEU_n{n}_m{m,\d+}/focal_{f,\d+}/t1_{t1}_knots_{knots}/ceu_t1_{t1}_knots_{knots}_filt.csv"
    run:
        n_t1 = ''; n_knot = '';
        if wildcards.t1 != 'None':
            n_t1 = '--t1 ' + wildcards.t1
        if wildcards.knots != 'None':
            n_knots = '--knots ' + wildcards.knots
        shell("smc++ estimate -o data/smcpp_output_sim/CEU_n{wildcards.n}_m{wildcards.m}/focal_{wildcards.f}/t1_{wildcards.t1}_knots_{wildcards.knots}/ -v {mu} {n_t1} {n_knots} data/smcpp_input_sim/CEU_{wildcards.m}/focal_{wildcards.f}/ceu_chr[0-9]*_n{wildcards.n}.filt.*.smcpp.gz")
        shell("smc++ plot {output.smc_pdf} {output.smc_out} --csv")

