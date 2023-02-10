#!python3


"""
    Bootstrap Analyses
"""

import gzip as gz
import random

CHROMS = [i for i in range(22, 0, -1)]

BOOTS = [i for i in range(20, 0, -1)]

REF_GENOME = (
    "/home/abiddanda/novembre_lab/data/external_public/reference_genomes/hs37d5.fa"
)

DATA_PATH = "/home/abiddanda/novembre_lab/abiddanda/smcpp_sardinia/data/smcpp_forArjun/"

VCF_PATH = "/home/abiddanda/novembre_lab/share/smcpp_input/20170816/"


mu = 1.25e-8

# ---- Useful Functions ---- #

base = lambda x: os.path.splitext(x)[0]


# --- Rules --------------------------- #


"""
    Creates files for bootstrapping over when conducting replicates
"""


rule create_bootstrap_files_supp:
    input:
        smc_out="data/smcpp_output_mult_final_supp/{POP}_t1_{t1}_knots_{knots}_filt/model.final.json",
    output:
        expand(
            VCF_PATH
            + "smcpp_format/bootstrap/{{POP}}_t1_{{t1}}_knots_{{knots}}_filt_{{n}}/run_1/bootstrap_chr{i}.gz",
            i=CHROMS,
        ),
    run:
        smcpp_dir = VCF_PATH + "smcpp_format/{wildcards.POP}/"
        out_dir = (
            VCF_PATH
            + "smcpp_format/bootstrap/{wildcards.POP}_t1_{wildcards.t1}_knots_{wildcards.knots}_filt_{wildcards.n}/run"
        )
        shell(
            "python3 src/bootstrap_smcpp.py --nr_bootstraps 1  --out_dir_prefix %s %s/*.smcpp.gz"
            % (out_dir, smcpp_dir)
        )

        """
        Estimate using bootstrap data
        """


rule bootstrap_estimate_supp:
    input:
        expand(
            VCF_PATH
            + "smcpp_format/bootstrap/{{POP}}_t1_{{t1}}_knots_{{knots}}_filt_{{n}}/run_1/bootstrap_chr{i}.gz",
            i=CHROMS,
        ),
    output:
        smc_bootstrap_out="data/smcpp_output_mult_final_supp/bootstraps/{POP}/{POP}_t1_{t1}_knots_{knots}_filt_boots_{n}/model.final.json",
        smc_bootstrap_pdf="data/smcpp_output_mult_final_supp/bootstraps/{POP}/{POP}_t1_{t1}_knots_{knots}_filt_boots_{n}/{POP}_t1_{t1}_knots_{knots}_filt_boots_{n}.pdf",
        smc_bootstrap_csv="data/smcpp_output_mult_final_supp/bootstraps/{POP}/{POP}_t1_{t1}_knots_{knots}_filt_boots_{n}/{POP}_t1_{t1}_knots_{knots}_filt_boots_{n}.csv",
    run:
        n_t1 = ""
        n_knots = ""
        if wildcards.t1 != "None":
            n_t1 = "--t1 " + wildcards.t1
        if wildcards.knots != "None":
            n_knots = "--knots " + wildcards.knots
        out_dir = "data/smcpp_output_mult_final_supp/bootstraps/{wildcards.POP}/{wildcards.POP}_t1_{wildcards.t1}_knots_{wildcards.knots}_filt_boots_{wildcards.n}/"
        shell("smc++ estimate -o %s {n_t1} {n_knots} -v {mu} {input}" % out_dir)
        shell("smc++ plot {output.smc_bootstrap_pdf} {output.smc_bootstrap_out} --csv")


        """
        Create bootstrap estimates of the split times
        """


rule create_bootstrap_split_supp:
    input:
        pop1_pop2=expand(
            VCF_PATH
            + "smcpp_format/split/{{POP1}}_{{POP2}}/{{POP1}}.{{POP2}}.chr{CHROM}.split.smcpp.gz",
            CHROM=CHROMS,
        ),
    output:
        pop1_pop2=expand(
            VCF_PATH
            + "smcpp_format/bootstrap/split/{{POP1}}_{{POP2}}/t1_{{t1}}_knots_{{knots}}_filt_{{n}}/run_1/bootstrap_chr{i}.gz",
            i=CHROMS,
        ),
    run:
        smcpp_dir = VCF_PATH + "smcpp_format/split/{wildcards.POP1}_{wildcards.POP2}/"
        out_dir1 = (
            VCF_PATH
            + "smcpp_format/bootstrap/split/{wildcards.POP1}_{wildcards.POP2}/t1_{wildcards.t1}_knots_{wildcards.knots}_filt_{wildcards.n}/run"
        )
        shell(
            "python3 src/bootstrap_smcpp.py --nr_bootstraps 1 --out_dir_prefix %s {input.pop1_pop2}"
            % (out_dir1)
        )

        """
        Estimating split times for the supplementary datasets
        """


rule bootstrap_split_supp:
    input:
        smcpp_input1=expand(
            VCF_PATH
            + "smcpp_format/bootstrap/split/{{POP1}}_{{POP2}}/t1_{{t1}}_knots_{{knots}}_filt_{{n}}/run_1/bootstrap_chr{i}.gz",
            i=CHROMS,
        ),
        smcpp_input2=expand(
            VCF_PATH
            + "smcpp_format/bootstrap/split/{{POP2}}_{{POP1}}/t1_{{t1}}_knots_{{knots}}_filt_{{n}}/run_1/bootstrap_chr{i}.gz",
            i=CHROMS,
        ),
        smc_model1="data/smcpp_output_mult_final_supp/{POP1}_t1_{t1}_knots_{knots}_filt/model.final.json",
        smc_model2="data/smcpp_output_mult_final_supp/{POP2}_t1_{t1}_knots_{knots}_filt/model.final.json",
    output:
        smc_bootstrap_out="data/smcpp_output_mult_final_supp/split/bootstraps/{POP1}_{POP2}/t1_{t1}_knots_{knots}_filt/boots_{n}/model.final.json",
    run:
        out_dir = "data/smcpp_output_mult_final_supp/split/bootstraps/{wildcards.POP1}_{wildcards.POP2}/t1_{wildcards.t1}_knots_{wildcards.knots}_filt/boots_{wildcards.n}/"
        shell(
            "smc++ split -o %s  {input.smc_model1} {input.smc_model2} {input.smcpp_input1} {input.smcpp_input2}"
            % out_dir
        )


        """
            Extract the likelihood of each split-time Estimate
        rule bootstrap_likelihood_split:
            input:

                grep "Loglik\|Split" */.debug.txt | awk -F":" '{print $3}' | paste -d "\t" - -

        """
