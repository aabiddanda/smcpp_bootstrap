#!python3

"""
    Supplementary Analyses
"""

import gzip as gz

CHROMS = [i for i in range(22, 0, -1)]

REF_GENOME = (
    "/home/abiddanda/novembre_lab/data/external_public/reference_genomes/hs37d5.fa"
)

DATA_PATH = "/home/abiddanda/novembre_lab/abiddanda/smcpp_sardinia/data/smcpp_forArjun/"

VCF_PATH_SUPP = "/home/abiddanda/novembre_lab/share/smcpp_input/20170816/"


mu = 1.25e-8

# ---- Useful Functions ---- #

base = lambda x: os.path.splitext(x)[0]


# ---- Definitions ---- #

CEU_INDIVS = ["NA12878", "NA12889"]

TSI_INDIVS = ["NA20502", "NA20509"]

ARZANA_INDIVS = ["Sample_7972", "Sample_33419"]

LANUSEI_INDIVS = ["Sample_8053", "Sample_12320"]

CAGLIARI_INDIVS = ["Sample_33645", "Sample_8700"]

ILBONO_INDIVS = ["Sample_10359", "Sample_8141"]

POP_2_INDIV = {
    "CEU": CEU_INDIVS,
    "TSI": TSI_INDIVS,
    "LANUSEI": LANUSEI_INDIVS,
    "ARZANA": ARZANA_INDIVS,
    "ILBONO": ILBONO_INDIVS,
    "CAGLIARI": CAGLIARI_INDIVS,
}


# --- Rules --------------------------- #


"""
    Converts files for multiple focal individuals at once
"""


rule conv_vcf2smc_1pop_mult_supp:
    input:
        vcf=VCF_PATH_SUPP + "{POP}.chr{CHROM}.smcpp.supplemental.vcf.gz",
        panel=VCF_PATH_SUPP + "aux.indlist.sampleid.supplemental.smcpp.{POP}.txt",
        inclusion_mask=DATA_PATH + "mappability_masks/hs37d5_chr{CHROM}.mask.bed.gz",
        contig_length=DATA_PATH + "contig_lengths.txt",
    output:
        vcf_filt=temp(
            VCF_PATH_SUPP + "smcpp_format/{POP}/{POP}.chr{CHROM,\d+}.{focal}.vcf.gz"
        ),
        vcf_filt_idx=temp(
            VCF_PATH_SUPP
            + "smcpp_format/{POP}/{POP}.chr{CHROM,\d+}.{focal}.vcf.gz.tbi"
        ),
        smc_out=VCF_PATH_SUPP
        + "smcpp_format/{POP}/{POP}.chr{CHROM,\d+}.{focal}.smcpp.gz",
    run:
        shell(
            "bcftools view -R {input.inclusion_mask} {input.vcf} | bgzip -c > {output.vcf_filt}; tabix {output.vcf_filt}"
        )
        ind_list = []
        pop_str = ""
        with open(input.panel, "r") as f:
            for line in f:
                id1, id2, s = line.split()
                pop_str += id1 + ","
                if s == "DEEPSEQ":
                    ind_list.append(id1)
        pop_str = pop_str.rstrip(",")
        contig_len = ""
        with open(input.contig_length, "r") as f:
            for line in f:
                l_splt = line.split()
                print(l_splt)
                if l_splt[0] == wildcards.CHROM:
                    contig_len = l_splt[1]
                    break
        shell(
            "smc++ vcf2smc --length {contig_len} -c 50000 -d {wildcards.focal} {wildcards.focal} {output.vcf_filt} {output.smc_out} {wildcards.CHROM} {wildcards.POP}:{pop_str}"
        )


        """
        Estimating for CEU
        """


rule run_smcpp_estimate_ceu_supp:
    input:
        smc_files=expand(
            VCF_PATH_SUPP + "smcpp_format/CEU/CEU.chr{CHROM}.{focal}.smcpp.gz",
            CHROM=CHROMS,
            focal=POP_2_INDIV["CEU"],
        ),
    output:
        smc_out="data/smcpp_output_mult_final_supp/CEU_t1_{t1}_knots_{knots}_filt/model.final.json",
    run:
        n_t1 = ""
        n_knots = ""
        if wildcards.t1 != "None":
            n_t1 = "--t1 " + wildcards.t1
        if wildcards.knots != "None":
            n_knots = "--knots " + wildcards.knots
        shell(
            "smc++ estimate -o data/smcpp_output_mult_final_supp/CEU_t1_{wildcards.t1}_knots_{wildcards.knots}_filt/ {n_t1} {n_knots} -v {mu}  {VCF_PATH_SUPP}/smcpp_format/CEU/CEU.*.smcpp.gz"
        )


        """
        Estimating for TSI
        """


rule run_smcpp_estimate_tsi_supp:
    input:
        smc_files=expand(
            VCF_PATH_SUPP + "smcpp_format/TSI/TSI.chr{CHROM}.{focal}.smcpp.gz",
            CHROM=CHROMS,
            focal=POP_2_INDIV["TSI"],
        ),
    output:
        smc_out="data/smcpp_output_mult_final_supp/TSI_t1_{t1}_knots_{knots}_filt/model.final.json",
    run:
        n_t1 = ""
        n_knots = ""
        if wildcards.t1 != "None":
            n_t1 = "--t1 " + wildcards.t1
        if wildcards.knots != "None":
            n_knots = "--knots " + wildcards.knots
        shell(
            "smc++ estimate -o data/smcpp_output_mult_final_supp/TSI_t1_{wildcards.t1}_knots_{wildcards.knots}_filt/ {n_t1} {n_knots} -v {mu}  {VCF_PATH_SUPP}/smcpp_format/TSI/TSI.*.smcpp.gz"
        )

        """
        Estimating for LANUSEI
        """


rule run_smcpp_estimate_lanusei_supp:
    input:
        smc_files=expand(
            VCF_PATH_SUPP + "smcpp_format/LANUSEI/LANUSEI.chr{CHROM}.{focal}.smcpp.gz",
            CHROM=CHROMS,
            focal=POP_2_INDIV["LANUSEI"],
        ),
    output:
        smc_out="data/smcpp_output_mult_final_supp/LANUSEI_t1_{t1}_knots_{knots}_filt/model.final.json",
    run:
        n_t1 = ""
        n_knots = ""
        if wildcards.t1 != "None":
            n_t1 = "--t1 " + wildcards.t1
        if wildcards.knots != "None":
            n_knots = "--knots " + wildcards.knots
        shell(
            "smc++ estimate -o data/smcpp_output_mult_final_supp/LANUSEI_t1_{wildcards.t1}_knots_{wildcards.knots}_filt/ {n_t1} {n_knots} -v {mu}  {VCF_PATH_SUPP}/smcpp_format/LANUSEI/LANUSEI.*.smcpp.gz"
        )


        """
        Estimating for ARZANA
        """


rule run_smcpp_estimate_arzana_supp:
    input:
        smc_files=expand(
            VCF_PATH_SUPP + "smcpp_format/{POP}/{POP}.chr{CHROM}.{focal}.smcpp.gz",
            POP="ARZANA",
            CHROM=CHROMS,
            focal=POP_2_INDIV["ARZANA"],
        ),
    output:
        smc_out="data/smcpp_output_mult_final_supp/ARZANA_t1_{t1}_knots_{knots}_filt/model.final.json",
    run:
        n_t1 = ""
        n_knots = ""
        if wildcards.t1 != "None":
            n_t1 = "--t1 " + wildcards.t1
        if wildcards.knots != "None":
            n_knots = "--knots " + wildcards.knots
        shell(
            "smc++ estimate -o data/smcpp_output_mult_final_supp/ARZANA_t1_{wildcards.t1}_knots_{wildcards.knots}_filt/ {n_t1} {n_knots} -v {mu}  {VCF_PATH_SUPP}/smcpp_format/ARZANA/ARZANA.*.smcpp.gz"
        )


        """
        Estimating for ILBONO
        """


rule run_smcpp_estimate_ilbono_supp:
    input:
        smc_files=expand(
            VCF_PATH_SUPP + "smcpp_format/{POP}/{POP}.chr{CHROM}.{focal}.smcpp.gz",
            POP="ILBONO",
            CHROM=CHROMS,
            focal=POP_2_INDIV["ILBONO"],
        ),
    output:
        smc_out="data/smcpp_output_mult_final_supp/ILBONO_t1_{t1}_knots_{knots}_filt/model.final.json",
    run:
        n_t1 = ""
        n_knots = ""
        if wildcards.t1 != "None":
            n_t1 = "--t1 " + wildcards.t1
        if wildcards.knots != "None":
            n_knots = "--knots " + wildcards.knots
        shell(
            "smc++ estimate -o data/smcpp_output_mult_final_supp/ILBONO_t1_{wildcards.t1}_knots_{wildcards.knots}_filt/ {n_t1} {n_knots} -v {mu}  {VCF_PATH_SUPP}/smcpp_format/ILBONO/ILBONO.*.smcpp.gz"
        )


        """
        Estimating for CAGLIARI
        """


rule run_smcpp_estimate_cagliari_supp:
    input:
        smc_files=expand(
            VCF_PATH_SUPP + "smcpp_format/{POP}/{POP}.chr{CHROM}.{focal}.smcpp.gz",
            POP="CAGLIARI",
            CHROM=CHROMS,
            focal=POP_2_INDIV["CAGLIARI"],
        ),
    output:
        smc_out="data/smcpp_output_mult_final_supp/CAGLIARI_t1_{t1}_knots_{knots}_filt/model.final.json",
    run:
        n_t1 = ""
        n_knots = ""
        if wildcards.t1 != "None":
            n_t1 = "--t1 " + wildcards.t1
        if wildcards.knots != "None":
            n_knots = "--knots " + wildcards.knots
        shell(
            "smc++ estimate -o data/smcpp_output_mult_final_supp/CAGLIARI_t1_{wildcards.t1}_knots_{wildcards.knots}_filt/ {n_t1} {n_knots} -v {mu}  {VCF_PATH_SUPP}/smcpp_format/CAGLIARI/CAGLIARI.*.smcpp.gz"
        )


        """
        Generates CSVs of the final results
        """


rule generate_csvs_supp:
    input:
        smc_out="data/smcpp_output_mult_final_supp/{POP}_t1_{t1}_knots_{knots}_filt/model.final.json",
    output:
        smc_csv="data/smcpp_output_mult_final_supp/{POP}_t1_{t1}_knots_{knots}_filt/{POP}_t1_{t1}_knots_{knots}_filt.csv",
        smc_pdf="data/smcpp_output_mult_final_supp/{POP}_t1_{t1}_knots_{knots}_filt/{POP}_t1_{t1}_knots_{knots}_filt.pdf",
    run:
        shell("smc++ plot {output.smc_pdf} {input.smc_out} --csv")


# ---- Estimating Split Times ---- #

"""
    Converts VCF to SMC++ format for a pair of populations (to calculate splits between populations)
"""


rule conv_vcf2smc_2pop_supp:
    input:
        vcf1=VCF_PATH_SUPP + "{pop1}.chr{CHROM}.smcpp.supplemental.vcf.gz",
        vcf2=VCF_PATH_SUPP + "{pop2}.chr{CHROM}.smcpp.supplemental.vcf.gz",
        panel1=VCF_PATH_SUPP + "aux.indlist.sampleid.supplemental.smcpp.{pop1}.txt",
        panel2=VCF_PATH_SUPP + "aux.indlist.sampleid.supplemental.smcpp.{pop2}.txt",
        inclusion_mask=DATA_PATH + "mappability_masks/hs37d5_chr{CHROM}.mask.bed.gz",
        contig_length=DATA_PATH + "contig_lengths.txt",
    output:
        temp_vcf_merged=temp(
            VCF_PATH_SUPP
            + "data/smcpp_input/{pop1}_{pop2}/{pop1}.{pop2}.chr{CHROM, \d+}.{pop1}.split.vcf.gz"
        ),
        temp_vcf_merged_index=temp(
            VCF_PATH_SUPP
            + "data/smcpp_input/{pop1}_{pop2}/{pop1}.{pop2}.chr{CHROM, \d+}.{pop1}.split.vcf.gz.tbi"
        ),
        temp_pop_file1=temp(
            VCF_PATH_SUPP
            + "data/smcpp_input/{pop1}_{pop2}/{pop1}.chr{CHROM, \d+}.{pop1}.split.list"
        ),
        temp_pop_file2=temp(
            VCF_PATH_SUPP
            + "data/smcpp_input/{pop1}_{pop2}/{pop2}.chr{CHROM, \d+}.{pop2}.split.list"
        ),
        smc_pop1_out1=VCF_PATH_SUPP
        + "smcpp_format/split/{pop1}_{pop2}/{pop1}.{pop2}.chr{CHROM, \d+}.split.smcpp.gz",
        smc_pop2_out2=VCF_PATH_SUPP
        + "smcpp_format/split/{pop1}_{pop2}/{pop2}.{pop1}.chr{CHROM, \d+}.split.smcpp.gz",
    run:
        shell(
            "bcftools merge -R {input.inclusion_mask} {input.vcf1} {input.vcf2} | bgzip -c > {output.temp_vcf_merged}; tabix {output.temp_vcf_merged}"
        )
        shell("tac {input.panel1} | awk '{{print $1}}' > {output.temp_pop_file1}")
        shell("tac {input.panel2} | awk '{{print $1}}' > {output.temp_pop_file2}")
        pop_str1 = ""
        with open(output.temp_pop_file1, "r") as f:
            for line in f:
                pop_str1 += line.rstrip() + ","
        pop_str1 = pop_str1.rstrip(",")
        print(pop_str1)
        pop_str2 = ""
        with open(output.temp_pop_file2, "r") as f:
            for line in f:
                pop_str2 += line.rstrip() + ","
        pop_str2 = pop_str2.rstrip(",")
        print(pop_str2)
        contig_len = ""
        with open(input.contig_length, "r") as f:
            for line in f:
                l_splt = line.split()
                print(l_splt)
                if l_splt[0] == wildcards.CHROM:
                    contig_len = l_splt[1]
                    break
        shell(
            "smc++ vcf2smc --length {contig_len} -c 50000 {output.temp_vcf_merged} {output.smc_pop1_out1} {wildcards.CHROM} {wildcards.pop1}:{pop_str1} {wildcards.pop2}:{pop_str2}"
        )
        shell(
            "smc++ vcf2smc --length {contig_len} -c 50000 {output.temp_vcf_merged} {output.smc_pop2_out2} {wildcards.CHROM} {wildcards.pop2}:{pop_str2} {wildcards.pop1}:{pop_str1}"
        )


        """
        SMC++ for estimating split times between populations
        """


rule run_smcpp_split_supp:
    input:
        smc_model1="data/smcpp_output_mult_final_supp/{pop1}_t1_{t1}_knots_{knots}_filt/model.final.json",
        smc_model2="data/smcpp_output_mult_final_supp/{pop2}_t1_{t1}_knots_{knots}_filt/model.final.json",
        smc_joint_files1=expand(
            VCF_PATH_SUPP
            + "smcpp_format/split/{{pop1}}_{{pop2}}/{{pop1}}.{{pop2}}.chr{CHROM}.split.smcpp.gz",
            CHROM=CHROMS,
        ),
        smc_joint_files2=expand(
            VCF_PATH_SUPP
            + "smcpp_format/split/{{pop1}}_{{pop2}}/{{pop2}}.{{pop1}}.chr{CHROM}.split.smcpp.gz",
            CHROM=CHROMS,
        ),
    output:
        smc_split_out="data/smcpp_output_mult_final_supp/split/{pop1}_{pop2}/t1_{t1}_knots_{knots}/model.final.json",
    run:
        shell(
            "smc++ split -o data/smcpp_output_mult_final_supp/split/{wildcards.pop1}_{wildcards.pop2}/t1_{wildcards.t1}_knots_{wildcards.knots}/ {input.smc_model1} {input.smc_model2}  {VCF_PATH_SUPP}/smcpp_format/split/{wildcards.pop1}_{wildcards.pop2}/*.smcpp.gz"
        )


        """
        Plots the SMC++ split results
        """


rule plot_smcpp_split_supp:
    input:
        smc_split="data/smcpp_output_mult_final_supp/split/{pop1}_{pop2}/t1_{t1}_knots_{knots}/model.final.json",
    output:
        smc_split_csv="data/smcpp_output_mult_final_supp/split/{pop1}_{pop2}/t1_{t1}_knots_{knots}/{pop1}_{pop2}_t1_{t1}_knots_{knots}.csv",
        smc_split_pdf="data/smcpp_output_mult_final_supp/split/{pop1}_{pop2}/t1_{t1}_knots_{knots}/{pop1}_{pop2}_t1_{t1}_knots_{knots}.pdf",
    run:
        shell("smc++ plot {output.smc_split_pdf} {input.smc_split} --csv")
