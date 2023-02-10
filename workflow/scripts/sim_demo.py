#!/usr/local/bin/python3

"""
    Simulating ground-truth demographies to test
        parameter values of SMC++
"""

import argparse as arg

import msprime as msp

"""
    Simulating from the Tennessen et al (2012) demography
"""


def tennessen_et_al(n=100, recomb_map=None):
    # Defining intial sizes and ancestral sizes for the population
    N_CEU0 = 5.12e6
    N_A = 14600

    population_configurations = [
        msp.PopulationConfiguration(sample_size=n, initial_size=N_CEU0),
    ]
    demographic_events = [
        msp.PopulationParametersChange(time=0, initial_size=N_CEU0, growth_rate=0.0195),
        msp.PopulationParametersChange(
            time=204, initial_size=19930.67, growth_rate=0.00307
        ),
        msp.PopulationParametersChange(time=920, initial_size=3722, growth_rate=0.0),
        msp.PopulationParametersChange(time=2040, initial_size=28948, growth_rate=0),
        msp.PopulationParametersChange(time=5920, growth_rate=0, initial_size=N_A),
    ]
    if recomb_map is not None:
        tree_sequence = msp.simulate(
            sample_size=n,
            recombination_map=recomb_map,
            demographic_events=demographic_events,
            mutation_rate=1.2e-8,
        )
    else:
        tree_sequence = msp.simulate(
            sample_size=n,
            length=10e6,
            demographic_events=demographic_events,
            mutation_rate=1.2e-8,
        )
    return tree_sequence


"""
    Simulating from the Tennessen et al (2012) demography
"""


def gravel_et_al(n=100, recomb_map=None):
    # Defining intial sizes and ancestral sizes for the population
    N_CEU0 = 67627
    N_A = 14600

    population_configurations = [
        msp.PopulationConfiguration(sample_size=n, initial_size=N_CEU0),
    ]
    demographic_events = [
        msp.PopulationParametersChange(time=0, initial_size=N_CEU0, growth_rate=0.0038),
        msp.PopulationParametersChange(time=920, initial_size=3722, growth_rate=0.0),
        msp.PopulationParametersChange(time=2040, initial_size=28948, growth_rate=0),
        msp.PopulationParametersChange(time=5920, growth_rate=0, initial_size=N_A),
    ]
    if recomb_map is not None:
        tree_sequence = msp.simulate(
            sample_size=n,
            recombination_map=recomb_map,
            demographic_events=demographic_events,
            mutation_rate=1.2e-8,
        )
    else:
        tree_sequence = msp.simulate(
            sample_size=n,
            length=10e6,
            demographic_events=demographic_events,
            mutation_rate=1.2e-8,
        )
    return tree_sequence


if __name__ == "__main__":
    # Parse all arguments given
    parser = arg.ArgumentParser()
    parser.add_argument(
        "-m",
        "--model",
        required=True,
        default=1,
        type=int,
        help="Choosing the demographic model to use",
    )
    parser.add_argument("-n", required=True, type=int, help="Diploid Sample Size")
    parser.add_argument(
        "-r", "--recomb", required=True, help="Recombination Map (Hapmap format)"
    )
    parser.add_argument(
        "-o", "--out", required=True, help="output VCF file for the variation data"
    )
    args = parser.parse_args()

    # Reading in the recombination rate file
    hapmap_recomb_map = msp.RecombinationMap.read_hapmap(args.recomb)
    print("Finished reading in ")

    # Simulating a tree sequence under the Tennessen et al demography
    if args.model == 1:
        tree_seq = tennessen_et_al(n=args.n, recomb_map=hapmap_recomb_map)
    if args.model == 2:
        tree_seq = gravel_et_al(n=args.n, recomb_map=hapmap_recomb_map)

    if args.model not in [1, 2]:
        print("No appropriate model chosen!")
        quit()
    print("Finished generating tree sequence!")

    with open(args.out, "w+") as vcf_file:
        tree_seq.write_vcf(vcf_file, ploidy=2)
