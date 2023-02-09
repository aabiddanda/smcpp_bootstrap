#!/usr/local/bin/python3
"""
    Created on Fri Jul 28 15:22:36 2017
"""

import argparse as arg
import os
import sys
import random
import gzip
import time


'''
    Function to bootstrap over chunks of chromosomes
'''
def main(nr_bootstraps, chunk_size, chunks_per_chromosome, nr_chromosomes, seed, out_dir_prefix, files):
    chunks = []
    offset = 0
    chunks_in_chrom = []
    if not seed:
        seed = int(time.time())
    random.seed(seed)
    print('seed: %s' % seed)
    for file in files:
        with gzip.open(file, 'rb') as f:
            header = f.readline().decode()
            for line in f:
                line = line.decode()
                line = [int(x) for x in line.strip().split()]
                pos = line[0] + offset
                offset += line[0]
                chunk_index = (pos - 1) // chunk_size
                if chunk_index > len(chunks_in_chrom)-1:
                    for _ in range(chunk_index - len(chunks_in_chrom)+1):
                        chunks_in_chrom.append([])
                chunks_in_chrom[chunk_index].append(line)
        chunks.extend(chunks_in_chrom)

    for bootstrap_id in range(1, nr_bootstraps +1):
        for chr_ in range(1, nr_chromosomes + 1):
            chr_dir = "{}_{}".format(out_dir_prefix, bootstrap_id)
            if not os.path.exists(chr_dir):
                os.makedirs(chr_dir)
            chr_file = "{}/bootstrap_chr{}.gz".format(chr_dir, chr_)
            print("writing", chr_file, file=sys.stderr)
            with gzip.open(chr_file, 'wb') as f:
                f.write(header.encode())
                for i in range(chunks_per_chromosome):
                    chunk_id = random.randrange(len(chunks))
                    for line in chunks[chunk_id]:
                        line = ' '.join([str(x) for x in line]) + '\n'
                        f.write(line.encode())



if __name__ == "__main__":
    parser = arg.ArgumentParser()
    parser.add_argument('--nr_bootstraps', type=int, help="nr of bootstraps [10]", default=10)
    parser.add_argument("--chunk_size", type=int, help="size of bootstrap chunks [10000000]", default=10000000)
    parser.add_argument("--chunks_per_chromosome", type=int, help="nr of chunks to put on one chromosome in the bootstrap [13]", default=13)
    parser.add_argument("--nr_chromosomes", type=int, help="nr of chromosomes to write [22]", default=22)
    parser.add_argument("--seed", type=int, help="initialize the random number generator", default=None)
    parser.add_argument("--out_dir_prefix", required=True, help="prefix for the outdirectory")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    main(args.nr_bootstraps, args.chunk_size, args.chunks_per_chromosome, args.nr_chromosomes, args.seed, args.out_dir_prefix, args.files)
