#!/bin/bash

#SBATCH -p short
#SBATCH -t 8:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -n 1
#SBATCH -e snakemake.err
#SBATCH -o snakemake.log
#SBATCH -J ChIP-exo-snakemake

snakemake -p --latency-wait 300 --rerun-incomplete --cluster-config cluster.yaml --use-conda --jobs 999 --cluster "sbatch -p {cluster.queue} -n {cluster.n} -t {cluster.time} --mem-per-cpu={cluster.mem} -J {cluster.name} -e {cluster.err} -o {cluster.log}"