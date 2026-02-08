#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:09:00
#SBATCH --output=sbatch/problem1_v2_i1/SLURM-mandelbrot-C32.log
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=32
#SBATCH --job-name=mandelbrot-C32
#SBATCH --distribution=block:block
#SBATCH --constraint=skylake
#SBATCH -p medium
./prog1_mandelbrot_threads/mandelbrot -t 32 -i 1 -v 2 2> logs/mandelbrot-C32.log
