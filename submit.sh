#! /bin/bash
#BSUB -J filter
#BSUB -n 8
#BSUB -R "span[ptile=8]"
#BSUB -R "rusage[mem=16]"
#BSUB -W 12:00
#BSUB -o %J.stdout
#BSUB -eo %J.stderr

cd /data/leslie/hy395/deconv/capselex
Rscript scripts/step3_filter.R
