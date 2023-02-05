#!/bin/bash

#SBATCH -A PAS0471
#SBATCH --time=300
#SBATCH --output=slurm-copy-fastq-%j.out

from=/fs/scratch/PAS0472/data-transfer/soledad_b/210430_Pearlly_GSL-PY-2114-transfer/
to=/fs/project/PAS0471/rawalranjana44/data

date
echo "## Starting to copy..."

cp -rv "$from" "$to"

echo -e "\n## Size of source dir:"
du -sh "$from"

echo -e "## Size of target dir:"
du -sh "$to"

echo -e "\n ## Checking MD5 sums:"
cd "$to" || exit
md5sum -c md5.txt

echo -e "\n## Done with script."