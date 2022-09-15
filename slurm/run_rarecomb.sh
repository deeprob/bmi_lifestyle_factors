#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=rarecomb_lf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=300G
#SBATCH --chdir /data5/deepro/ukbiobank/analysis/lifestyle_factors/src
#SBATCH -o /data5/deepro/ukbiobank/analysis/lifestyle_factors/slurm/logs/out_rarecomb_multi.log
#SBATCH -e /data5/deepro/ukbiobank/analysis/lifestyle_factors/slurm/logs/err_rarecomb_multi.log
#SBATCH --nodelist=qingyu

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data5/deepro/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data5/deepro/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/data5/deepro/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data5/deepro/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<


conda activate lifestyle

echo `date` starting job on $HOSTNAME

python /data5/deepro/ukbiobank/analysis/lifestyle_factors/src/3_run_rarecomb.py

echo `date` ending job
