#!/bin/bash
# set -ue # don't set ue because it interrrupts with conda activate

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

# activate conda environment
conda activate rarecomb

input_table=$1
primary_entities=$2
secondary_entities=$3
output_table=$4
combos=$5

Rscript /data5/deepro/ukbiobank/analysis/lifestyle_factors/src/scripts/run_rarecomb.R $input_table $primary_entities $secondary_entities $output_table $combos
