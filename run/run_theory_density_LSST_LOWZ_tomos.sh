#!/bin/bash
NAME_RUN="LSSTxLOWZ_DENSITY"
DIR_OUTSCREEN='/data/AMARINS/LSSxHI-CODES/screen_outputs'
DIR_SCRIPTS='/data/AMARINS/LSSxHI-CODES/scripts'

FILEOUT="$DIR_OUTSCREEN/$NAME_RUN"
TIMEI=$(date +%Y-%m-%d-%H:%M:%S)
FILEOUT="$FILEOUT.$TIMEI"
FILEOUT="$FILEOUT.out"

echo $TIMEI | tee $FILEOUT

#######
#TERMINAL INFO
VERBOSE=1
#
FIELD1='hi'
FIELD2='density'
PATHOUT='/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI'
PREFIX="LSST_bin"
SUFFIX="LOWZ"
#
#BINNING=True
PATH_TO_COUNT_PARAMS_FILE="/data/AMARINS/LSSxHI-CODES/scripts/parameters/lsst_desc_parameters.yaml"
SURVEY_COUNT="LSST"
DATASET_COUNT="Y10"
MODE_COUNT='photometric'
BIAS_MODEL='SQRT'
BINNING=1
#
FREQ_MIN=980
FREQ_MAX=1260
NCHANNELS=30

##############
RUN_NAME='generating_APS_theory'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"
#############
eval "$(conda shell.bash hook)"
conda activate  /data/AMARINS/anaconda3/amarins3.10
############
for SIGMA_Z in 0.05 0.01 0.001; do 
    for NBINS in 5 10 40; do 
        max_j=$((NBINS - 1))  # Calculate i + 5
        for j in $(seq 0 $max_j); do
            PREFIX_new="${PREFIX}${j}_${NBINS}tomos_${SIGMA_Z}sigmaz"
            python $PATHFILE --verbose $VERBOSE --pathout $PATHOUT --prefix $PREFIX_new --suffix $SUFFIX\
                             --field_1 $FIELD1 --field_2 $FIELD2\
                             --survey_count $SURVEY_COUNT --dataset_count $DATASET_COUNT --mode_count $MODE_COUNT\
                             --bias_model $BIAS_MODEL --binning_count $BINNING\
                             --bin_to_use_count ${j} --path_to_count_params_file $PATH_TO_COUNT_PARAMS_FILE\
                             --freq_min $FREQ_MIN --freq_max $FREQ_MAX  --nchannels $NCHANNELS\
                             --nbins_count ${NBINS} --sigmaz_count ${SIGMA_Z}| tee -a $FILEOUT
        done
    done
done    
echo $TIMEF | tee -a $FILEOUT