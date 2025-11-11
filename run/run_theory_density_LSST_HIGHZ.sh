#!/bin/bash
NAME_RUN="HIGHZ_DENSITY_LSSTY10_256_nch30"
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
PREFIX="LSSTY10_bin"
SUFFIX="nch70_350_1050"
#
#BINNING=True
PATH_TO_COUNT_PARAMS_FILE="/data/AMARINS/CMBWLxHI-CODES/scripts/parameters/lsst_desc_parameters.yaml"
SURVEY_COUNT="LSST"
DATASET_COUNT="Y10"
MODE_COUNT='photometric'
BIAS_MODEL='SQRT'
BINNING=1
#
FREQ_MIN=350
FREQ_MAX=1050
NCHANNELS=70

##############
RUN_NAME='generating_APS_theory'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"

#############
eval "$(conda shell.bash hook)"
conda activate  /home/amarins/.conda/envs/amarins_camb

############
cd $DIR_SCRIPTS
bin_min=0
bin_max=4
for i in $(seq $bin_min $bin_max)
do
    PREFIX_new="${PREFIX}${i}"
    python3 $PATHFILE --verbose $VERBOSE --pathout $PATHOUT --prefix $PREFIX_new --suffix $SUFFIX\
                      --field_1 $FIELD1 --field_2 $FIELD2\
                      --survey_count $SURVEY_COUNT --dataset_count $DATASET_COUNT --mode_count $MODE_COUNT\
                      --bias_model $BIAS_MODEL --binning_count $BINNING\
                      --bin_to_use_count ${i} --path_to_count_params_file $PATH_TO_COUNT_PARAMS_FILE\
                      --freq_min $FREQ_MIN --freq_max $FREQ_MAX  --nchannels $NCHANNELS | tee -a $FILEOUT
done
echo $TIMEF | tee -a $FILEOUT



