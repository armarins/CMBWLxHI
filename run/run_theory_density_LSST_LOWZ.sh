#!/bin/bash
NAME_RUN="LSST_40tomos_0.001sigmaz"
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
SIGMAZ_COUNT=0.001
NTOMOS=40
#
FREQ_MIN=980
FREQ_MAX=1260
NCHANNELS=30
#
GENERATE_CH=0
GENERATE_CG=1
GENERATE_CX=0
##############
RUN_NAME='generating_APS_theory'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"

#############
eval "$(conda shell.bash hook)"
conda activate  /data/AMARINS/anaconda3/PY10

############
cd $DIR_SCRIPTS
bin_min=9
bin_max=39
for i in $(seq $bin_min $bin_max)
do
    PREFIX_new="${PREFIX}${i}_${NTOMOS}_0001"
    python3 $PATHFILE --verbose $VERBOSE --pathout $PATHOUT --prefix $PREFIX_new --suffix $SUFFIX\
                      --field_1 $FIELD1 --field_2 $FIELD2\
                      --generate_cl_field_2 $GENERATE_CG \
                      --survey_count $SURVEY_COUNT --dataset_count $DATASET_COUNT --mode_count $MODE_COUNT\
                      --bias_model $BIAS_MODEL --binning_count $BINNING\
                      --sigmaz_count $SIGMAZ_COUNT --nbins_count $NTOMOS\
                      --bin_to_use_count ${i} --path_to_count_params_file $PATH_TO_COUNT_PARAMS_FILE\
                      --freq_min $FREQ_MIN --freq_max $FREQ_MAX  --nchannels $NCHANNELS | tee -a $FILEOUT
done
echo $TIMEF | tee -a $FILEOUT



