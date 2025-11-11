#!/bin/bash
NAME_RUN="HIGHZ_SHEAR_LSST_256_nch70"
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
FIELD2='shear'
PATHOUT='/data/AMARINS/LSSxHI-DATA/theoretical/GWLxHI'
PREFIX="LSSTY10_bin"
SUFFIX="nch70_350_1050"
#
PATH_TO_PARAMS_FILE="/data/AMARINS/CMBWLxHI-CODES/scripts/parameters/lsst_desc_parameters.yaml"
SURVEY="LSST"
DATASET="Y10"
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
bin_max=9
for i in $(seq $bin_min $bin_max)
do
    PREFIX_new="${PREFIX}${i}"
    python3 $PATHFILE --verbose $VERBOSE --pathout $PATHOUT --prefix $PREFIX_new --suffix $SUFFIX\
                      --field_1 $FIELD1 --field_2 $FIELD2\
                      --survey_shear $SURVEY --dataset_shear $DATASET --binning_shear $BINNING\
                      --bin_to_use_shear ${i} --path_to_shear_params_file $PATH_TO_PARAMS_FILE\
                      --freq_min $FREQ_MIN --freq_max $FREQ_MAX  --nchannels $NCHANNELS | tee -a $FILEOUT
done
echo $TIMEF | tee -a $FILEOUT


