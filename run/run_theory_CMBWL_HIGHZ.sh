#!/bin/bash
NAME_RUN="HIGHZ_CMBWL_256_nch70"
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
FIELD2='cmbwl'
PATHOUT='/data/AMARINS/LSSxHI-DATA/theoretical/CMBWLxHI'
PREFIX="CMBWL"
SUFFIX="nch70_350_1050"
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

python3 $PATHFILE --verbose $VERBOSE --pathout $PATHOUT --prefix ${PREFIX} --suffix $SUFFIX\
                  --field_1 $FIELD1 --field_2 $FIELD2\
                  --freq_min $FREQ_MIN --freq_max $FREQ_MAX  --nchannels $NCHANNELS | tee -a $FILEOUT
echo $TIMEF | tee -a $FILEOUT


