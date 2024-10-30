#!/bin/bash
NAME_RUN="LOWZ_256_nch30"
DIR_OUTSCREEN='/data/AMARINS/CMBWLxHI-CODES/screen_outputs'
DIR_SCRIPTS='/data/AMARINS/CMBWLxHI-CODES/scripts'

FILEOUT="$DIR_OUTSCREEN/$NAME_RUN"
TIMEI=$(date +%Y-%m-%d-%H:%M:%S)
FILEOUT="$FILEOUT.$TIMEI"
FILEOUT="$FILEOUT.out"

echo $TIMEI | tee $FILEOUT

#######
#TERMINAL INFO
VERBOSE=1
PROJECT='lowz_nch30_980_1260'
FILEPATH_FIELD1='/data/AMARINS/CMBWLxHI-CODES/theoretical/lowz_HI_cl_nch30_980_1260.txt'
FILEPATH_FIELD2='/data/AMARINS/CMBWLxHI-CODES/theoretical/lowz_CMBWL_cl_nch30_980_1260.txt'
FILEPATH_CROSS='/data/AMARINS/CMBWLxHI-CODES/theoretical/lowz_CMBWLxHI_cl_nch30_980_1260.txt'
#
PATHOUT='/data/AMARINS/CMBWLxHI-DATA/simulations/'
#
NREALIZATIONS=100
SEED0=9000
LIMITED_CORRELATED_CHANNELS=1
CHANNEL_TAX=1010
CHANNEL_MIN_CORR=1
CHANNEL_MAX_CORR=30

##############
RUN_NAME='generating_APS_simulations'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"

#############
eval "$(conda shell.bash hook)"
conda activate amarins_camb

############
cd $DIR_SCRIPTS
python3 $PATHFILE --verbose $VERBOSE --project $PROJECT --pathout $PATHOUT \
                 --filepath_field1 $FILEPATH_FIELD1 --filepath_field2 $FILEPATH_FIELD2 --filepath_cross $FILEPATH_CROSS \
                 --limited_correlated_channels $LIMITED_CORRELATED_CHANNELS  --channel_min_corr $CHANNEL_MIN_CORR --channel_max_corr $CHANNEL_MAX_CORR\
                 --seed0 $SEED0 --nrealizations $NREALIZATIONS --channel_tax $CHANNEL_TAX | tee -a $FILEOUT


TIMEF=$(date +%Y-%m-%d-%H:%M:%S)
echo $TIMEF | tee -a $FILEOUT



