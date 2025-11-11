#!/bin/bash
NAME_RUN="CMB_linear_HIGHZ_256_nch70_65_70"
DIR_OUTSCREEN='/data/AMARINS/LSSxHI-CODES/screen_outputs'
DIR_SCRIPTS='/data/AMARINS/LSSxHI-CODES/scripts'

FILEOUT="$DIR_OUTSCREEN/$NAME_RUN"
TIMEI=$(date +%Y-%m-%d-%H:%M:%S)
FILEOUT="$FILEOUT.$TIMEI"
FILEOUT="$FILEOUT.out"

echo $TIMEI 
#| tee $FILEOUT

#######
#TERMINAL INFO
VERBOSE=1
PROJECT='CMB_linear_HIGHZ_nch70_350_1050_65_70'
FILEPATH_FIELD1='/data/AMARINS/LSSxHI-DATA/theoretical/CMBWLxHI/CMBWL_linear_HI_cl_nch70_350_1050.txt'
FILEPATH_FIELD2='/data/AMARINS/LSSxHI-DATA/theoretical/CMBWLxHI/CMBWL_linear_CMBWL_cl_nch70_350_1050.txt'
FILEPATH_CROSS='/data/AMARINS/LSSxHI-DATA/theoretical/CMBWLxHI/CMBWL_linear_HIxCMBWL_cl_nch70_350_1050.txt'
#
PATHOUT='/data/AMARINS/LSSxHI-DATA/simulations/CMBWLxHI/'
#
NREALIZATIONS=100
SEED0=9000
LIMITED_CORRELATED_CHANNELS=1
CHANNEL_TAX=1010
CHANNEL_MIN_CORR=65
CHANNEL_MAX_CORR=70

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
                 --seed0 $SEED0 --nrealizations $NREALIZATIONS --channel_tax $CHANNEL_TAX 
                 #| tee -a $FILEOUT


TIMEF=$(date +%Y-%m-%d-%H:%M:%S)
echo $TIMEF 
#| tee -a $FILEOUT



