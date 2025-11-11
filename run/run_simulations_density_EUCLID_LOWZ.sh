#!/bin/bash
NAME_RUN="EUCLID_nch30_980_1260"
DIR_OUTSCREEN='/data/AMARINS/LSSxHI-CODES/screen_outputs'
DIR_SCRIPTS='/data/AMARINS/LSSxHI-CODES/scripts'

FILEOUT="$DIR_OUTSCREEN/$NAME_RUN"
TIMEI=$(date +%Y-%m-%d-%H:%M:%S)
FILEOUT="$FILEOUT.$TIMEI"
FILEOUT="$FILEOUT.out"

echo $TIMEI | tee $FILEOUT

#######
#TERMINAL INFO
FIELD1="HI"
FIELD2="DENSITY"
VERBOSE=1
PROJECT='EUCLID_nch30_980_1260'
#
PATHOUT="/data/AMARINS/LSSxHI-DATA/simulations/GCxHI"
#
NREALIZATIONS=2
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
##############
eval "$(conda shell.bash hook)"
conda activate /data/AMARINS/anaconda3/amarins
############
bin_min=0
bin_max=9
for i in $(seq $bin_min $bin_max)
do
    PROJECT_BIN="${PROJECT}/bin${i}"
    FILEPATH_FIELD1="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/EUCLID_bin${i}_${FIELD1}_cl_nch30_980_1260.txt"
    FILEPATH_FIELD2="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/EUCLID_bin${i}_${FIELD2}_cl_nch30_980_1260.txt"  
    FILEPATH_CROSS="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/EUCLID_bin${i}_${FIELD1}x${FIELD2}_cl_nch30_980_1260.txt"
    python3 $PATHFILE --verbose $VERBOSE --project $PROJECT_BIN --pathout $PATHOUT \
                      --filepath_field1 $FILEPATH_FIELD1 --filepath_field2 $FILEPATH_FIELD2 --filepath_cross $FILEPATH_CROSS \
                      --limited_correlated_channels $LIMITED_CORRELATED_CHANNELS  --channel_min_corr $CHANNEL_MIN_CORR --channel_max_corr $CHANNEL_MAX_CORR\
                      --seed0 $SEED0 --nrealizations $NREALIZATIONS --channel_tax $CHANNEL_TAX | tee -a $FILEOUT    
done