#!/bin/bash
NAME_RUN="LSSTY10_nch70_350_1050"
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
FIELD1="HI"
FIELD2="COUNT"
PROJECT='LSSTxHIGHZ_40_0001'
PATHOUT="/data/AMARINS/LSSxHI-DATA/simulations/GCxHI"
SUFFIX="HIGHZ"
#
NREALIZATIONS=100
SEED0=9000
LIMITED_CORRELATED_CHANNELS='no'
CHANNEL_TAX=1010
CHANNEL_MIN_CORR=1
CHANNEL_MAX_CORR=70
jTOMO=40
SIGMAZ="0001"
##############
RUN_NAME='generating_APS_simulations'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"
#############
eval "$(conda shell.bash hook)"
conda activate /data/AMARINS/anaconda3/PY10
############
bin_min=0
bin_max=39
for i in $(seq $bin_min $bin_max)
do
    PROJECT_BIN="${PROJECT}/bin${i}"
    #/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin4_40tomos_0.001sigmaz_HI_cl_LOWZ.txt
    FILEPATH_FIELD1="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin${i}_${jTOMO}_${SIGMAZ}_${FIELD1}_cl_${SUFFIX}.txt"
    FILEPATH_FIELD2="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin${i}_${jTOMO}_${SIGMAZ}_${FIELD2}_cl_${SUFFIX}.txt"  
    FILEPATH_CROSS="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin${i}_${jTOMO}_${SIGMAZ}_${FIELD1}x${FIELD2}_cl_${SUFFIX}.txt"
    
    python3 $PATHFILE --verbose $VERBOSE --project $PROJECT_BIN --pathout $PATHOUT \
                      --filepath_field1 $FILEPATH_FIELD1 --filepath_field2 $FILEPATH_FIELD2 --filepath_cross $FILEPATH_CROSS \
                      #--limited_correlated_channels $LIMITED_CORRELATED_CHANNELS\
                      --channel_min_corr $CHANNEL_MIN_CORR --channel_max_corr $CHANNEL_MAX_CORR\
                      --seed0 $SEED0 --nrealizations $NREALIZATIONS --channel_tax $CHANNEL_TAX | tee -a $FILEOUT    
done