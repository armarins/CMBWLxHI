#!/bin/bash
NAME_RUN="LSST_LOWZ_DENSITY_256_postprocessing"
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
FIELD1="HI"
FIELD2="DENSITY"
PROJECT='LSSTY10_nch30_980_1260'
#
DIRPATH_FOREGROUNDS='/data/AMARINS/MAPS/FG256'
FILENAME_FOREGROUNDS='FG_I_256_980mhz1260mhz_30bins_full_nonfrps_L0.fits'
#
NREALIZATIONS=40
#PREFIX
#NS=3
##############
RUN_NAME='generating_APS_PostProcessing'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"
##############
eval "$(conda shell.bash hook)"
conda activate /data/AMARINS/anaconda3/amarins
############
bin_min=0
bin_max=4
ns_min=3
ns_max=4
for i in $(seq $bin_min $bin_max); do
    for j in $(seq $ns_min $ns_max); do
        PROJECT_BIN="${PROJECT}/bin${i}"
        FILEPATH_FIELD1="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSSTY10_bin${i}_${FIELD1}_cl_nch30_980_1260.txt"
        FILEPATH_FIELD2="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSSTY10_bin${i}_${FIELD2}_cl_nch30_980_1260.txt"  
        FILEPATH_CROSS="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSSTY10_bin${i}_${FIELD1}x${FIELD2}_cl_nch30_980_1260.txt"
        DIRPATH_SIMS="/data/AMARINS/LSSxHI-DATA/simulations/GCxHI/LSSTY10_nch30_980_1260/bin${i}"        
        DIRPATH_ESTIMATED="/data/AMARINS/LSSxHI-DATA/FGremoval/GCxHI/LSSTY10_nch30_980_1260/bin${i}/fullsky"
        DIRPATH_OUT="/data/AMARINS/LSSxHI-DATA/postprocessed/GCxHI/LSSTY10_nch30_980_1260/bin${i}/fullsky"
        python3 $PATHFILE --verbose $VERBOSE --project $PROJECT_BIN --dirpath_out $DIRPATH_OUT \
                          --filepath_field1 $FILEPATH_FIELD1 --filepath_field2 $FILEPATH_FIELD2 --filepath_cross $FILEPATH_CROSS \
                          --dirpath_sims $DIRPATH_SIMS --dirpath_estimated $DIRPATH_ESTIMATED\
                          --dirpath_foregrounds $DIRPATH_FOREGROUNDS --filename_foregrounds $FILENAME_FOREGROUNDS \
                          --ns ${j} --nrealizations $NREALIZATIONS | tee -a $FILEOUT
    done
done

TIMEF=$(date +%Y-%m-%d-%H:%M:%S)
echo $TIMEF | tee -a $FILEOUT