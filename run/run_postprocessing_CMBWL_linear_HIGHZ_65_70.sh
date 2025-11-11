#!/bin/bash
NAME_RUN="CMB_linear_HIGHZ_256_nch70_postprocessing_65_70"
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
FIELD1="HI"
FIELD2="CMBWL"
PROJECT='CMB_linear_HIGHZ_nch70_350_1050_65_70'
DIRPATH_FOREGROUNDS='/data/AMARINS/MAPS/FG256'
FILENAME_FOREGROUNDS='FG_I_256_350mhz1050mhz_70bins_full_L0.fits'
NREALIZATIONS=100
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
ns_min=3
ns_max=5
for j in $(seq $ns_min $ns_max); do
    FILEPATH_FIELD1="/data/AMARINS/LSSxHI-DATA/theoretical/CMBWLxHI/CMBWL_linear_${FIELD1}_cl_nch70_350_1050.txt"
    FILEPATH_FIELD2="/data/AMARINS/LSSxHI-DATA/theoretical/CMBWLxHI/CMBWL_linear_${FIELD2}_cl_nch70_350_1050.txt"
    FILEPATH_CROSS="/data/AMARINS/LSSxHI-DATA/theoretical/CMBWLxHI/CMBWL_linear_${FIELD1}x${FIELD2}_cl_nch70_350_1050.txt"
    DIRPATH_SIMS="/data/AMARINS/LSSxHI-DATA/simulations/CMBWLxHI/CMB_linear_HIGHZ_nch70_350_1050_65_70/"        
    DIRPATH_ESTIMATED="/data/AMARINS/LSSxHI-DATA/FGremoval/CMBWLxHI/CMB_linear_HIGHZ_nch70_350_1050_65_70/fullsky"
    DIRPATH_OUT="/data/AMARINS/LSSxHI-DATA/postprocessed/CMBWLxHI/CMB_linear_HIGHZ_nch70_350_1050_65_70/fullsky"
    python3 $PATHFILE --verbose $VERBOSE --project $PROJECT --dirpath_out $DIRPATH_OUT \
                          --filepath_field1 $FILEPATH_FIELD1 --filepath_field2 $FILEPATH_FIELD2 --filepath_cross $FILEPATH_CROSS \
                          --dirpath_sims $DIRPATH_SIMS --dirpath_estimated $DIRPATH_ESTIMATED\
                          --dirpath_foregrounds $DIRPATH_FOREGROUNDS --filename_foregrounds $FILENAME_FOREGROUNDS \
                          --ns ${j} --nrealizations $NREALIZATIONS 
                          #| tee -a $FILEOUT
done
TIMEF=$(date +%Y-%m-%d-%H:%M:%S)
echo $TIMEF
#| tee -a $FILEOUT