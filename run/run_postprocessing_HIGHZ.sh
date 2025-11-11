#!/bin/bash
NAME_RUN="HIGHZ_256_nch70"
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
PROJECT='highz_nch70_350_1050'
FILEPATH_FIELD1='/data/AMARINS/CMBWLxHI-CODES/theoretical/highz_HI_cl_nch70_350_1050.txt'
FILEPATH_FIELD2='/data/AMARINS/CMBWLxHI-CODES/theoretical/highz_CMBWL_cl_nch70_350_1050.txt'
FILEPATH_CROSS='/data/AMARINS/CMBWLxHI-CODES/theoretical/highz_CMBWLxHI_cl_nch70_350_1050.txt'
#
DIRPATH_SIMS='/data/AMARINS/CMBWLxHI-DATA/simulations/highz_nch70_350_1050'
DIRPATH_ESTIMATED='/data/AMARINS/CMBWLxHI-DATA/FGremoval/highz_nch70_350_1050/fullsky'
DIRPATH_FOREGROUNDS='/data/AMARINS/CMBWLxHI-DATA/MAPS/FG256'
FILENAME_FOREGROUNDS='FG_I_256_350mhz1050mhz_70bins_full_nonfrps_L0.fits'
DIRPATH_OUT='/data/AMARINS/CMBWLxHI-DATA/postprocessed/fullsky/highz_nch70_350_1050'
#
NREALIZATIONS=100
#PREFIX
#NS=3
##############
RUN_NAME='generating_APS_PostProcessing'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"

#############
eval "$(conda shell.bash hook)"
conda activate amarins_camb

############
cd $DIR_SCRIPTS
#for i in 3 4 5 
for i in 2
do
    python3 $PATHFILE --verbose $VERBOSE --project $PROJECT --dirpath_out $DIRPATH_OUT \
                      --filepath_field1 $FILEPATH_FIELD1 --filepath_field2 $FILEPATH_FIELD2 --filepath_cross $FILEPATH_CROSS \
                      --dirpath_sims $DIRPATH_SIMS --dirpath_estimated $DIRPATH_ESTIMATED\
                      --dirpath_foregrounds $DIRPATH_FOREGROUNDS --filename_foregrounds $FILENAME_FOREGROUNDS \
                      --ns $i --nrealizations $NREALIZATIONS | tee -a $FILEOUT
done

TIMEF=$(date +%Y-%m-%d-%H:%M:%S)
echo $TIMEF | tee -a $FILEOUT