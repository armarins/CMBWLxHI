#!/bin/bash
NAME_RUN="LOWZ_256_nch30_fgremoval"
DIR_OUTSCREEN='/data/AMARINS/CMBWLxHI-CODES/screen_outputs'
DIR_SCRIPTS='/data/AMARINS/CMBWLxHI-CODES/scripts'

FILEOUT="$DIR_OUTSCREEN/$NAME_RUN"
TIMEI=$(date +%Y-%m-%d-%H:%M:%S)
FILEOUT="$FILEOUT.$TIMEI"
FILEOUT="$FILEOUT.out"

echo $TIMEI | tee $FILEOUT

#######
#TERMINAL INFO
PROJECT='lowz_nch30_980_1260'
DIRPATH_SAVEDATA='/data/AMARINS/CMBWLxHI-DATA/FGremoval'
DIRPATH_SIMS='/data/AMARINS/CMBWLxHI-DATA/simulations/lowz_nch30_980_1260'
FILEPATH_FIELD1='/data/AMARINS/CMBWLxHI-CODES/theoretical/lowz_HI_cl_nch30_980_1260.txt'
FILEPATH_FIELD2='/data/AMARINS/CMBWLxHI-CODES/theoretical/lowz_CMBWL_cl_nch30_980_1260.txt'
FILEPATH_CROSS='/data/AMARINS/CMBWLxHI-CODES/theoretical/lowz_CMBWLxHI_cl_nch30_980_1260.txt'
DIRPATH_FOREGROUNDS='/data/AMARINS/CMBWLxHI-DATA/MAPS/FG256'
FILENAME_FOREGROUNDS='FG_I_256_980mhz1260mhz_30bins_full_nonfrps_L0.fits'
#APPLY_MASK=1
#DIRPATH_MASK='/data/AMARINS/CMBWLxHI-DATA/MAPS/MASK'
#FILENAME_MASK=
#METHOD='ICA'
#WTRANSFORM='identity'
VERBOSE=1
NSIMS=100
NSIDE=256
NUMIN=980
NUMAX=1260
NBANDS=30
##############
RUN_NAME='generating_APS_FGremoval'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"

##############
eval "$(conda shell.bash hook)"
conda activate amarins

############
#for i in 3 4 5 
for i in 2
do
    python $PATHFILE --verbose $VERBOSE --project $PROJECT --dirpath_sims $DIRPATH_SIMS \
                     --filepath_field1 $FILEPATH_FIELD1 --filepath_field2 $FILEPATH_FIELD2 --filepath_cross $FILEPATH_CROSS \
                     --dirpath_savedata $DIRPATH_SAVEDATA \
                     --dirpath_foregrounds $DIRPATH_FOREGROUNDS --filename_foregrounds $FILENAME_FOREGROUNDS \
                     --nside $NSIDE --numin $NUMIN --numax $NUMAX --nbands $NBANDS \
                     --ns $i --nrealizations $NSIMS | tee -a $FILEOUT

done

TIMEF=$(date +%Y-%m-%d-%H:%M:%S)
echo $TIMEF | tee -a $FILEOUT
