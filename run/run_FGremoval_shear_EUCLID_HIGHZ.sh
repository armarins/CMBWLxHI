#!/bin/bash
NAME_RUN="EUCLID_HIGHZ_256_fgremoval"
DIR_OUTSCREEN='/data/AMARINS/CMBWLxHI-CODES/screen_outputs'
DIR_SCRIPTS='/data/AMARINS/CMBWLxHI-CODES/scripts'

FILEOUT="$DIR_OUTSCREEN/$NAME_RUN"
TIMEI=$(date +%Y-%m-%d-%H:%M:%S)
FILEOUT="$FILEOUT.$TIMEI"
FILEOUT="$FILEOUT.out"

echo $TIMEI | tee $FILEOUT

#######
#TERMINAL INFO
FIELD1="HI"
FIELD2="SHEAR"
VERBOSE=1
PROJECT='EUCLID_nch70_350_1050'
DIRPATH_SAVEDATA='/data/AMARINS/LSSxHI-DATA/FGremoval/GWLxHI'
DIRPATH_FOREGROUNDS='/data/AMARINS/MAPS/FG256'
FILENAME_FOREGROUNDS='FG_I_256_350mhz1050mhz_70bins_full_L0.fits'
#APPLY_MASK=1
#DIRPATH_MASK='/data/AMARINS/CMBWLxHI-DATA/MAPS/MASK'
#FILENAME_MASK=
#METHOD='ICA'
#WTRANSFORM='identity'
NSIMS=40
NSIDE=256
NUMIN=350
NUMAX=1050
NBANDS=70
##############
RUN_NAME='generating_APS_FGremoval'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"
##############
eval "$(conda shell.bash hook)"
conda activate /data/AMARINS/anaconda3/amarins
############
bin_min=0
bin_max=9
ns_min=3
ns_max=4
for i in $(seq $bin_min $bin_max); do
    for j in $(seq $ns_min $ns_max); do
        PROJECT_BIN="${PROJECT}/bin${i}"
        FILEPATH_FIELD1="/data/AMARINS/LSSxHI-DATA/theoretical/GWLxHI/EUCLID_bin${i}_${FIELD1}_cl_nch70_350_1050.txt"
        FILEPATH_FIELD2="/data/AMARINS/LSSxHI-DATA/theoretical/GWLxHI/EUCLID_bin${i}_${FIELD2}_cl_nch70_350_1050.txt"  
        FILEPATH_CROSS="/data/AMARINS/LSSxHI-DATA/theoretical/GWLxHI/EUCLID_bin${i}_${FIELD1}x${FIELD2}_cl_nch70_350_1050.txt"
        DIRPATH_SIMS="/data/AMARINS/LSSxHI-DATA/simulations/GWLxHI/EUCLID_nch70_350_1050/bin${i}"
        python $PATHFILE --verbose $VERBOSE --project $PROJECT_BIN --dirpath_sims $DIRPATH_SIMS \
                         --filepath_field1 $FILEPATH_FIELD1 --filepath_field2 $FILEPATH_FIELD2 --filepath_cross $FILEPATH_CROSS \
                         --dirpath_savedata $DIRPATH_SAVEDATA \
                         --dirpath_foregrounds $DIRPATH_FOREGROUNDS --filename_foregrounds $FILENAME_FOREGROUNDS \
                         --nside $NSIDE --numin $NUMIN --numax $NUMAX --nbands $NBANDS \
                         --ns ${j} --nrealizations $NSIMS | tee -a $FILEOUT

    done
done

TIMEF=$(date +%Y-%m-%d-%H:%M:%S)
echo $TIMEF | tee -a $FILEOUT


