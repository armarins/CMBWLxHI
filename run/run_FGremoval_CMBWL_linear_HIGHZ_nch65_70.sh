#!/bin/bash
NAME_RUN="CMB_linear_HIGHZ_256_nch70_fgremoval_65_70"
DIR_OUTSCREEN='/data/AMARINS/LSSxHI-CODES/screen_outputs'
DIR_SCRIPTS='/data/AMARINS/LSSxHI-CODES/scripts'

FILEOUT="$DIR_OUTSCREEN/$NAME_RUN"
TIMEI=$(date +%Y-%m-%d-%H:%M:%S)
FILEOUT="$FILEOUT.$TIMEI"
FILEOUT="$FILEOUT.out"

echo $TIMEI | tee $FILEOUT


#######
#TERMINAL INFO
PROJECT='CMB_linear_HIGHZ_nch70_350_1050_65_70'
DIRPATH_SAVEDATA='/data/AMARINS/LSSxHI-DATA/FGremoval/CMBWLxHI'
DIRPATH_SIMS='/data/AMARINS/LSSxHI-DATA/simulations/CMBWLxHI/CMB_linear_HIGHZ_nch70_350_1050_65_70'
FILEPATH_FIELD1='/data/AMARINS/LSSxHI-DATA/theoretical/CMBWLxHI/CMBWL_linear_HI_cl_nch70_350_1050.txt'
FILEPATH_FIELD2='/data/AMARINS/LSSxHI-DATA/theoretical/CMBWLxHI/CMBWL_linear_CMBWL_cl_nch70_350_1050.txt'
FILEPATH_CROSS='/data/AMARINS/LSSxHI-DATA/theoretical/CMBWLxHI/CMBWL_linear_HIxCMBWL_cl_nch70_350_1050.txt'
DIRPATH_FOREGROUNDS='/data/AMARINS/MAPS/FG256'
FILENAME_FOREGROUNDS='FG_I_256_350mhz1050mhz_70bins_full_L0.fits'
DIRPATH_CHISEL='/data/AMARINS/scripts'
#APPLY_MASK=1
#DIRPATH_MASK='/data/AMARINS/CMBWLxHI-DATA/MAPS/MASK'
#FILENAME_MASK=
#METHOD='ICA'
#WTRANSFORM='identity'
VERBOSE=1
NSIMS=100
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
ns_min=3
ns_max=5
for j in $(seq $ns_min $ns_max);
do
    python $PATHFILE --verbose $VERBOSE --project $PROJECT --dirpath_sims $DIRPATH_SIMS \
                     --filepath_field1 $FILEPATH_FIELD1 --filepath_field2 $FILEPATH_FIELD2 --filepath_cross $FILEPATH_CROSS \
                     --dirpath_savedata $DIRPATH_SAVEDATA \
                     --dirpath_foregrounds $DIRPATH_FOREGROUNDS --filename_foregrounds $FILENAME_FOREGROUNDS \
                     --nside $NSIDE --numin $NUMIN --numax $NUMAX --nbands $NBANDS \
                     --ns ${j} --nrealizations $NSIMS | tee -a $FILEOUT

done

TIMEF=$(date +%Y-%m-%d-%H:%M:%S)
echo $TIMEF | tee -a $FILEOUT


