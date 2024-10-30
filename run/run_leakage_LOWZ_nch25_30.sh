#!/bin/bash
NAME_RUN="LOWZ_256_100_nch25_30"
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
PROJECT='lowz_nch30_980_1260_nch25_30'
FILEPATH_CROSS='/data/AMARINS/CMBWLxHI-CODES/theoretical/lowz_CMBWLxHI_cl_nch30_980_1260.txt'
#
DIRPATH_SIMS='/data/AMARINS/CMBWLxHI-DATA/simulations/lowz_nch30_980_1260_nch25_30'
DIRPATH_ESTIMATED='/data/AMARINS/CMBWLxHI-DATA/FGremoval/lowz_nch30_980_1260_nch25_30/fullsky'
DIRPATH_POSTPROC='/data/AMARINS/CMBWLxHI-DATA/postprocessed/fullsky/lowz_nch30_980_1260_nch25_30'
DIRPATH_OUT='/data/AMARINS/CMBWLxHI-DATA/leakage/fullsky/lowz_nch30_980_1260_nch25_30'

#
NREALIZATIONS=100
NUMIN=24
NUMAX=29
#PREFIX
#NS=3
##############
RUN_NAME='generating_APS_leakage_estimation'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"

#############
eval "$(conda shell.bash hook)"
conda activate amarins_camb

############
#cd $DIRPATH_OUT
#rm -rf *
#mkdir 

############
cd $DIR_SCRIPTS
#for i in 3 4 5
for i in 2
do
    python3 $PATHFILE --verbose $VERBOSE --project $PROJECT --dirpath_out $DIRPATH_OUT \
                      --filepath_cross $FILEPATH_CROSS --dirpath_postprocessing $DIRPATH_POSTPROC\
                      --dirpath_sims $DIRPATH_SIMS --dirpath_estimated $DIRPATH_ESTIMATED\
                      --ns $i --nrealizations $NREALIZATIONS --nu_min_correlated $NUMIN --nu_max_correlated $NUMAX | tee -a $FILEOUT
done

TIMEF=$(date +%Y-%m-%d-%H:%M:%S)
echo $TIMEF | tee -a $FILEOUT