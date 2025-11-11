#!/bin/bash
NAME_RUN="LSSTY10_nch30_980_1260__40_0001"
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
FIELD2="DENSITY"
PROJECT='LSST_bin4_40_0001'
PATHOUT="/data/AMARINS/LSSxHI-DATA/tfunction/LSST_40_0001"
SUFFIX="LOWZ"
#
NREALIZATIONS=500
INCLUDE_WN='true'
INCLUDE_BEAM='true'
WN_SIGMApix=0.154
BEAM_FWHM=0.011635
##############
TF_HIxHI="true"
TF_GxHI="true"
CROSS_G_BIN=5
AUTO_HI_BINS=-1
CROSS_HI_BINS=-1
##############
SAVE_SIMULATED='true'
SAVE_THEORETICAL='true'
SAVE_HDF5='true'
SAVE_TXT="no"
#############
RUN_NAME='generating_APS_TransferFunction'
INI_FILE="$RUN_NAME.ini"
RUN_FILE="$RUN_NAME.py"
PATHFILE="$DIR_SCRIPTS/$RUN_FILE"
#############
eval "$(conda shell.bash hook)"
conda activate /data/AMARINS/anaconda3/PY10
############
bin_min=1
bin_max=7
for i in $(seq $bin_min $bin_max)
do
    PROJECT_BIN="${PROJECT}/bin${i}"
    #/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin4_40tomos_0.001sigmaz_HI_cl_LOWZ.txt
    FILEPATH_FIELD1="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin${i}_${jTOMO}_${SIGMAZ}_${FIELD1}_cl_${SUFFIX}.txt"
    FILEPATH_FIELD2="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin${i}_${jTOMO}_${SIGMAZ}_${FIELD2}_cl_${SUFFIX}.txt"  
    FILEPATH_CROSS="/data/AMARINS/LSSxHI-DATA/theoretical/GCxHI/LSST_bin${i}_${jTOMO}_${SIGMAZ}_${FIELD1}x${FIELD2}_cl_${SUFFIX}.txt"
    DIRPATH_SIMS="/data/AMARINS/LSSxHI-DATA/simulations/GCxHI/LSSTxLOWZ_40_0001/bin${i}"
    
    python3 $PATHFILE --verbose $VERBOSE --pathout $PATHOUT --nrealizations $NREALIZATIONS\
                      #--inclde_wn $INCLUDE_WN --include_beam $INCLUDE_BEAM \
                      #--wn_sigma_pix $WN_SIGMApix --beam_fwhm $BEAM_FWHM \
                      #--tf_hixhi $TF_HIxHI --tf_gxhi $TF_GxHI \
                      --cross_g_bin $CROSS_G_BIN \ #--auto_hi_bins $AUTO_HI_BINS --cross_hi_bins $CROSS_HI_BINS \
                      --save_simulated $SAVE_SIMULATED --save_theoretical $SAVE_THEORETICAL --save_hdf5 $SAVE_HDF5 --save_txt $SAVE_TXT\
                      --filepath_field1 $FILEPATH_FIELD1 --filepath_field2 $FILEPATH_FIELD2 --filepath_cross $FILEPATH_CROSS \
                      --dirpath_sims $DIRPATH_SIMS | tee -a $FILEOUT    
                      #--seed0 $SEED0 --channel_tax $CHANNEL_TAX | tee -a $FILEOUT    
done