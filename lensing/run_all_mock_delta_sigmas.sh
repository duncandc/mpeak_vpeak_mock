#!/bin/bash

: <<'END' #use this sequence to comment out a section
##########################################################################################
#threshold samples
##########################################################################################

#sigma = 0.2
python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c-1.0_250' 9.5 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c-1.0_250' 10.0 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c-1.0_250' 10.5 &
wait

python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c-0.75_250' 9.5 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c-0.75_250' 10.0 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c-0.75_250' 10.5 &
wait

python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c-0.5_250' 9.5 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c-0.5_250' 10.0 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c-0.5_250' 10.5 &
wait

python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c-0.25_250' 9.5 &
wait
python gal_gal_lensin_thresholdg.py 'sm_9.5_s0.2_sfr_c-0.25_250' 10.0 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c-0.25_250' 10.5 &
wait

python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c0.0_250' 9.5 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c0.0_250' 10.0 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.2_sfr_c0.0_250' 10.5 &
wait

#sigma = 0.0
python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c-1.0_250' 9.5 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c-1.0_250' 10.0 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c-1.0_250' 10.5 &
wait

python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c-0.75_250' 9.5 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c-0.75_250' 10.0 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c-0.75_250' 10.5 &
wait

python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c-0.5_250' 9.5 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c-0.5_250' 10.0 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c-0.5_250' 10.5 &
wait

python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c-0.25_250' 9.5 &
wait
python gal_gal_lensin_thresholdg.py 'sm_9.5_s0.0_sfr_c-0.25_250' 10.0 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c-0.25_250' 10.5 &
wait

python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c0.0_250' 9.5 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c0.0_250' 10.0 &
wait
python gal_gal_lensing_threshold.py 'sm_9.5_s0.0_sfr_c0.0_250' 10.5 &
wait 

END

##########################################################################################
#binned samples
##########################################################################################

#sigma = 0.3
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-1.0_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-1.0_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-1.0_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-1.0_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-0.75_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-0.75_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-0.75_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-0.75_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-0.5_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-0.5_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-0.5_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-0.5_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-0.25_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-0.25_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-0.25_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c-0.25_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c0.0_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c0.0_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c0.0_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.3_sfr_c0.0_250' 11.0 11.5 &
wait

#sigma = 0.2
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-1.0_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-1.0_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-1.0_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-1.0_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-0.75_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-0.75_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-0.75_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-0.75_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-0.5_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-0.5_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-0.5_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-0.5_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-0.25_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-0.25_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-0.25_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c-0.25_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c0.0_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c0.0_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c0.0_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.2_sfr_c0.0_250' 11.0 11.5 &
wait

#sigma = 0.1
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-1.0_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-1.0_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-1.0_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-1.0_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-0.75_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-0.75_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-0.75_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-0.75_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-0.5_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-0.5_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-0.5_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-0.5_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-0.25_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-0.25_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-0.25_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c-0.25_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c0.0_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c0.0_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c0.0_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.1_sfr_c0.0_250' 11.0 11.5 &
wait

#sigma = 0.0
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-1.0_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-1.0_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-1.0_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-1.0_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-0.75_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-0.75_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-0.75_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-0.75_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-0.5_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-0.5_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-0.5_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-0.5_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-0.25_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-0.25_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-0.25_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c-0.25_250' 11.0 11.5 &
wait

python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c0.0_250' 9.5 10.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c0.0_250' 10.0 10.5 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c0.0_250' 10.5 11.0 &
wait
python gal_gal_lensing_bins.py 'sm_9.5_s0.0_sfr_c0.0_250' 11.0 11.5 &

