#!/bin/bash

#: <<'END' #use this sequence to comment out a section
#END

##########################################################################################
#Bolshoi
##########################################################################################

#0.3 scatter
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-1.0_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-1.0_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-1.0_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-1.0_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-0.75_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-0.75_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-0.75_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-0.75_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-0.5_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-0.5_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-0.5_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-0.5_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-0.25_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-0.25_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-0.25_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c-0.25_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c0.0_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c0.0_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c0.0_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.3_sfr_c0.0_250' 11.0 11.5 &
wait

#0.2 scatter
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-1.0_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-1.0_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-1.0_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-1.0_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-0.75_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-0.75_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-0.75_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-0.75_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-0.5_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-0.5_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-0.5_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-0.5_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-0.25_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-0.25_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-0.25_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c-0.25_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c0.0_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c0.0_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c0.0_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.2_sfr_c0.0_250' 11.0 11.5 &
wait

#0.1 scatter
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-1.0_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-1.0_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-1.0_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-1.0_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-0.75_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-0.75_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-0.75_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-0.75_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-0.5_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-0.5_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-0.5_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-0.5_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-0.25_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-0.25_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-0.25_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c-0.25_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c0.0_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c0.0_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c0.0_250' 10.5 11.0 &
wait 
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.1_sfr_c0.0_250' 11.0 11.5 &
wait 

#0.0 scatter
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-1.0_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-1.0_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-1.0_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-1.0_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-0.75_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-0.75_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-0.75_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-0.75_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-0.5_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-0.5_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-0.5_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-0.5_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-0.25_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-0.25_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-0.25_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c-0.25_250' 11.0 11.5 &
wait

python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c0.0_250' 9.5 10.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c0.0_250' 10.0 10.5 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c0.0_250' 10.5 11.0 &
wait
python stellar_mass_binned_central_wp.py 'sm_9.5_s0.0_sfr_c0.0_250' 11.0 11.5 &
wait

<<COMMENT
##########################################################################################
#Chinchilla
##########################################################################################

#0.3 scatter
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-1.0_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-1.0_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-1.0_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-1.0_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-0.75_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-0.75_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-0.75_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-0.75_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-0.5_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-0.5_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-0.5_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-0.5_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-0.25_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-0.25_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-0.25_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c-0.25_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c0.0_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c0.0_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c0.0_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.3_sfr_c0.0_125' 10.5 11.0 &
wait

#0.2 scatter
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-1.0_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-1.0_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-1.0_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-1.0_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-0.75_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-0.75_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-0.75_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-0.75_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-0.5_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-0.5_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-0.5_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-0.5_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-0.25_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-0.25_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-0.25_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c-0.25_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c0.0_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c0.0_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c0.0_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.2_sfr_c0.0_125' 10.5 11.0 &
wait

#0.1 scatter
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-1.0_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-1.0_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-1.0_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-1.0_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-0.75_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-0.75_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-0.75_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-0.75_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-0.5_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-0.5_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-0.5_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-0.5_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-0.25_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-0.25_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-0.25_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c-0.25_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c0.0_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c0.0_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c0.0_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.1_sfr_c0.0_125' 10.5 11.0 &
wait

#0.0 scatter
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-1.0_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-1.0_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-1.0_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-1.0_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-0.75_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-0.75_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-0.75_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-0.75_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-0.5_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-0.5_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-0.5_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-0.5_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-0.25_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-0.25_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-0.25_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c-0.25_125' 10.5 11.0 &
wait

python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c0.0_125' 9.0 9.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c0.0_125' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c0.0_125' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_8.5_s0.0_sfr_c0.0_125' 10.5 11.0 &
wait
COMMENT

