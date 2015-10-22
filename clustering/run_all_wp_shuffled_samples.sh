#!/bin/bash

#: <<'END' #use this sequence to comment out a section
#END

##########################################################################################
#Bolshoi central shuffle
##########################################################################################

#0.3 scatter
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-1.0_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-1.0_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-1.0_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-1.0_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.75_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.75_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.75_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.75_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.5_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.5_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.5_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.5_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.25_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.25_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.25_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.25_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c0.0_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c0.0_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c0.0_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c0.0_250_cen_shuffle' 11.0 11.5 &
wait

#0.2 scatter
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-1.0_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-1.0_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-1.0_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-1.0_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.75_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.75_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.75_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.75_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.5_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.5_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.5_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.5_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.25_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.25_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.25_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.25_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c0.0_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c0.0_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c0.0_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c0.0_250_cen_shuffle' 11.0 11.5 &
wait

#0.1 scatter
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-1.0_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-1.0_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-1.0_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-1.0_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.75_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.75_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.75_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.75_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.5_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.5_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.5_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.5_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.25_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.25_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.25_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.25_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c0.0_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c0.0_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c0.0_250_cen_shuffle' 10.5 11.0 &
wait 
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c0.0_250_cen_shuffle' 11.0 11.5 &
wait 

#0.0 scatter
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-1.0_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-1.0_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-1.0_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-1.0_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.75_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.75_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.75_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.75_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.5_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.5_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.5_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.5_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.25_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.25_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.25_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.25_250_cen_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c0.0_250_cen_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c0.0_250_cen_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c0.0_250_cen_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c0.0_250_cen_shuffle' 11.0 11.5 &
wait


##########################################################################################
#Bolshoi central + satellite shuffle
##########################################################################################

#0.3 scatter
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-1.0_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-1.0_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-1.0_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-1.0_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.75_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.75_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.75_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.75_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.5_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.5_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.5_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.5_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.25_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.25_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.25_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c-0.25_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c0.0_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c0.0_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c0.0_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.3_sfr_c0.0_250_cen_sat_shuffle' 11.0 11.5 &
wait

#0.2 scatter
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-1.0_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-1.0_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-1.0_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-1.0_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.75_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.75_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.75_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.75_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.5_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.5_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.5_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.5_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.25_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.25_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.25_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c-0.25_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c0.0_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c0.0_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c0.0_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.2_sfr_c0.0_250_cen_sat_shuffle' 11.0 11.5 &
wait

#0.1 scatter
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-1.0_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-1.0_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-1.0_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-1.0_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.75_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.75_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.75_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.75_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.5_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.5_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.5_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.5_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.25_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.25_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.25_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c-0.25_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c0.0_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c0.0_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c0.0_250_cen_sat_shuffle' 10.5 11.0 &
wait 
python stellar_mass_binned_wp.py 'sm_9.5_s0.1_sfr_c0.0_250_cen_sat_shuffle' 11.0 11.5 &
wait 

#0.0 scatter
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-1.0_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-1.0_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-1.0_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-1.0_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.75_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.75_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.75_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.75_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.5_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.5_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.5_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.5_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.25_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.25_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.25_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c-0.25_250_cen_sat_shuffle' 11.0 11.5 &
wait

python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c0.0_250_cen_sat_shuffle' 9.5 10.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c0.0_250_cen_sat_shuffle' 10.0 10.5 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c0.0_250_cen_sat_shuffle' 10.5 11.0 &
wait
python stellar_mass_binned_wp.py 'sm_9.5_s0.0_sfr_c0.0_250_cen_sat_shuffle' 11.0 11.5 &
wait
