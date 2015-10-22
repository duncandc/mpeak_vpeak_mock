#!/bin/bash

##########################################################################################
#Bolshoi mocks
##########################################################################################
#sigma=0.3
python shuffle_mock.py 'sm_9.5_s0.3_sfr_c-1.0_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.3_sfr_c-0.75_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.3_sfr_c-0.5_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.3_sfr_c-0.25_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.3_sfr_c0.0_250' &
wait

#sigma=0.2
python shuffle_mock.py 'sm_9.5_s0.2_sfr_c-1.0_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.2_sfr_c-0.75_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.2_sfr_c-0.5_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.2_sfr_c-0.25_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.2_sfr_c0.0_250' &
wait

#sigma=0.1
python shuffle_mock.py 'sm_9.5_s0.1_sfr_c-1.0_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.1_sfr_c-0.75_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.1_sfr_c-0.5_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.1_sfr_c-0.25_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.1_sfr_c0.0_250' &
wait

#sigma=0.0
python shuffle_mock.py 'sm_9.5_s0.0_sfr_c-1.0_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.0_sfr_c-0.75_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.0_sfr_c-0.5_250'i &
wait
python shuffle_mock.py 'sm_9.5_s0.0_sfr_c-0.25_250' &
wait
python shuffle_mock.py 'sm_9.5_s0.0_sfr_c0.0_250' &
