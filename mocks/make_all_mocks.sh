#!/bin/bash

#: <<'END' #use this sequence to comment out a section
##########################################################################################
#Bolshoi mocks
##########################################################################################
#sigma=0.3
python make_mock.py -1.0 9.5 0.3 Bolshoi &
wait
python make_mock.py -0.75 9.5 0.3 Bolshoi &
wait
python make_mock.py -0.5 9.5 0.3 Bolshoi &
wait
python make_mock.py -0.25 9.5 0.3 Bolshoi &
wait
python make_mock.py 0.0 9.5 0.3 Bolshoi &
wait

#sigma=0.2
python make_mock.py -1.0 9.5 0.2 Bolshoi &
wait
python make_mock.py -0.75 9.5 0.2 Bolshoi &
wait
python make_mock.py -0.5 9.5 0.2 Bolshoi &
wait
python make_mock.py -0.25 9.5 0.2 Bolshoi &
wait
python make_mock.py 0.0 9.5 0.2 Bolshoi &
wait

#sigma=0.1
python make_mock.py -1.0 9.5 0.1 Bolshoi &
wait
python make_mock.py -0.75 9.5 0.1 Bolshoi &
wait
python make_mock.py -0.5 9.5 0.1 Bolshoi &
wait
python make_mock.py -0.25 9.5 0.1 Bolshoi &
wait
python make_mock.py 0.0 9.5 0.1 Bolshoi &
wait

#sigma=0.0
python make_mock.py -1.0 9.5 0.0 Bolshoi &
wait
python make_mock.py -0.75 9.5 0.0 Bolshoi &
wait
python make_mock.py -0.5 9.5 0.0 Bolshoi &
wait
python make_mock.py -0.25 9.5 0.0 Bolshoi &
wait
python make_mock.py 0.0 9.5 0.0 Bolshoi &
wait

#END

##########################################################################################
#Chinchilla mocks
##########################################################################################
#sigma=0.3
python make_mock.py -1.0 8.5 0.3 Chinchilla &
wait
python make_mock.py -0.75 8.5 0.3 Chinchilla &
wait
python make_mock.py -0.5 8.5 0.3 Chinchilla &
wait
python make_mock.py -0.25 8.5 0.3 Chinchilla &
wait
python make_mock.py 0.0 8.5 0.3 Chinchilla &

#sigma=0.2
python make_mock.py -1.0 8.5 0.2 Chinchilla &
wait
python make_mock.py -0.75 8.5 0.2 Chinchilla &
wait
python make_mock.py -0.5 8.5 0.2 Chinchilla &
wait
python make_mock.py -0.25 8.5 0.2 Chinchilla &
wait
python make_mock.py 0.0 8.5 0.2 Chinchilla &

#sigma=0.1
python make_mock.py -1.0 8.5 0.1 Chinchilla &
wait
python make_mock.py -0.75 8.5 0.1 Chinchilla &
wait
python make_mock.py -0.5 8.5 0.1 Chinchilla &
wait
python make_mock.py -0.25 8.5 0.1 Chinchilla &
wait
python make_mock.py 0.0 8.5 0.1 Chinchilla &

#sigma=0.0
python make_mock.py -1.0 8.5 0.0 Chinchilla &
wait
python make_mock.py -0.75 8.5 0.0 Chinchilla &
wait
python make_mock.py -0.5 8.5 0.0 Chinchilla &
wait
python make_mock.py -0.25 8.5 0.0 Chinchilla &
wait
python make_mock.py 0.0 8.5 0.0 Chinchilla &