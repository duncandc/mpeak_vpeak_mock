#!/bin/bash

python stellar_mass_binned_wp_fiducial_manual.py 9.5 10.0 &

wait

python stellar_mass_binned_wp_fiducial_manual.py 10.0 10.5 &

wait

#python stellar_mass_binned_wp_fiducial_manual.py 10.5 11.0 &

