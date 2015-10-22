#!/bin/bash

rsync -avz /Users/duncan/Documents/projects/output/processed_data/campbell_mocks/sm_9.5* astro_ftp:~/assembly_bias_mocks/
rsync -avz /Users/duncan/Documents/projects/output/processed_data/campbell_mocks/README.md astro_ftp:~/assembly_bias_mocks/
rsync -avz ./example_open_mock.py astro_ftp:~/assembly_bias_mocks/