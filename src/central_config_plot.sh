#!/bin/bash

echo "Plotting phi-t figure..."

../src/simple_plot.py ../data/central_config.dat --columns 0 2 \
                      --labels 't' '\phi{}' --tex \
                      --figname 'central_config_phi' \
                      --xlim 0 20 --ylim 0 12 --figtype 'png'

echo "Plotting phidot-t figure..."

../src/simple_plot.py ../data/central_config.dat --columns 0 3 \
                      --labels 't' 'd\phi{}/dt' --tex \
                      --figname 'central_config_phidot' \
                      --xlim 0 20 --ylim -4 1 --figtype 'png'

echo "Plotting phidot-phi figure..."

../src/simple_plot.py ../data/central_config.dat --columns 2 3 \
                      --labels '\phi{}' 'd\phi{}/dt' --tex \
                      --figname 'central_config_phiphase' \
                      --xlim 0 8 --ylim -4 1 --figtype 'png'

echo "Done!"
