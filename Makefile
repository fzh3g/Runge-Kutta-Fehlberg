# dirs
DIR_BIN = ./bin
DIR_INC = ./include
DIR_SRC = ./src
DIR_DAT = ./data
DIR_IMG = ./img

# The default compiler is g++
CC = g++

# Flags for the compiler. Ask for warnings.
CFLAGS = -Wall -std=c++14 -I${DIR_INC}

# shell
SHELL = /bin/bash

all: ${DIR_BIN}/twobody ${DIR_BIN}/central_config


${DIR_BIN}/twobody: ${DIR_SRC}/twobody.cpp \
		${DIR_INC}/rkf78.hpp ${DIR_INC}/orbit_ellipse_2d.hpp
	$(CC) $(CFLAGS) $< -o $@

${DIR_BIN}/central_config: ${DIR_SRC}/central_config.cpp \
		${DIR_INC}/rkf78.hpp
	$(CC) $(CFLAGS) $< -o $@

.PHONY: all clean run_twobody run_centconf plot_twobody plot_centconf

# Remove files compiled
clean:
	rm *.o

# Run the program
run_twobody:
	cd ${DIR_DAT} && .${DIR_BIN}/twobody && cd ..

run_centconf:
	cd ${DIR_DAT} && .${DIR_BIN}/central_config && cd ..

# Plot orbit trace using python
plot_twobody:
	cd ${DIR_IMG} && .${DIR_SRC}/simple_plot.py \
	.${DIR_DAT}/twobody_output.dat --show \
	--columns 2 3 --equal --del-header 1 --title 'Orbit Trace' \
	--figname 'orbit_trace' --figtype 'png' && cd ..

plot_centconf:
	cd ${DIR_IMG} && .${DIR_SRC}/central_config_plot.sh && cd ..
