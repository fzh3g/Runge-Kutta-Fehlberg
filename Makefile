# dirs
DIR_BIN = ./bin
DIR_INC = ./include
DIR_SRC = ./src
DIR_DAT = ./data
DIR_IMG = ./img

# The default compiler is g++
CC = g++

# Flags for the compiler.
CFLAGS = -Wall -O2 -march=native -pipe -std=c++11 -I${DIR_INC}

# shell
SHELL = /bin/bash

all: ${DIR_BIN}/twobody ${DIR_BIN}/central_config ${DIR_BIN}/pcr3b \
		${DIR_BIN}/poincare_section ${DIR_BIN}/sitnikov \
		${DIR_BIN}/lorenz ${DIR_BIN}/gaussperturb ${DIR_BIN}/sor


${DIR_BIN}/twobody: ${DIR_SRC}/twobody.cpp \
		${DIR_INC}/rkf78.hpp ${DIR_INC}/orbit_ellipse_2d.hpp
	$(CC) $(CFLAGS) $< -o $@

${DIR_BIN}/central_config: ${DIR_SRC}/central_config.cpp \
		${DIR_INC}/rkf78.hpp
	$(CC) $(CFLAGS) $< -o $@

${DIR_BIN}/gaussperturb: ${DIR_SRC}/gaussperturb.cpp \
		${DIR_INC}/rkf78.hpp
	$(CC) $(CFLAGS) $< -o $@

${DIR_BIN}/pcr3b: ${DIR_SRC}/pcr3b.cpp ${DIR_INC}/rkf78.hpp
	$(CC) $(CFLAGS) $< -o $@

${DIR_BIN}/poincare_section: ${DIR_SRC}/poincare_section.cpp \
		${DIR_INC}/rkf78.hpp
	$(CC) $(CFLAGS) $< -o $@

${DIR_BIN}/sitnikov: ${DIR_SRC}/sitnikov.cpp \
		${DIR_INC}/rkf78.hpp
	$(CC) $(CFLAGS) $< -o $@

${DIR_BIN}/lorenz: ${DIR_SRC}/lorenz.cpp \
		${DIR_INC}/rkf78.hpp
	$(CC) $(CFLAGS) $< -o $@

${DIR_BIN}/sor: ${DIR_SRC}/spin_orbit_resonance.cpp \
		${DIR_INC}/rkf78.hpp
	$(CC) $(CFLAGS) $< -o $@

.PHONY: all clean run_twobody plot_twobody run_centconf plot_centconf \
		run_pcr3b plot_pcr3b run_poincsec run_sitnikov plot_sitnikov \
		run_lorenz run_sor plot_sor

# Remove files compiled
clean:
	rm *.o

# Run the program
run_twobody:
	cd ${DIR_DAT} && .${DIR_BIN}/twobody && cd ..

run_centconf:
	cd ${DIR_DAT} && .${DIR_BIN}/central_config && cd ..

run_gaussperturb:
	cd ${DIR_DAT} && .${DIR_BIN}/gaussperturb && cd ..

run_pcr3b:
	cd ${DIR_DAT} && .${DIR_BIN}/pcr3b && cd ..

run_poincsec:
	cd ${DIR_DAT} && .${DIR_BIN}/poincare_section && cd ..

run_sitnikov:
	cd ${DIR_DAT} && .${DIR_BIN}/sitnikov && cd ..

run_lorenz:
	cd ${DIR_DAT} && .${DIR_BIN}/lorenz && cd ..

run_sor:
	cd ${DIR_DAT} && .${DIR_BIN}/sor && cd ..

# Plot orbit trace using python
plot_twobody:
	cd ${DIR_IMG} && .${DIR_SRC}/simple_plot.py \
	.${DIR_DAT}/twobody_output.dat --show \
	--columns 2 3 --equal --del-header 1 --title 'Orbit Trace' \
	--figname 'orbit_trace' --figtype 'png' && cd ..

plot_centconf:
	cd ${DIR_IMG} && .${DIR_SRC}/central_config_plot.sh && cd ..

plot_pcr3b:
	cd ${DIR_IMG} && .${DIR_SRC}/simple_plot.py \
	.${DIR_DAT}/pcr3b.dat --show --sci \
	--columns 2 3 --equal --title 'Planar Circular Restricted 3 Body' \
	--figname 'pcr3b' --figtype 'png' --line && cd ..

plot_gaussperturb:
	cd ${DIR_IMG} && .${DIR_SRC}/gaussperturb.py && cd ..

plot_sitnikov:
	cd ${DIR_IMG} && .${DIR_SRC}/simple_plot.py \
	.${DIR_DAT}/sitnikov.dat --show --sci --del-header 1 \
	--columns 2 3 --title 'Phase Diagram of Sitnikov Problem' \
	--labels 'x' 'dx/dt' --figname 'sitnikov' --figtype 'png' \
	&& cd ..

plot_sor:
	cd ${DIR_IMG} && .${DIR_SRC}/simple_plot.py \
	.${DIR_DAT}/spin_orbit_resonance.dat --show --tex \
	--labels '\theta{}\,(2\pi)' '\dot{\theta}' --xlim 0 1 --ylim -0.5 2.5 \
	--figname 'sor' --figtype 'png'
