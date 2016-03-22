# The default compiler is g++
CC = g++

# Flags for the compiler. Ask for warnings.
CFLAGS = -Wall -std=c++14

twobody: twobody.cpp rkf78.hpp orbit_ellipse_2d.hpp
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean run plot

# Remove files compiled
clean:
	rm twobody

# Run the program
run:
	./twobody

# Plot orbit trace using python
plot:
	./plot_trace.py
