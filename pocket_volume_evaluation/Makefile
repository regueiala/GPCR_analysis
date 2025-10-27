# Makefile for running EPOCK volume computation pipeline

# Phony targets (not actual files)
.PHONY: all convert extract process epock clean check_jobs

# Default target
all: convert extract process epock

# 1. Convert DCD to XTC
convert:
	@echo "=== Step 1: Converting DCD to XTC ==="
	python3 dcd2xtc.py

# 2. Extract geometry with VMD in text mode
extract:
	@echo "=== Step 2: Running VMD extract_geom.tcl (wait 20 seconds) ==="
	@vmd -dispdev text -e extract_geom.tcl > vmd_output.log 2>&1 &
	@sleep 20  # wait for initialization
	@VMD_PID=$$!; \
	echo "VMD started with PID $$VMD_PID"; \
	wait $$VMD_PID || true; \
	if ps -p $$VMD_PID > /dev/null; then \
		echo "Killing VMD process $$VMD_PID"; \
		kill -9 $$VMD_PID; \
	fi

# 3. Process data
process:
	@echo "=== Step 3: Processing .dat files ==="
	python3 treat_dat.py

# 4. Run EPOCK
epock:
	@echo "=== Step 4: Running EPOCK ==="
	./launch_epock.sh

# Optional cleanup
clean:
	rm -f vmd_output.log

# Optional cleanup
check_jobs:
        tsp -C
	tsp
