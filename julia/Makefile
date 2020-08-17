SEED = $(word 3, $(shell grep seed 2d.ini))
OUTPUT := output$(SEED)

all: $(OUTPUT)

$(OUTPUT): temperatures.txt Jij.txt utriangles.txt test 2d.ini
	@sh local_job.sh

%.txt: mk_%.jl
	@julia $< 

