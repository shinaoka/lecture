SEED = $(word 3, $(shell grep seed 2d.ini))
OUTPUT := output$(SEED)

all: $(OUTPUT)

$(OUTPUT): temperatures.txt Jij.txt utriangles.txt 2d_$(SEED).ini
	@sh local_job.sh $(SEED)

%.txt: mk_%.jl
	@julia $< 

2d_$(SEED).ini: 2d.ini
	@cp 2d.ini 2d_$(SEED).ini

