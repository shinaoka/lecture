SEED = $(word 3, $(shell grep seed 2d.ini))
OUTPUT := output$(SEED)

all: $(OUTPUT)

$(OUTPUT): temperatures.txt Jij.txt utriangles.txt 2d_$(SEED).ini init_non_eq_state.txt
	@sh local_job.sh $(SEED)

%.txt: mk_%.jl
	@julia $< 

2d_$(SEED).ini: 2d.ini
	@cp 2d.ini 2d_$(SEED).ini

clean: 
	@cp $(OUTPUT) bak_$(OUTPUT)
	@rm -f $(OUTPUT)


