
SEED = $(word 3, $(shell grep seed 2d.ini))

all:
	@echo $(SEED)

