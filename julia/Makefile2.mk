
INPUTS  = $(shell ls output*)
OUTPUTS = $(INPUTS:output%=output%.dat)

all: $(OUTPUTS)

output%.dat:
	@sh local_job.sh

2d_hoge.ini: temperatures_hoge.txt

temperatures_hoge.txt: init_temperatures_hoge.txt
	@julia insert_temperatures.jl

init_temperatures_hoge.txt:
	@head -n 1 init_temperatures.txt  > init_temperatures_hoge.txt
	@grep sh hoge | awk '{print $2}' >> init_temperatures_hoge.txt


