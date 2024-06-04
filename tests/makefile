all: build

build:
	./gen_ideal_test.py

cp:
	cp ../src/icar .

run: cp
	./icar icar_options.nml

post:
	../helpers/aggregate_parallel_files.py

clean:
	rm -f forcing.nc init.nc icar_options.nml *~
cleanall: clean
	rm -f icar *.nc
