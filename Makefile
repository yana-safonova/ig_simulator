all:
	mkdir -p bin
	g++ src/ideal_repertoire_constructor/main_ideal_repertoire.cpp -std=c++11 -o bin/ideal_repertoire_constructor
	g++ src/ig_simulator/main.cpp -std=c++11 -o bin/ig_simulator
	g++ src/paired_read_merger/main.cpp -std=c++11 -o bin/paired_read_merger
	g++ src/art_src_GreatSmokyMountains_Linux/art_illumina_src/*.cpp -std=c++11 -o bin/art_illumina

