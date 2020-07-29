OBJECTS= constructs.o 
CPPFLAGS=-g3 -pthread -Wall 


all: lamodel tests

run:
	make && ./lamodel && python vis.py

constructs.o: constructs.h

lamodel: $(OBJECTS) lamodel.o 
	g++  $(OBJECTS) lamodel.o -o ../lamodel -lm    -L/usr/local/lib/ # -lGLU -lglut 


tests: $(OBJECTS) tests.o 
	g++  $(OBJECTS) tests.o -o ../tests -lm    -L/usr/local/lib/ # -lGLU -lglut 


remake: clean all

cleanup:
	rm -f ../submit_lamodel.sh.e*
	rm -f ../submit_lamodel.sh.o*

clean: cleanup
	rm -f *.o lamodel.o*  ../lamodel ../tests

cdata:
	sh run_cluster.sh

data:
	sh run_serial.sh


