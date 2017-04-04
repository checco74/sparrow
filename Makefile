CCOMP = g++
ANN_DIR = annlight
SOURCES = $(wildcard src/*.C)
DEPENDFILE = .depend
OPT = -O3 -funroll-loops -Wall
#OPT = -m64 -march=k8 -O3 -fprefetch-loop-arrays -mfpmath=sse -funroll-loops -Wall
LOPT_STATIC = -Wl,-Bstatic
LOPT_DYNAMIC = -Llib -lz -Wl,-Bdynamic
objx = CSVLogger.o scoreGear.o myClassifier.o LLt.o toolz.o aatype.o
vpath %.h include
vpath %.o bin
vpath %.C src
.SUFFIXES:

predMOD = LSMpred


all : depend sparrow computeConfInd ann

depend: $(SOURCES)
	$(CCOMP) -Iinclude -I/opt/acml3.6.0/gnu64/include -MM $(SOURCES) > $(DEPENDFILE)

-include $(DEPENDFILE)
	
sparrow : $(predMOD).o $(objx)
	$(CCOMP) -o sparrow bin/$(predMOD).o $(patsubst %.o,bin/%.o,$(objx)) $(OPT)

computeConfInd : computeConfInd.o $(objx)
	$(CCOMP) -o computeConfInd bin/computeConfInd.o $(patsubst %.o,bin/%.o,$(objx)) $(OPT)

bin/$(predMOD).o : src/$(predMOD).C
	$(CCOMP) -c src/$(predMOD).C -o bin/$(predMOD).o -Iinclude $(OPT)

bin/toolz.o : src/toolz.C
	$(CCOMP) -c src/toolz.C -Iinclude -o bin/toolz.o $(OPT)

bin/LLt.o : src/LLt.C
	$(CCOMP) -c src/LLt.C -Iinclude -o bin/LLt.o $(OPT)

bin/CSVLogger.o : src/CSVLogger.C
	$(CCOMP) -c src/CSVLogger.C -Iinclude -o bin/CSVLogger.o $(OPT)

bin/aatype.o : src/aatype.C
	$(CCOMP) -c src/aatype.C -Iinclude -o bin/aatype.o $(OPT)

bin/scoreGear.o : src/scoreGear.C
	$(CCOMP) -c src/scoreGear.C -Iinclude -o bin/scoreGear.o $(OPT)
        
bin/myClassifier.o : src/myClassifier.C
	$(CCOMP) -c src/myClassifier.C -Iinclude -o bin/myClassifier.o $(OPT)

bin/computeConfInd.o : src/computeConfInd.C
	$(CCOMP) -c src/computeConfInd.C -o bin/computeConfInd.o -Iinclude $(OPT)

ann :
	cd $(ANN_DIR) ; make ; cd ../

archive : clean reset
	tar -czvf sparrow.tar.gz include src bin $(ANN_DIR) Makefile *.sh

clean :
	rm -f *.log
	rm -f *.dat
	rm -f *.csv
	rm -f *.txt
	rm -f *.out
	rm -f *.niv
	rm -f *.tex
	rm -f dbg*
	cd $(ANN_DIR) ; make clean ; cd ../

reset :
	rm -f bin/*.o
	rm -f sparrow
	rm -f $(DEPENDFILE)
