CCOMP = g++
OPT = -g -O3 -Wall
#OPT = -O3 -Wall
#OPT = -O3 -Wall -pg -g
#OPT = -m64 -march=k8 -O3 -fprefetch-loop-arrays -mfpmath=sse -funroll-loops -fmove-all-movables -Wall
#OPT = -m64 -march=k8 -O3 -fprefetch-loop-arrays -mfpmath=sse -funroll-loops -Wall
#LINK_OPT = -pg
AR = ar

INCLUDE = include
SRC = src
BIN = bin
DOC = docs
vpath %.h $(INCLUDE)
vpath %.o $(BIN)
vpath %.cpp $(SRC)
.SUFFIXES:

VERSION = 1.0
ANN_OBJ = ann.o annDefs.o ParameterParser.o ParserException.o DefaultParameters.o ANN.o ANNException.o DiscriminantParameters.o KMeansClusterer.o
DEPENDFILE = .depend
SOURCES = $(wildcard $(SRC)/*.cpp)
ARCHIVE_NAME = ANN$(VERSION)-src.zip

all: depend ann

depend: $(SOURCES)
	$(CCOMP) -I$(INCLUDE) -MM $(SOURCES) > $(DEPENDFILE)

-include $(DEPENDFILE)

ann: $(ANN_OBJ)
	$(CCOMP) -o ann $(patsubst %.o,$(BIN)/%.o,$(ANN_OBJ)) $(LINK_OPT)

%.o: %.cpp
	$(CCOMP) -I$(INCLUDE) -o $(BIN)/$@ -c $< $(OPT)
doc:
	doxygen

clean :
	rm -f $(DEPENDFILE)
	rm -f $(BIN)/*.o
	for i in . $(SRC) $(INCLUDE); do rm -f $$i/*~; done
	for i in $(PROGRAMS); do rm -f $$i; done
	rm -f *.log
	rm -f *.out
	rm -f $(ARCHIVE_NAME)

docclean :
	rm -f -r $(DOC)/*

fullclean: clean docclean

archive: fullclean
	zip -r $(ARCHIVE_NAME) *
