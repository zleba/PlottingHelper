CXX = g++ -g -Wall -std=c++1y -O3
CFLAGS = $(shell root-config --cflags)
LIBS := $(shell root-config --libs)

all: libPlottingHelper.so plottingHelper_C.so

plottingHelper_C.so: plottingHelper.C
	root -l -b -q $<++g

libPlottingHelper.so: plottingHelper.o
	${CXX}  -shared -Wl,-soname,$@ -o $@   $< ${LIBS} 

plottingHelper.o: plottingHelper.C
	${CXX}  ${CFLAGS} -fPIC -c $< ${LIBS} 

clean: 
	rm -f libPlottingHelper.so plottingHelper.o plottingHelper_C.d plottingHelper_C_ACLiC_dict_rdict.pcm plottingHelper_C.so



#libplottingHelperNew.so: plottingHelper.C
#	clang++ -O0  ${CFLAGS} -fPIC -c $<
#	clang++ -O0 -shared -Wl,-soname,libplottingHelperNew.so -o libplottingHelperNew.so   plottingHelper.o

#test: test.cpp 
#	${CXX} -Wall ${CFLAGS} $<  -Wl,-rpath,. plottingHelper_C.so -o test  ${LIBS} 
#
#testFast: test.cpp libplottingHelper.so
#	${CXX} -Wall ${CFLAGS} $<  -Wl,-rpath,. -L. -lplottingHelper -o testFast  ${LIBS} 
