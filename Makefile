#CFLAGS := -g -lm -Wall -pthread
CFLAGS += -msse -O3 -fstrict-aliasing -fomit-frame-pointer -Wall -pthread -ggdb
CXXFLAGS += -msse -O3 -fstrict-aliasing -fomit-frame-pointer -Wall -pthread -ggdb
LDLIBS += -lm -lz
LDFLAGS += -pthread
#CFLAGS := -lm -msse -O3 -fPIC -fopenmp -fstrict-aliasing -fomit-frame-pointer -Wall

all: LIBLINEAR LIBOCAS SOAP predictor/svm_libocas_classify qualrecal

qualrecal:
	make -C qualRecal

LIBOCAS: 
	make -C libocas_v093

LIBLINEAR:
	make -C liblinear-1.8mod

SOAP:
	make -C soap_1.11_patched

predictor/svm_libocas_classify: \
			predictor/svm_libocas_classify.o \
			predictor/cifinput.o \
			predictor/erf.o \
			predictor/firecrest_input.o \
			predictor/input.o \
			predictor/ipar_input.o \
			predictor/locs_input.o \
			gzstream/gzstream.o \
			libocas_v093/lib_svmlight_format.o \
			libocas_v093/libocas.o \
			libocas_v093/libqp_splx.o \
			predictor/BAMWriter.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)


libocas_v093/libocas.o: libocas_v093/libqp.h
libocas_v093/libqp_splx.o: libocas_v093/libocas.h 

predictor/svm_libocas_classify.o: predictor/input.h \
								  predictor/erf.h \
				  predictor/BAMWriter.h \
                                  libocas_v093/lib_svmlight_format.h \
                                  libocas_v093/libocas.h
predictor/cifinput.o : predictor/input.h
predictor/erf.o : predictor/erf.h
predictor/BAMWriter.o : predictor/BAMWriter.h
predictor/firecrest_input.o : predictor/input.h \
	                          gzstream/gzstream.h
predictor/input.o : predictor/input.h
predictor/ipar_input.o : predictor/input.h
predictor/locs_input.o : predictor/input.h

clean: 
		rm -f *.o libocas_v093/*.o predictor/*.o gzstream/*.o
		rm -f predictor/svm_libocas_classify
		make -C soap_1.11_patched clean
		make -C libocas_v093 clean
		make -C liblinear-1.8mod clean
