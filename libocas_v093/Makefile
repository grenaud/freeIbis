CFLAGS += -msse -O3 -fPIC -fstrict-aliasing -fomit-frame-pointer -Wall -pthread
LDLIBS += -lm
LDFLAGS += -pthread

all: libocas.so svmocas msvmocas linclass

svmocas:	svmocas.c lib_svmlight_format.c sparse_mat.c ocas_helper.c ocas_helper.h libocas.h sparse_mat.h libocas.c
		$(CC) $(CFLAGS) -o $@ svmocas.c lib_svmlight_format.c sparse_mat.c ocas_helper.c libocas.c libqp_splx.c $(LDLIBS)

msvmocas:	msvmocas.c lib_svmlight_format.c sparse_mat.c ocas_helper.c ocas_helper.h libocas.h sparse_mat.h libocas.c
		$(CC) $(CFLAGS) -o $@ msvmocas.c lib_svmlight_format.c sparse_mat.c ocas_helper.c libocas.c libqp_splx.c $(LDLIBS)

linclass:	linclass.c lib_svmlight_format.c libocas.h 
		$(CC) $(CFLAGS) -o $@ linclass.c lib_svmlight_format.c   $(LDLIBS)


libocas.so:	libocas.c libocas.h libqp_splx.c libqp.h 
		$(CC) $(CFLAGS) -shared -o $@ libocas.c libqp_splx.c

clean: 
		rm -f *~ svmocas.$(MEXSUFFIX) svmlight_linclass.$(MEXSUFFIX) svmocas msvmocas linclass libocas.so msvmocas.$(MEXSUFFIX) lbpfilter.$(MEXSUFFIX) lbppyr_features.$(MEXSUFFIX) svmocas_lbp.$(MEXSUFFIX) lbppyr.$(MEXSUFFIX) compute_errors.$(MEXSUFFIX) svm_libocas_classify
