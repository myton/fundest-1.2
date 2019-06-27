#
# Unix/Linux GCC Makefile for robust fmatrix estimation
#
# The homest library can be optionally linked with homest_demo to enable
# the estimation of a homography and then of the corresponding GRIC
# 

CC=gcc
LEVMAR_PATH=/home/lourakis/levmar-src/levmar-2.6/ # CHANGE THIS TO POINT TO YOUR COMPILED COPY OF LEVMAR
#HOMEST_PATH=/home/lourakis/homest-src/homest-1.5/ # CHANGE THIS TO POINT TO YOUR COMPILED COPY OF HOMEST
#HOMEST_CC=-DHAVE_HOMEST -I$(HOMEST_PATH)
#HOMEST_LD=-L$(HOMEST_PATH) -lhomest
INCLUDES=-I$(LEVMAR_PATH)
CFLAGS=$(HOMEST_CC) $(INCLUDES) -O3 -funroll-loops -Wall #-g #-pg
LAPACKLIBS_PATH=/usr/local/lib # WHEN USING LAPACK, CHANGE THIS TO WHERE YOUR COMPILED LIBS ARE!
LDFLAGS=-L. -L$(LEVMAR_PATH) -L$(LAPACKLIBS_PATH)

LIBOBJS=calc_fund_coeffs.o fundest.o lqs.o buckets.o norm.o linalg.o gric.o misc.o
LIBSRCS=calc_fund_coeffs.c fundest.c lqs.c buckets.c norm.c linalg.c gric.c misc.c

DEMOOBJS=fundest_demo.o
DEMOSRCS=fundest_demo.c
AR=ar
RANLIB=ranlib
RM=rm -f
MAKE=make
MAPLE=maple

LAPACKLIBS=-llapack -lblas -lf2c

all: libfundest.a fundest_demo

libfundest.a: $(LIBOBJS)
	$(AR) crv libfundest.a $(LIBOBJS)
	$(RANLIB) libfundest.a

calc_fund_coeffs.o: calc_fund_coeffs.h
fundest.o: util.h fundest.h calc_fund_coeffs.h lqs.h
lqs.o: lqs.h compiler.h
linalg.o: compiler.h util.h
norm.o: compiler.h util.h
buckets.o: util.h lqs.h

fundest_demo: $(DEMOOBJS) libfundest.a
	$(CC) $(LDFLAGS) $(DEMOOBJS) -o fundest_demo -lfundest $(HOMEST_LD) -llevmar $(LAPACKLIBS) -lm 

fundest_demo.o: fundest.h

#calc_fund_coeffs.c: clc-fundcoeffs.mpl
#	$(MAPLE) <  $<

clean:
	-$(RM) $(LIBOBJS) $(DEMOOBJS) gmon.out

cleanall: clean
	-$(RM) libfundest.a

depend:
	makedepend $(INCLUDES) -f Makefile $(LIBSRCS) $(DEMOSRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

