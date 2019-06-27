#
# Unix/Linux Makefile Makefile for MATLAB interface to fundest
#

MEX=mex
MEXCFLAGS=-I.. -O #-g
# WHEN USING LAPACK, CHANGE THE NEXT TWO LINES TO WHERE YOUR COMPILED LAPACK/BLAS & F2C LIBS ARE!
LAPACKBLASLIBS_PATH=/usr/lib
F2CLIBS_PATH=/usr/local/lib
LEVMARLIBSPATH=/home/lourakis/levmar-src/levmar-2.6

INTFACESRCS=fundest.c
INTFACESOBJ=fundest.mexglx
LEVMARLIB=$(LEVMARLIBSPATH)/liblevmar.a
FUNDESTLIB=../libfundest.a
LAPACKLIBS=$(LAPACKBLASLIBS_PATH)/liblapack.a $(LAPACKBLASLIBS_PATH)/libblas.a $(F2CLIBS_PATH)/libf2c.a

LIBS=$(FUNDESTLIB) $(LEVMARLIB) $(LAPACKLIBS)

dummy: $(INTFACESOBJ)

$(INTFACESOBJ): $(INTFACESRCS) $(FUNDESTLIB)
	$(MEX) $(MEXCFLAGS) $(INTFACESRCS) $(LIBS)

clean:
	@rm -f $(INTFACESOBJ)

realclean cleanall: clean

depend:
	makedepend -f Makefile $(INTFACESRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

