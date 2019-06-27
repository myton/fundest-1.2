* install levmar
git clone https://github.com/myton/levmar
cd levmar
mkdir build
cmake ..
sudo make install
* install fundest1.2
git clone https://github.com/myton/fundest-1.2.git
cd fundest-1.2
mkdir build
cd build
cmake ..
make 
./fundest_demo  ../test/matches.txt


    **************************************************************
                                FUNDEST
                              version 1.2
                          By Manolis Lourakis

                     Institute of Computer Science
            Foundation for Research and Technology - Hellas
                       Heraklion, Crete, Greece
    **************************************************************


GENERAL
This is fundest, a copylefted C library for non-linear, robust fundamental
matrix estimation from matched image point features. fundest requires my levmar
Levenberg-Marquardt non-linear least squares library and LAPACK/BLAS, available
from http://www.ics.forth.gr/~lourakis/levmar and http://www.netlib.org/clapack,
respectively. Precompiled LAPACK/BLAS libraries for Windows can be found here:
http://ylzhao.blogspot.gr/2013/10/blas-lapack-precompiled-binaries-for.html

Fundest also implements robust regression techniques for coping with outliers.
Note that fundest does *not* include any means for detecting and matching point
features between images. Such functionality can be supplied by other software
such as D. Lowe's SIFT or (in some cases) S. Birchfield's KLT.

Briefly, the approach implemented by fundest is the following:

1) Normalization of point coordinates to improve conditioning as
   described in
   R.I. Hartley "In Defense of the Eight-Point Algorithm",
   IEEE Trans. on PAMI, Vol. 19, No. 6, pp. 580-593, June 1997.

2) Least Median of Squares (LMedS) linear fitting (8-point algorithm)
   to detect outliers, see P.J. Rousseeuw, "Least Median of Squares Regression",
   Journal of the American Statistics Association, Vol. 79, No. 388, pp. 871-880,
   Dec. 1984. To ensure adequate spatial distribution of point quadruples
   over the image, LMedS random sampling employs the bucketing technique
   proposed in
   Z. Zhang, R. Deriche, O. Faugeras, Q.T. Luong, 
   "A Robust Technique for Matching Two Uncalibrated Images Through the
   Recovery of the Unknown Epipolar Geometry", INRIA RR-2273, May 1994 
   The rank-2 contraint is enforced a posteriori on the estimated F by
   setting its smallest singular value to zero.

3) Non-linear refinement of the fundamental matrix linear estimate by
   minimizing either:
   i)  the algebraic error
   ii) the symmetric distance of points from their epipolar lines in both images.
   iii) the Sampson error (see P.D. Sampson, "Fitting Conic Sections to ``very
        scattered'' Data: An Iterative Refinement of the Bookstein Algorithm",
        CGIP, Vol. 18, Is. 1, pp. 97-108, Jan. 1982.)
   Algorithm i) intrinsically enforces the rank-2 constraint. In both ii) and iii),
   the rank-2 constraint is enforced by a parameterization that expresses one row
   of the fundamental matrix as a linear combination of the other two. In all cases,
   the minimization is performed using the Levenberg-Marquardt algorithm, see
   M.I.A. Lourakis, "levmar: Levenberg-Marquardt Nonlinear Least Squares
   Algorithms in C/C++", http://www.ics.forth.gr/~lourakis/levmar/

   To identify degenerate data, fundest includes an implementation of GRIC
   which can assess the relative quality of a fundamental matrix and homography
   fit to the data. The homography can be fitted using my homest library,
   http://www.ics.forth.gr/~lourakis/homest/

For more details on fundamental matrices, projective geometry and such,
refer to R. Hartley and A. Zisserman, "Multiple View Geometry in
Computer Vision", Cambridge University Press, 2000-3


COMPILATION
 - Non-linear minimization in fundest requires the levmar library. Thus,
   before building it, you should make sure that levmar and LAPACK/BLAS
   are installed on your system. Detailed installation instructions can be
   found in levmar's distribution at http://www.ics.forth.gr/~lourakis/levmar

 - The preferred way to build fundest is through the CMake cross-platform
   build system. The included CMakeLists.txt file can be used to generate
   makefiles for Unix systems or project files for Windows systems.
   CMakeLists.txt defines some configuration variables that control
   certain aspects of fundest and can be modified from CMake's user
   interface. The values of these variables are automatically propagated
   to fundest.h after CMake runs. More information on how to use CMake
   can be found at http://www.cmake.org

 - fundest can also be built using the supplied makefiles. Platform-specific
   instructions are given below. Before compiling, you should edit the
   appropriate makefile to specify the location of your compiled levmar and
   LAPACK/BLAS libraries. You also might consider changing a few configuration
   options found at the top of fundest.h. See the accompanying comments for
   more details.

   -- On a Linux/Unix system, typing "make" will build both fundest and the
      demo program using gcc. Alternatively, if Intel's C++ compiler is
      installed, it can be used by typing "make -f Makefile.icc".

   -- Under Windows and if Visual C is installed & configured for command line
      use, type "nmake /f Makefile.vc" in a cmd window to build fundest and the
      demo program. In case of trouble, read the comments on top of Makefile.vc

MATLAB INTERFACE
fundest's distribution includes a matlab mex interface.
See the 'matlab' subdirectory for more information and examples of use.

USE
fundest() is the main routine for estimating a fundamental matrix. See the
comments in fundest.c for an explanation of its arguments; fundest_demo.c illustrates
an example of using fundest() with matched point pairs that are read from a text
file. An example of such point pairs is included in subdirectory test.

Typing 
./fundest_demo  test/matches.txt 

should produce output similar to

Fundamental matrix estimation using 314 image matches

Estimated F matrix [16 outliers, 5.10%]
-1.067272e-07 -2.342201e-06 0.001673826 
5.646679e-06 3.764683e-07 -0.02341995 
-0.003070482 0.02210014 0.9748294 

Elapsed time: 0.06 seconds, 60.00 msecs

Symmetric epipolar RMS and RMedS errors for input points: 0.30155 0.184734
        25%, 50% and 75% quartiles: 0.0835032 0.184903 0.340236


CONTACT
Send your comments/bug reports to _mysurname_ **at** ics dot forth dot gr
