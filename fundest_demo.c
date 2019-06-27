/////////////////////////////////////////////////////////////////////////////////
// 
//  Fundamental matrix estimation from 2D point matches
//  Copyright (C) 2002 - 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

/******************************************************************************** 
 * fundest demo. The program accepts a text file containing quadruples of matching
 * point coordinates (i.e., x0 y0 x1 y1 and estimates the underlying fundamental
 * matrix
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "fundest.h"

#define MAXSTRLEN 1024

/* read matching points from a file */
static int readMatchingPoints(char *fname, double (**pts0)[2], double (**pts1)[2])
{
register int i;
int ncoords, nmatches;
double coords[4];
FILE *fp;
char buf[MAXSTRLEN];

  if((fp=fopen(fname, "r"))==NULL){
    fprintf(stderr, "cannot open file %s\n", fname);
    exit(1);
  }

  fgets(buf, MAXSTRLEN, fp);
  if(ferror(fp)){
    fprintf(stderr, "File %s: error reading first line\n", fname);
    exit(1);
  }

  ncoords=sscanf(buf, "%lf%lf%lf%lf", coords, coords+1, coords+2, coords+3);
  if(ncoords==4){ /* no lines number */
    for(nmatches=1; !feof(fp); nmatches++){
      fscanf(fp, "%*g%*g%*g%*g\n");
      if(ferror(fp)){
        fprintf(stderr, "File %s: error reading point coordinates, line %d\n", fname, nmatches + 1);
        exit(1);
      }
    }

    rewind(fp);
  }
  else{
    sscanf(buf, "%d", &nmatches);
  }

  *pts0=(double (*)[2])malloc(nmatches*sizeof(double[2]));
  *pts1=(double (*)[2])malloc(nmatches*sizeof(double[2]));
  if(!pts0 || !pts1){
    fprintf(stderr, "Memory allocation request failed in readMatchingPoints()\n");
    exit(1);
  }

  /* read in points and store them */
  for(i=0; !feof(fp); i++){
    ncoords=fscanf(fp, "%lf%lf%lf%lf\n", (*pts0)[i], (*pts0)[i]+1, (*pts1)[i], (*pts1)[i]+1);
    if(ncoords==EOF) break;

    if(ncoords!=4){
      fprintf(stderr, "File %s: line %d contains only %d coordinates\n", fname, i + 1, ncoords);
      exit(1);
    }

    if(ferror(fp)){
      fprintf(stderr, "File %s: error reading point coordinates, line %d\n", fname, i + 1);
      exit(1);
    }

  }
  fclose(fp);

  if(i!=nmatches){
    fprintf(stderr, "number of actuall points in file %s does not agree with that in first line (%d != %d)!\n",
                     fname, i, nmatches);
    exit(1);
  }

  return nmatches;
}


#define INL_PCENT 0.65

#ifdef HAVE_HOMEST
#include <homest.h>
#endif

int main(int argc, char *argv[])
{
double (*pts0)[2], (*pts1)[2];
register int i;
int npts, noutl, *outidx=NULL;
char *matchesfile;
double F[NUM_FPARAMS], rms, rmeds;
int donorm=1;

clock_t start_time, end_time;

  /* arguments parsing */
  if(argc!=2){
    fprintf(stderr, "Usage: %s <matched points>\n", argv[0]);
    exit(1);
  }

  matchesfile=argv[1];

  npts=readMatchingPoints(matchesfile, &pts0, &pts1);

#if 0
  for(i=0; i<npts; ++i){
    printf("%g %g  %g %g\n", pts1[i][0], pts1[i][1], pts0[i][0], pts0[i][1]);
  }
#endif

#ifdef NEED_OUTLIERS
  if((outidx=(int *)malloc(npts*sizeof(int)))==NULL){
    fprintf(stderr, "Memory allocation request failed in main()\n");
    exit(1);
  }
#endif /* NEED_OUTLIERS */

  fprintf(stdout, "Fundamental matrix estimation using %d image matches\n", npts);

  start_time=clock();
  {
    int cstfunc;

    //cstfunc=FUNDEST_8PT; // just linear estimation
    //cstfunc=FUNDEST_ALGMIN; // algebraic error minimization
    //cstfunc=FUNDEST_EPIP_DIST; // use the symmetric distances from the epipolar lines
    cstfunc=FUNDEST_SAMPSON_ERROR; // use the Sampson approximation to geometric error
    i=fundest(pts0, pts1, npts, INL_PCENT, F, donorm, cstfunc, outidx, &noutl, 1);
  }
  end_time=clock();

  if(i!=FUNDEST_OK){
    fprintf(stderr, "Fundamental matrix estimation failed!\n"
                    "Make sure that the matches do not originate from a planar surface or a pure rotation.\n");
    goto cleanup;
  }

  // sample fundest_wie() call:
  // fundest_wie(pts0, pts1, npts, -0.8, F0, FUNDEST_8PT, outidx, &noutl, 1); // keeps 80% inliers, F0 should be initialized


  fprintf(stdout, "\nEstimated F matrix [%d outliers, %.2lf%%]\n", noutl, (double)(100.0*noutl)/npts);
  for(i=0; i<NUM_FPARAMS; ++i){
    if(i && !(i%3)) putchar('\n');
    fprintf(stdout, "%.7g ", F[i]);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "\nElapsed time: %.2lf seconds, %.2lf msecs\n", ((double) (end_time - start_time)) / CLOCKS_PER_SEC,
                        ((double) (end_time - start_time)) / (CLOCKS_PER_SEC/1000.0));


  fundest_RMS_RMedS(pts0, pts1, npts, F, &rms, &rmeds);
  fprintf(stdout, "\nSymmetric epipolar RMS and RMedS errors for input points: %g %g\n", rms, rmeds);
  fflush(stdout);

  { double q1, q2, q3;

    fundest_quartiles(pts0, pts1, npts, F, &q1, &q2, &q3);
    fprintf(stdout, "\t25%%, 50%% and 75%% quartiles: %g %g %g\n", q1, q2, q3);
    fflush(stdout);
  }

#ifdef NEED_OUTLIERS
  fprintf(stdout, "Indices of the %d outlying pairs:\n", noutl);
  for(i=0; i<noutl; ++i)
    fprintf(stdout, "%d ", outidx[i]);
  fputc('\n', stdout);
#endif /* NEED_OUTLIERS */

#ifdef HAVE_HOMEST
  {
  double H[9];
  int noutl, donorm=1;
  double gF, gH, s=0.8;

  homest(pts0, pts1, npts, INL_PCENT, H, donorm, HOMEST_SYM_XFER_ERROR, NULL, &noutl, 0);
  fundest_GRIC(pts0, pts1, npts, F, H, s, &gF, &gH);
  fprintf(stdout, "\nGRICs for fundamental matrix and homography are resp. %g, %g\n", gF, gH);
  }
#else
  {
  double gF, s=0.8;

  fundest_GRIC(pts0, pts1, npts, F, NULL, s, &gF, NULL);
  fprintf(stdout, "\nGRIC for fundamental matrix %g\n", gF);
  }
#endif

cleanup:
  if(outidx) free(outidx);

  free(pts0);
  free(pts1);

  return 0;
}
