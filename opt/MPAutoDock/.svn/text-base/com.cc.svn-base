/*

 $Id: com.cc,v 1.3 2010/07/30 23:51:12 pkcoff Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "ranlib.h"
#include "structs.h"
#include "constants.h"

extern FILE *logFile;
extern int local_search_begun;

//pkcoff - openmp for bluegene change section - start

// There are many global variables used in many global functions here, due to the nature of how they have been declared in this file they couldn't be made threadprivate, therefore there were made into arrays of MAX_NUM_THREADS length and indexed by threadID
extern int threadID;

extern void slaveExit(int rc);

#pragma omp threadprivate (threadID)

//pkcoff - openmp for bluegene change section - end


void advnst(FourByteLong k)
/*
**********************************************************************
     void advnst(FourByteLong k)
               ADV-a-N-ce ST-ate
     Advances the state  of  the current  generator  by 2^K values  and
     resets the initial seed to that value.
     This is  a  transcription from   Pascal to  Fortran    of  routine
     Advance_State from the paper
     L'Ecuyer, P. and  Cote, S. "Implementing  a  Random Number Package
     with  Splitting   Facilities."  ACM  Transactions  on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     k -> The generator is advanced by2^K values
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(FourByteLong getset,FourByteLong *qvalue);
extern void gscgn(FourByteLong getset,FourByteLong *g);
extern FourByteLong Xm1[],Xm2[],Xa1[],Xa2[],Xcg1[MAX_NUM_THREADS][32],Xcg2[MAX_NUM_THREADS][32];
static FourByteLong g,i,ib1,ib2;
static FourByteLong qrgnin;


//pkcoff - openmp for bluegene change section - start
// these static variables in this global function must be thread private

#pragma omp threadprivate (g,i,ib1,ib2,qrgnin)

//pkcoff - openmp for bluegene change section - end


/*
     Abort unless random number generator initialized
*/
    gsrgs(0L,&qrgnin);
    if(qrgnin) goto S10;
    fputs(" ADVNST called before random generator initialized - ABORT",stderr);
    slaveExit(1);
S10:
    gscgn(0L,&g);
    ib1 = Xa1[threadID];
    ib2 = Xa2[threadID];
    for (i=1; i<=k; i++) {
        ib1 = mltmod(ib1,ib1,Xm1[threadID]);
        ib2 = mltmod(ib2,ib2,Xm2[threadID]);
    }
    setsd(mltmod(ib1,*(Xcg1[threadID]+g-1),Xm1[threadID]),mltmod(ib2,*(Xcg2[threadID]+g-1),Xm2[threadID]));
/*
     NOW, IB1 = A1**K AND IB2 = A2**K
*/
#undef numg
}
void getsd(FourByteLong *iseed1,FourByteLong *iseed2)
/*
**********************************************************************
     void getsd(FourByteLong *iseed1,FourByteLong *iseed2)
               GET SeeD
     Returns the value of two integer seeds of the current generator
     This  is   a  transcription from  Pascal   to  Fortran  of routine
     Get_State from the paper
     L'Ecuyer, P. and  Cote,  S. "Implementing a Random Number  Package
     with   Splitting Facilities."  ACM  Transactions   on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 <- First integer seed of generator G
     iseed2 <- Second integer seed of generator G
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(FourByteLong getset,FourByteLong *qvalue);
extern void gscgn(FourByteLong getset,FourByteLong *g);
extern FourByteLong Xcg1[MAX_NUM_THREADS][32],Xcg2[MAX_NUM_THREADS][32];
static FourByteLong g;
static FourByteLong qrgnin;


//pkcoff - openmp for bluegene change section - start
// these static variables in this global function must be thread private

#pragma omp threadprivate (g,qrgnin)

//pkcoff - openmp for bluegene change section - end

/*
     Abort unless random number generator initialized
*/
    gsrgs(0L,&qrgnin);
    if(qrgnin) goto S10;
    fprintf(stderr,"%s\n",
      " GETSD called before random number generator  initialized -- abort!");
    slaveExit(0);
S10:
    gscgn(0L,&g);
    *iseed1 = *(Xcg1[threadID]+g-1);
    *iseed2 = *(Xcg2[threadID]+g-1);
#undef numg
}
FourByteLong ignlgi(void)
/*
**********************************************************************
     FourByteLong ignlgi(void)
               GeNerate LarGe Integer
     Returns a random integer following a uniform distribution over
     (1, 2147483562) using the current generator.
     This is a transcription from Pascal to Fortran of routine
     Random from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(FourByteLong getset,FourByteLong *qvalue);
extern void gssst(FourByteLong getset,FourByteLong *qset);
extern void gscgn(FourByteLong getset,FourByteLong *g);
extern void inrgcm(void);
extern FourByteLong Xm1[],Xm2[],Xa1[],Xa2[],Xcg1[MAX_NUM_THREADS][32],Xcg2[MAX_NUM_THREADS][32];
extern FourByteLong Xqanti[MAX_NUM_THREADS][32];
static FourByteLong ignlgi,curntg,k,s1,s2,z;
static FourByteLong qqssd,qrgnin;


//pkcoff - openmp for bluegene change section - start
// these static variables in this global function must be thread private

#pragma omp threadprivate (ignlgi,curntg,k,s1,s2,z,qqssd,qrgnin )

//pkcoff - openmp for bluegene change section - end

/*
     IF THE RANDOM NUMBER PACKAGE HAS NOT BEEN INITIALIZED YET, DO SO.
     IT CAN BE INITIALIZED IN ONE OF TWO WAYS : 1) THE FIRST CALL TO
     THIS ROUTINE  2) A CALL TO SETALL.
*/
    gsrgs(0L,&qrgnin);
    if(!qrgnin) inrgcm();
    gssst(0,&qqssd);
    if(!qqssd) setall(1234567890L,123456789L);
/*
     Get Current Generator
*/
    gscgn(0L,&curntg);
    s1 = *(Xcg1[threadID]+curntg-1);
    s2 = *(Xcg2[threadID]+curntg-1);
    k = s1/53668L;
    s1 = Xa1[threadID]*(s1-k*53668L)-k*12211;
    if(s1 < 0) s1 += Xm1[threadID];
    k = s2/52774L;
    s2 = Xa2[threadID]*(s2-k*52774L)-k*3791;
    if(s2 < 0) s2 += Xm2[threadID];
    *(Xcg1[threadID]+curntg-1) = s1;
    *(Xcg2[threadID]+curntg-1) = s2;
    z = s1-s2;
    if(z < 1) z += (Xm1[threadID]-1);
    if(*(Xqanti[threadID]+curntg-1)) z = Xm1[threadID]-z;
    ignlgi = z;
    
//    if (((threadID == 1) || (threadID == 0) || (threadID == 2)) && local_search_begun)
//      printf("threadID: %d ignlgi: %d\n",threadID,ignlgi);
    
    
    return ignlgi;
#undef numg
}
void initgn(FourByteLong isdtyp)
/*
**********************************************************************
     void initgn(FourByteLong isdtyp)
          INIT-ialize current G-e-N-erator
     Reinitializes the state of the current generator
     This is a transcription from Pascal to Fortran of routine
     Init_Generator from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     isdtyp -> The state to which the generator is to be set
          isdtyp = -1  => sets the seeds to their initial value
          isdtyp =  0  => sets the seeds to the first value of
                          the current block
          isdtyp =  1  => sets the seeds to the first value of
                          the next block
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(FourByteLong getset,FourByteLong *qvalue);
extern void gscgn(FourByteLong getset,FourByteLong *g);
extern FourByteLong Xm1[],Xm2[],Xa1w[],Xa2w[],Xig1[MAX_NUM_THREADS][32],Xig2[MAX_NUM_THREADS][32],Xlg1[MAX_NUM_THREADS][32],Xlg2[MAX_NUM_THREADS][32],Xcg1[MAX_NUM_THREADS][32],Xcg2[MAX_NUM_THREADS][32];
static FourByteLong g;
static FourByteLong qrgnin;


//pkcoff - openmp for bluegene change section - start
// these static variables in this global function must be thread private

#pragma omp threadprivate (g,qrgnin )

//pkcoff - openmp for bluegene change section - end

/*
     Abort unless random number generator initialized
*/
    gsrgs(0L,&qrgnin);
    if(qrgnin) goto S10;
    fprintf(stderr,"%s\n",
      " INITGN called before random number generator  initialized -- abort!");
    slaveExit(1);
S10:
    gscgn(0L,&g);
    if(-1 != isdtyp) goto S20;
    *(Xlg1[threadID]+g-1) = *(Xig1[threadID]+g-1);
    *(Xlg2[threadID]+g-1) = *(Xig2[threadID]+g-1);
    goto S50;
S20:
    if(0 != isdtyp) goto S30;
    goto S50;
S30:
/*
     do nothing
*/
    if(1 != isdtyp) goto S40;
    *(Xlg1[threadID]+g-1) = mltmod(Xa1w[threadID],*(Xlg1[threadID]+g-1),Xm1[threadID]);
    *(Xlg2[threadID]+g-1) = mltmod(Xa2w[threadID],*(Xlg2[threadID]+g-1),Xm2[threadID]);
    goto S50;
S40:
    fprintf(stderr,"%s\n","isdtyp not in range in INITGN");
    slaveExit(1);
S50:
    *(Xcg1[threadID]+g-1) = *(Xlg1[threadID]+g-1);
    *(Xcg2[threadID]+g-1) = *(Xlg2[threadID]+g-1);
#undef numg
}
void inrgcm(void)
/*
**********************************************************************
     void inrgcm(void)
          INitialize Random number Generator CoMmon
                              Function
     Initializes common area  for random number  generator.  This saves
     the  nuisance  of  a  BLOCK DATA  routine  and the  difficulty  of
     assuring that the routine is loaded with the other routines.
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(FourByteLong getset,FourByteLong *qvalue);
extern FourByteLong Xm1[],Xm2[],Xa1[],Xa2[],Xa1w[],Xa2w[],Xa1vw[],Xa2vw[];
extern FourByteLong Xqanti[MAX_NUM_THREADS][32];
static FourByteLong T1;
static FourByteLong i;

//pkcoff - openmp for bluegene change section - start
// these static variables in this global function must be thread private

#pragma omp threadprivate (T1,i )

//pkcoff - openmp for bluegene change section - end

/*
     V=20;                            W=30;
     A1W = MOD(A1**(2**W),M1)         A2W = MOD(A2**(2**W),M2)
     A1VW = MOD(A1**(2**(V+W)),M1)    A2VW = MOD(A2**(2**(V+W)),M2)
   If V or W is changed A1W, A2W, A1VW, and A2VW need to be recomputed.
    An efficient way to precompute a**(2*j) MOD m is to start with
    a and square it j times modulo m using the function MLTMOD.
*/



for (int iter1=0;iter1<MAX_NUM_THREADS;iter1++) {
    Xm1[iter1] = 2147483563L;
    Xm2[iter1] = 2147483399L;
    Xa1[iter1] = 40014L;
    Xa2[iter1] = 40692L;
    Xa1w[iter1] = 1033780774L;
    Xa2w[iter1] = 1494757890L;
    Xa1vw[iter1] = 2082007225L;
    Xa2vw[iter1] = 784306273L;
    for (i=0; i<numg; i++) *(Xqanti[iter1]+i) = 0;


}

/*
    Xm1[threadID] = 2147483563L;
    Xm2[threadID] = 2147483399L;
    Xa1[threadID] = 40014L;
    Xa2[threadID] = 40692L;
    Xa1w[threadID] = 1033780774L;
    Xa2w[threadID] = 1494757890L;
    Xa1vw[threadID] = 2082007225L;
    Xa2vw[threadID] = 784306273L;
    for (i=0; i<numg; i++) *(Xqanti[threadID]+i) = 0;
*/

    T1 = 1;
/*
     Tell the world that common has been initialized
*/
    gsrgs(1L,&T1);
#undef numg
}
void setall(FourByteLong iseed1,FourByteLong iseed2)
/*
**********************************************************************
     void setall(FourByteLong iseed1,FourByteLong iseed2)
               SET ALL random number generators
     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
     initial seeds of the other generators are set accordingly, and
     all generators states are set to these seeds.
     This is a transcription from Pascal to Fortran of routine
     Set_Initial_Seed from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 -> First of two integer seeds
     iseed2 -> Second of two integer seeds
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(FourByteLong getset,FourByteLong *qvalue);
extern void gssst(FourByteLong getset,FourByteLong *qset);
extern void gscgn(FourByteLong getset,FourByteLong *g);
extern FourByteLong Xm1[],Xm2[],Xa1vw[],Xa2vw[],Xig1[MAX_NUM_THREADS][32],Xig2[MAX_NUM_THREADS][32];
static FourByteLong T1;
static FourByteLong g,ocgn;
static FourByteLong qrgnin;


//pkcoff - openmp for bluegene change section - start
// these static variables in this global function must be thread private

#pragma omp threadprivate (T1,g,ocgn,qrgnin )

//pkcoff - openmp for bluegene change section - end

    T1 = 1;
/*
     TELL IGNLGI, THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
      HAS BEEN CALLED.
*/
    gssst(1,&T1);
    gscgn(0L,&ocgn);
/*
     Initialize Common Block if Necessary
*/
    gsrgs(0L,&qrgnin);
    if(!qrgnin) inrgcm();
    *Xig1[threadID] = iseed1;
    *Xig2[threadID] = iseed2;
    initgn(-1L);
    for (g=2; g<=numg; g++) {
        *(Xig1[threadID]+g-1) = mltmod(Xa1vw[threadID],*(Xig1[threadID]+g-2),Xm1[threadID]);
        *(Xig2[threadID]+g-1) = mltmod(Xa2vw[threadID],*(Xig2[threadID]+g-2),Xm2[threadID]);
        gscgn(1L,&g);
        initgn(-1L);
    }
    gscgn(1L,&ocgn);
#undef numg
}
void setant(FourByteLong qvalue)
/*
**********************************************************************
     void setant(FourByteLong qvalue)
               SET ANTithetic
     Sets whether the current generator produces antithetic values.  If
     X   is  the value  normally returned  from  a uniform [0,1] random
     number generator then 1  - X is the antithetic  value. If X is the
     value  normally  returned  from a   uniform  [0,N]  random  number
     generator then N - 1 - X is the antithetic value.
     All generators are initialized to NOT generate antithetic values.
     This is a transcription from Pascal to Fortran of routine
     Set_Antithetic from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     qvalue -> nonzero if generator G is to generating antithetic
                    values, otherwise zero
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(FourByteLong getset,FourByteLong *qvalue);
extern void gscgn(FourByteLong getset,FourByteLong *g);
extern FourByteLong Xqanti[MAX_NUM_THREADS][32];
static FourByteLong g;
static FourByteLong qrgnin;


//pkcoff - openmp for bluegene change section - start
// these static variables in this global function must be thread private

#pragma omp threadprivate (g,qrgnin )

//pkcoff - openmp for bluegene change section - end

/*
     Abort unless random number generator initialized
*/
    gsrgs(0L,&qrgnin);
    if(qrgnin) goto S10;
    fprintf(stderr,"%s\n",
      " SETANT called before random number generator  initialized -- abort!");
    slaveExit(1);
S10:
    gscgn(0L,&g);
    Xqanti[threadID][g-1] = qvalue;
#undef numg
}
void setsd(FourByteLong iseed1,FourByteLong iseed2)
/*
**********************************************************************
     void setsd(FourByteLong iseed1,FourByteLong iseed2)
               SET S-ee-D of current generator
     Resets the initial  seed of  the current  generator to  ISEED1 and
     ISEED2. The seeds of the other generators remain unchanged.
     This is a transcription from Pascal to Fortran of routine
     Set_Seed from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 -> First integer seed
     iseed2 -> Second integer seed
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(FourByteLong getset,FourByteLong *qvalue);
extern void gscgn(FourByteLong getset,FourByteLong *g);
extern FourByteLong Xig1[MAX_NUM_THREADS][32],Xig2[MAX_NUM_THREADS][32];
static FourByteLong g;
static FourByteLong qrgnin;

//pkcoff - openmp for bluegene change section - start
// these static variables in this global function must be thread private

#pragma omp threadprivate (g,qrgnin )

//pkcoff - openmp for bluegene change section - end

/*
     Abort unless random number generator initialized
*/
    gsrgs(0L,&qrgnin);
    if(qrgnin) goto S10;
    fprintf(stderr,"%s\n",
      " SETSD called before random number generator  initialized -- abort!");
    slaveExit(1);
S10:
    gscgn(0L,&g);
    *(Xig1[threadID]+g-1) = iseed1;
    *(Xig2[threadID]+g-1) = iseed2;
    initgn(-1L);
#undef numg
}

//pkcoff - openmp for bluegene change section - start

// There are many global variables used in many global functions here, due to the nature of how they have been declared in this file they couldn't be made threadprivate, therefore there were made into arrays of MAX_NUM_THREADS length and indexed by threadID

FourByteLong Xm1[MAX_NUM_THREADS],Xm2[MAX_NUM_THREADS],Xa1[MAX_NUM_THREADS],Xa2[MAX_NUM_THREADS],Xcg1[MAX_NUM_THREADS][32],Xcg2[MAX_NUM_THREADS][32],Xa1w[MAX_NUM_THREADS],Xa2w[MAX_NUM_THREADS],Xig1[MAX_NUM_THREADS][32],Xig2[MAX_NUM_THREADS][32],Xlg1[MAX_NUM_THREADS][32],
    Xlg2[MAX_NUM_THREADS][32],Xa1vw[MAX_NUM_THREADS],Xa2vw[MAX_NUM_THREADS];
FourByteLong Xqanti[MAX_NUM_THREADS][32];

//pkcoff - openmp for bluegene change section - end
