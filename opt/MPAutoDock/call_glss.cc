
/*

 $Id: call_glss.cc,v 1.3 2010/07/30 23:51:12 pkcoff Exp $

 AutoDock  

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
 All Rights Reserved.

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

/********************************************************************
     Call_glss:  Invokes a GA-LS hybrid to try and solve the
                 docking problem.

                                rsh 9/95
********************************************************************/

#include <string.h>
#include "gs.h"
#include "ls.h"
#include "support.h"
#include "eval.h"
#include "hybrids.h"
#include "constants.h"
#include "structs.h"
#include "openfile.h"
#include "qmultiply.h"

extern FILE *logFile;
extern char *programname;

//pkcoff - openmp for bluegene change section - start

// this file is where most of the time goes during a run, hence the CALL_GLSS_PROFILING directives
//putting this with the rest of the globals in autoglobal.h
extern int global_ntor;

// evaluate is now an array and declared in autoglobal.h
// all references to global evaluate oject in this file have been replaced with evaluate[threadID] to use the correct array element corresponding to the thread ID
extern Eval evaluate[];
extern int threadID;

extern char searchResultsBuffer[MAX_NUM_THREADS][DLG_BUFFER_LENGTH], searchResultsLine[MAX_NUM_THREADS][DLG_LINE_LENGTH];

// these need to be thread private
#pragma omp threadprivate (global_ntor,threadID)

//pkcoff - openmp for bluegene change section - end


#ifdef CALL_GLSS_PROFILING
#include "timesyshms.h"
struct tms tms_callGlssStart;
struct tms tms_callGlssSubStart;
struct tms tms_callGlssSubSubStart;
struct tms tms_callGlssSubSubSubStart;
struct tms tms_callGlssEnd;
#endif

Representation **generate_R(int num_torsions, GridMapSetInfo *info)
{
   Representation **retval;
   Quat q;

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  about to create a new Representation with 5 elements, retval...\n");
#endif
   retval = new Representation *[5];
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  done creating   a new Representation with 5 elements, retval...\n");
#endif
   // Set the x-translation
   retval[0] = new RealVector( 1, info->lo[X], info->hi[X] );
   // Set the y-translation
   retval[1] = new RealVector( 1, info->lo[Y], info->hi[Y] );
   // Set the z-translation
   retval[2] = new RealVector( 1, info->lo[Z], info->hi[Z] );

   // Generate a uniformly-distributed random quaternion for a random rotation (UDQ)
   q = uniformQuat();
   q = convertQuatToRot( q );
#ifdef DEBUG
   printQuat( logFile, q );
#endif

   // Set the unit vector components (the "axis"), for the rotation about axis
   retval[3] = new RealVector( 3, -1., 1., q.nx, q.ny, q.nz ); // uniformly-distributed quaternion (UDQ)

   // Set the angle (the "rotation") for the rotation about axis, 
   // and any torsion angles
   retval[4] = new RealVector( num_torsions+1, -PI, PI, q.ang ); // uniformly-distributed quaternion (UDQ)
   // retval[4] = new RealVector( num_torsions+1, -PI, PI );  // rotation-about-axis angle is uniformly distributed, -PI to PI, not UDQ

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  done assigning each of the retval[0-5] elements...\n");
#endif

   return(retval);
}

Representation **generate_R_quaternion(int num_torsions, GridMapSetInfo *info)
{
   Representation **retval;
   Quat q;

#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/Representation **generate_R_quaternion()  about to create a new Representation with 5 elements, retval...\n");
#endif
   retval = new Representation *[5];
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R_quaternion()  done creating   a new Representation with 5 elements, retval...\n");
#endif
   // Set the x-translation
   retval[0] = new RealVector( 1, info->lo[X], info->hi[X] );
   // Set the y-translation
   retval[1] = new RealVector( 1, info->lo[Y], info->hi[Y] );
   // Set the z-translation
   retval[2] = new RealVector( 1, info->lo[Z], info->hi[Z] );

   // Generate a uniformly-distributed random quaternion for a random rotation (UDQ)
   q = uniformQuat();
   q = convertQuatToRot( q );
#ifdef DEBUG
   printQuat( logFile, q );
#endif

#ifdef DEBUG_QUAT
#ifdef DEBUG_QUAT_PRINT
    pr( logFile, "DEBUG_QUAT: generate_R_quaternion()\n" );
    (void) fflush(logFile);
#endif
    //  Make sure the quaternion is suitable for 3D rotation
    assertQuatOK( q );
#endif

   // Set the quaternion (x,y,z,w) genes
   retval[3] = new RealVector( 4, -1., 1., q.x, q.y, q.z, q.w ); // uniformly-distributed quaternion (UDQ)
   // TODO retval[3] = new ConstrainedRealVector( 4, -1., 1., q.x, q.y, q.z, q.w ); // uniformly-distributed quaternion (UDQ)

   // Set the torsion angles
   retval[4] = new RealVector( num_torsions, -PI, PI );

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R_quaternion()  done assigning each of the retval[0-5] elements...\n\n");
#endif

   return(retval);
}

Genotype generate_Gtype(int num_torsions, GridMapSetInfo *info)
{
#ifdef DEBUG
    // (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Gtype() about to call Genotype temp(5, generate_R())...\n");
    (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Gtype() about to call Genotype temp(5, generate_R_quaternion())...\n");
#endif
   // Genotype temp((unsigned int)5, generate_R(num_torsions, info));
   Genotype temp((unsigned int)5, generate_R_quaternion(num_torsions, info));
#ifdef DEBUG
   // (void)fprintf(logFile,"call_glss.cc/Genotype generate_Gtype() done calling  Genotype temp(5, generate_R())...\n\n");
   (void)fprintf(logFile,"call_glss.cc/Genotype generate_Gtype() done calling  Genotype temp(5, generate_R_quaternion())...\n\n");
#endif

   return(temp);
}

Phenotype generate_Ptype(int num_torsions, GridMapSetInfo *info) 
{
#ifdef DEBUG
    // (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Ptype() about to call Phenotype temp(5, generate_R())...\n");
    (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Ptype() about to call Phenotype temp(5, generate_R_quaternion())...\n");
#endif
   // Phenotype temp((unsigned int)5, generate_R(num_torsions, info));
   Phenotype temp((unsigned int)5, generate_R_quaternion(num_torsions, info));
#ifdef DEBUG
   // (void)fprintf(logFile,"call_glss.cc/Genotype generate_Ptype() done calling  Phenotype temp(5, generate_R())...\n\n");
   (void)fprintf(logFile,"call_glss.cc/Genotype generate_Ptype() done calling  Phenotype temp(5, generate_R_quaternion())...\n\n");
#endif

   return(temp);
}

Individual random_ind(int num_torsions,  GridMapSetInfo *info) 
{

#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/Individual random_ind()  About to generate_Gtype()...\n");
#endif
   Genotype temp_Gtype = generate_Gtype(num_torsions, info);
#ifdef DEBUG
   (void)fprintf(logFile,"call_glss.cc/Individual random_ind()  About to generate_Ptype()...\n");
#endif
   Phenotype temp_Ptype = generate_Ptype(num_torsions, info); 

#ifdef DEBUG
   (void)fprintf(logFile,"call_glss.cc/Individual random_ind()  About to Individual temp(temp_Gtype, temp_Ptype)...\n");
#endif
   //shotgun wedding: does not map genotype to phenotype
   Individual temp(temp_Gtype, temp_Ptype);
   temp.mapping();

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Individual random_ind()  Done     Individual temp(temp_Gtype, temp_Ptype)...\n\n");
#endif

   return(temp);
}

#ifdef FALSE
Individual set_ind(int num_torsions,  GridMapSetInfo *info, State state)
{
   Genotype temp_Gtype;
   Phenotype temp_Ptype;
   Quat q;
   int i;

   temp_Gtype = generate_Gtype(num_torsions, info);
   temp_Ptype = generate_Ptype(num_torsions, info);

   // use the state to generate a Genotype
   temp_Gtype.write(state.T.x, 0);
   temp_Gtype.write(state.T.y, 1);
   temp_Gtype.write(state.T.z, 2);

   q = convertRotToQuat( state.Q );

#ifdef DEBUG_QUAT
#ifdef DEBUG_QUAT_PRINT
    pr( logFile, "DEBUG_QUAT: set_ind()\n" );
    (void) fflush(logFile);
#endif
    //  Make sure the quaternion is suitable for 3D rotation
    assertQuatOK( q );
#endif

   temp_Gtype.write( q.x, 3);
   temp_Gtype.write( q.y, 4);
   temp_Gtype.write( q.z, 5);
   temp_Gtype.write( q.w, 6);

   for (i=0;i<state.ntor; i++) {
       temp_Gtype.write(state.tor[i], 7+i);
   };

   Individual temp(temp_Gtype, temp_Ptype);   

   // use mapping to generate a Phenotype
   //temp.phenotyp =  temp.mapping();
   temp.mapping();

   return(temp);
}
#endif

State call_glss(Global_Search *global_method, Local_Search *local_method, 
                State sInit, 
                unsigned int num_evals, unsigned int pop_size, 
                int outlev, 
                unsigned int extOutputEveryNgens, Molecule *mol, 
                Boole B_RandomTran0, Boole B_RandomQuat0, Boole B_RandomDihe0,
                GridMapSetInfo *info, char *FN_pop_file,
                int end_of_branch[MAX_TORS])
{

#ifdef CALL_GLSS_PROFILING

    Clock  callGlssStart,callGlssSubStart,callGlssSubSubStart,callGlssSubSubSubStart;
    Clock  callGlssEnd;

#endif



    register unsigned int i;
    register int j;
    int num_generations = 0, allEnergiesEqual = 1, numTries = 0;
    int indiv = 0; // Number of Individual in Population to set initial state variables for.
    int max_numTries = 1000;
    double firstEnergy = 0.0;
    EvalMode localEvalMode = Normal_Eval;
    FILE *pop_fileptr;

    global_method->reset(extOutputEveryNgens);
    local_method->reset();
    evaluate[threadID].reset();

#ifdef CALL_GLSS_PROFILING
    callGlssStart = times( &tms_callGlssStart );
#endif

//pkcoff - write to search results buffer
//    (void)fprintf( logFile, "\nCreating an initial population of %u individuals.\n", pop_size);
    sprintf(searchResultsLine[threadID],"\nCreating an initial population of %u individuals.\n", pop_size);
    strcat(searchResultsBuffer[threadID], searchResultsLine[threadID]);

    Population thisPop(pop_size);
    //  Pass in the end_of_branch tree for Branch Crossover Mode.
    thisPop.set_eob( end_of_branch );

    if (sInit.ntor > 0) {
//pkcoff - write to search results buffer
//        (void)fprintf( logFile, "\nAssigning a random translation, a random orientation and %d random torsions to each of the %u individuals.\n\n", sInit.ntor, pop_size);
	    sprintf(searchResultsLine[threadID],"\nAssigning a random translation, a random orientation and %d random torsions to each of the %u individuals.\n\n", sInit.ntor, pop_size);
	    strcat(searchResultsBuffer[threadID], searchResultsLine[threadID]);
    } else {
//pkcoff - write to search results buffer
//        (void)fprintf( logFile, "\nAssigning a random translation and a random orientation to each of the %u individuals.\n\n", pop_size);
	    sprintf(searchResultsLine[threadID],"\nAssigning a random translation and a random orientation to each of the %u individuals.\n\n", pop_size);
	    strcat(searchResultsBuffer[threadID], searchResultsLine[threadID]);
    }
    global_ntor = sInit.ntor; //DEBUG
#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/State call_glss():  {\n");
#endif

#ifdef CALL_GLSS_PROFILING
    callGlssSubStart = times( &tms_callGlssSubStart );
#endif


    do {
        ++numTries;
        // Create a population of pop_size random individuals...
        for (i=0; i<pop_size; i++) {
#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/State call_glss():  Creating individual thisPop[i=%d] using random_ind(%d,info)...\n", i, sInit.ntor);
#endif
            thisPop[i] = random_ind( sInit.ntor, info );
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/State call_glss(): Created  individual i= %d in thisPop[i]\n\n", i);
#endif
            thisPop[i].mol = mol;
            thisPop[i].age = 0L;
        }

        // If initial values were supplied, put them in thisPop[0] and remap
        if (!B_RandomTran0) {
            if (outlev > 1) { (void)fprintf(logFile, "Setting the initial translation (tran0) for individual number %d to %.2lf %.2lf %.2lf\n\n", indiv+1, sInit.T.x, sInit.T.y, sInit.T.z); }
            thisPop[indiv].genotyp.write( sInit.T.x, 0 );
            thisPop[indiv].genotyp.write( sInit.T.y, 1 );
            thisPop[indiv].genotyp.write( sInit.T.z, 2 );
            // Remember to keep the phenotype up-to-date
            thisPop[indiv].mapping();
        }
        if (!B_RandomQuat0) {
            if (outlev > 1) { 
                (void)fprintf(logFile, "Setting the initial orientation using axis-angle values for individual number %d to %.2lf %.2lf %.2lf  %.2lf deg\n\n", indiv+1, sInit.Q.nx, sInit.Q.ny, sInit.Q.nz, RadiansToDegrees(sInit.Q.ang)); 
                (void)fprintf(logFile, "which corresponds to the quaternion (x,y,z,w) values:  %.2lf %.2lf %.2lf %.2lf\n\n", sInit.Q.x, sInit.Q.y, sInit.Q.z, sInit.Q.w); 
            }
            thisPop[indiv].genotyp.write( sInit.Q.x, 3 );
            thisPop[indiv].genotyp.write( sInit.Q.y, 4 );
            thisPop[indiv].genotyp.write( sInit.Q.z, 5 );
            thisPop[indiv].genotyp.write( sInit.Q.w, 6 );
            // Remember to keep the phenotype up-to-date
            thisPop[indiv].mapping();
        }
        if (sInit.ntor > 0) {
            if (!B_RandomDihe0) {
                if (outlev > 1) { (void)fprintf(logFile, "Setting the initial torsions (dihe0) for individual number %d to ", indiv+1); }
                for (j=0; j<sInit.ntor; j++) {
                    thisPop[indiv].genotyp.write( sInit.tor[j], 7+j );
                    if (outlev > 1) { (void)fprintf(logFile, "%.2lf ", RadiansToDegrees(sInit.tor[j])); }
                };
                if (outlev > 1) { (void)fprintf(logFile, " deg\n\n"); }
                // Remember to keep the phenotype up-to-date
                thisPop[indiv].mapping();
            }
        }

#ifdef DEBUG
    (void)fprintf(logFile,"\n\ncall_glss.cc  // ensuring there is variation in the energies\n\n");
#endif
        // Now ensure that there is some variation in the energies...
        firstEnergy = thisPop[0].value(localEvalMode);
#ifdef DEBUG
    (void)fprintf(logFile,"\n\ncall_glss.cc  // ensuring there is variation in the energies, firstEnergy=%lf\n\n", firstEnergy);
#endif
        for (i=1; i<pop_size; i++) {
#ifdef DEBUG
    (void)fprintf(logFile,"\n\ncall_glss.cc  // ensuring there is variation in the energies, i=%d, thisPop[i].value=%lf\n\n", i, thisPop[i].value(localEvalMode));
#endif
             allEnergiesEqual = allEnergiesEqual && (thisPop[i].value(localEvalMode) == firstEnergy);
        }
        if ( pop_size>1 && allEnergiesEqual) {
            (void)fprintf(logFile,"NOTE: All energies are equal in population; re-initializing. (Try Number %d)\n", numTries);
        }
        if (numTries > max_numTries) {
            (void)fprintf(logFile,"WARNING: the number of tries has exceeded the maximum number of tries permitted.\nWARNING: AutoDock will attempt continue with the currently-generated random population.\n");
            break;
        }
    } while (pop_size>1 && allEnergiesEqual);


#ifdef CALL_GLSS_PROFILING

    callGlssEnd = times( &tms_callGlssEnd );
    pr( logFile, "\nCALL_GLSS_PROFILING: population creation for call_glss: ");
    timesyshms( callGlssEnd - callGlssSubStart, &tms_callGlssSubStart, &tms_callGlssEnd );

#endif

#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/State call_glss():  }\n");
    if (outlev > 2) { 
    thisPop.printPopulationAsCoordsEnergies( logFile, pop_size, sInit.ntor );
    }
#endif

    if (outlev > 2) { 
        (void)fprintf( logFile, "The initial population consists of the following %d individuals:\n\n", pop_size);
        (void)fprintf( logFile, "<generation t=\"%d\" after_performing=\"initialisation of population\">\n", num_generations);
        (void)fprintf( logFile, "</generation>\n\n\n");
    }

    if (pop_size > 1 && outlev > 3 ) { minmeanmax( logFile, thisPop, num_generations, info ); }

//We now have a mapped and evaluated population suitable for global search

//pkcoff - write to search results buffer
//        (void)fprintf( logFile, "Beginning Lamarckian Genetic Algorithm (LGA), with a maximum of %u\nenergy evaluations.\n\n", num_evals);
	    sprintf(searchResultsLine[threadID],"Beginning Lamarckian Genetic Algorithm (LGA), with a maximum of %u\nenergy evaluations.\n\n", num_evals);
	    strcat(searchResultsBuffer[threadID], searchResultsLine[threadID]);
    
#ifdef CALL_GLSS_PROFILING
    callGlssSubStart = times( &tms_callGlssSubStart );

#endif


    do {
        ++num_generations;

        if (outlev > 1) { (void)fprintf( logFile, "Global-Local Search Iteration: %d\n", num_generations); }
        
        if (outlev > 1) { (void)fprintf( logFile, "Performing Global Search.\n"); }

#ifdef DEBUG
    (void)fprintf(logFile,"\n\ncall_glss.cc do loop global search for thread ID %d\n",threadID);
#endif

        global_method->search(thisPop);

        if (outlev > 2) {
            (void)fprintf( logFile, "<generation t=\"%d\" after_performing=\"global search\">\n", num_generations);
            thisPop.printPopulationAsStates( logFile, pop_size, sInit.ntor );
            (void)fprintf( logFile, "</generation>\n\n\n");
        }

        if (pop_size > 1 && outlev > 3) { minmeanmax( logFile, thisPop, num_generations, info ); }

        if (outlev > 1) { (void)fprintf( logFile, "Performing Local Search.\n"); }

#ifdef CALL_GLSS_PROFILING
    callGlssSubSubStart = times( &tms_callGlssSubSubStart );

#endif

        for (i=0; i<pop_size; i++) {
            if (outlev > 1) {
                (void)fprintf( logFile, "LS: %d",num_generations); 
                (void)fprintf( logFile, " %d",i+1); 
                (void)fprintf( logFile, " %f",thisPop[i].value(localEvalMode)); 
            }

#ifdef DEBUG
    (void)fprintf(logFile,"\n call_glss.cc do loop local search for threadID %d\n",threadID);
#endif

            local_method->search(thisPop[i]);

#ifdef DEBUG
    (void)fprintf(logFile,"\n AFTER call_glss.cc do loop local search for thread id %d\n",threadID);
#endif


            if (outlev > 1) {
                (void)fprintf( logFile, " %f",thisPop[i].value(localEvalMode)); 
                (void)fprintf( logFile, " \n"); 
            }
        }

#ifdef CALL_GLSS_PROFILING

    callGlssEnd = times( &tms_callGlssEnd );
    pr( logFile, "\nCALL_GLSS_PROFILING: local loop search for call_glss: ");
    timesyshms( callGlssEnd - callGlssSubSubStart, &tms_callGlssSubSubStart, &tms_callGlssEnd );

#endif

        if (outlev > 2) {
            (void)fprintf( logFile, "<generation t=\"%d\" after_performing=\"local search\">\n", num_generations);
            thisPop.printPopulationAsStates( logFile, pop_size, sInit.ntor );
            (void)fprintf( logFile, "</generation>\n\n\n");
        }

        if (pop_size > 1 && outlev > 3) { minmeanmax( logFile, thisPop, num_generations, info ); }

        if (strcmp (FN_pop_file, "") != 0) { // YES, do print!
            if ((pop_fileptr = ad_fopen( FN_pop_file, "w")) == NULL) {
                pr(logFile, "\n%s: ERROR:  I'm sorry, I cannot create\"%s\".\n\n", programname, FN_pop_file);
            } else {
                thisPop.printPopulationAsCoordsEnergies( pop_fileptr, pop_size, sInit.ntor); 
                fclose( pop_fileptr );
            }
        }

#ifdef DEBUG
    (void)fprintf(logFile,"\n\ncall_glss.cc do loop END for threadid %d\n",threadID);
#endif

        //(void)fflush(logFile);
    } while ((evaluate[threadID].evals() < num_evals) && (!global_method->terminate()));


#ifdef CALL_GLSS_PROFILING

    callGlssEnd = times( &tms_callGlssEnd );
    pr( logFile, "\nCALL_GLSS_PROFILING: global search for call_glss: ");
    timesyshms( callGlssEnd - callGlssSubStart, &tms_callGlssSubStart, &tms_callGlssEnd );

#endif

    thisPop.msort(1);
//pkcoff - write to search results buffer
//        (void)fprintf(logFile,"Final-Value: %.3f\n", thisPop[0].value(Normal_Eval));
	    sprintf(searchResultsLine[threadID],"Final-Value: %.3f\n", thisPop[0].value(Normal_Eval));
	    strcat(searchResultsBuffer[threadID], searchResultsLine[threadID]);
    
#ifdef CALL_GLSS_PROFILING

    callGlssEnd = times( &tms_callGlssEnd );
    pr( logFile, "\nCALL_GLSS_PROFILING: total time taken for call_glss: ");
    timesyshms( callGlssEnd - callGlssStart, &tms_callGlssStart, &tms_callGlssEnd );
#endif

    return( thisPop[0].state(sInit.ntor) );
}
