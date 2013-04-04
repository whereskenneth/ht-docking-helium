/*

 $Id: autoglobal.h,v 1.4 2010/08/03 22:31:05 pkcoff Exp $

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

#ifndef _AUTOGLOBAL
#define _AUTOGLOBAL

#include <sys/types.h>
#include <string.h>
#include <stdio.h>

#include "structs.h"

/******************************************************************************/
/*      Name: autoglobal.h                                                    */
/*  Function: Global variables for Autodock modules.                          */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett M. Morris                                               */
/*                                                                            */
/*            e-mail: garrett@scripps.edu				      */
/*                                                                            */
/*            The Scripps Research Institute                                  */
/*            Department of Molecular Biology, MB5                            */
/*            10666 North Torrey Pines Road                                   */
/*            La Jolla, CA 92037.                                             */
/*                                                                            */
/*      Date: 03/18/93                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: None.                                                           */
/*   Returns: Nothing.                                                        */
/*   Globals: programname, AutoDockHelp, command_mode,          */
/*            command_in_fp, command_out_fp, GPF, logFile.                    */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 03/18/93 GMM     Created.                                                  */
/******************************************************************************/


/*----------------------------------------------------------------------------*/
/* Globals,                                                                   */
/*----------------------------------------------------------------------------*/

char    *programname;

char    dock_param_fn[PATH_MAX];
char    grid_param_fn[PATH_MAX];

int     command_mode = FALSE;
int     debug = 0;
int	    ElecMap = 0;
int	    DesolvMap = 0;
int     ignore_errors = FALSE;
int     keepresnum = 1;
int     parse_tors_mode = FALSE;
int	    true_ligand_atoms = 0;
int     write_stateFile = FALSE;
// For energy breakdown of non-bonded interactions
int     Nnb_array[3] = {0};    // number of nonbonds in the ligand, intermolecular and receptor groups

Real	idct = 1.0;
// For energy breakdown of non-bonded interactions
Real    nb_group_energy[3] = {0.0};  // total energy of each nonbond group (intra-ligand, inter, and intra-receptor)

FILE    *command_in_fp;
FILE    *command_out_fp;
FILE    *parFile;
FILE    *GPF;
FILE    *logFile;
FILE    *stateFile;

Linear_FE_Model AD3;
Linear_FE_Model AD4_wrt_3;
Linear_FE_Model AD4;

Unbound_Model ad4_unbound_model = Unbound_Default;


//pkcoff - openmp for bluegene change section - start
// global_ntor moved from call_glss.cc to here and made threadprivate
int global_ntor;

// evaluate won't link as threadprivate, so array has been created.  evaluate also moved from call_glss.cc to here
Eval evaluate[MAX_NUM_THREADS];

// threadid is threadprivate and set in mpi_slave_main function, used to index evaluate object array
int threadID = -1;

// numThreads is threadprivate and set in mpi_slave_main function
int numThreads = 0;

// search objects and rho arrays need to be threadprivate and created in each individual thread in mpi_slave_main function
Global_Search *GlobalSearchMethod = NULL;
Local_Search *LocalSearchMethod = NULL;
Real *rho_ptr = NULL;
Real *lb_rho_ptr = NULL;

#pragma omp threadprivate (threadID, global_ntor,GlobalSearchMethod, LocalSearchMethod, rho_ptr, lb_rho_ptr, nb_group_energy, Nnb_array)

// Buffer for writing the search results during the parallelization sections
char searchResultsBuffer[MAX_NUM_THREADS][DLG_BUFFER_LENGTH];

// temp space for writing 1 line for the search results buffer
char searchResultsLine[MAX_NUM_THREADS][DLG_LINE_LENGTH];

//pkcoff - openmp for bluegene change section - end

// pkcoff: support for unduplicated file-loading
// int num_maps = 0;
int num_maps = 0;
char loadedMapFileNames[MAX_MAPS][FILENAME_MAX] = {'\0'};
char atomTypeMap[MAX_MAPS][3] = {'\0'};
char loadedMapFileNamesCopy[FILENAME_MAX] = {'\0'};
int num_maps_loaded = 0;
GridMapSetInfo *info;  // this information is from the AVS field file
//map has to be threadprivate
MapType map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS];
MapType mapCopyOneRow[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS];
char atomTypeMapCopy[3] = {'\0'};
char dockingDirectoryBaseName[FILENAME_MAX];

//#pragma omp threadprivate(map)
// pkcoff: support for unduplicated file-loading end


#endif /*_AUTOGLOBAL*/
/* EOF */
