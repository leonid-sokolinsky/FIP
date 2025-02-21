/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: BIP (Block-lterative Projection) method (No MPI)
Module: _Problems-NetLib-LP.h (Problems from the NETLIB LP Test Problem Set)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
==============================================================================*/
#pragma once

#define PP_MPS_FORMAT

/*============================== adlittle LP problem =======================*
// Number of equations : 15
// Subspace dimension : 82
#define PP_PROBLEM_NAME		"adlittle"
#define PP_M 56	// Number of constraints in mps-file
#define PP_N 97	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 		-225494.96316238038228101176621492
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11		// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	PP_EPS_ZERO	// Precision for point to be in halfspace
#define PP_EPS_PROJECTION			PP_EPS_ZERO	// Precision to calculate projection
//--------------------------------------------------------------------------
// Elapsed time: 3
// Number of iterations: 96889
// Number of including inequality hyperplanes: 34
// Objective function: -652787
//--------------------------------------------------------------------------

/*============================== afiro LP problem ==========================*
// Number of equations : 8
// Subspace dimension : 24
#define PP_PROBLEM_NAME	"afiro"
#define PP_M 27		// Number of constraints in mps-file
#define PP_N 32		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 464.75314285714285714285714285714
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11		// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	PP_EPS_ZERO	// Precision for point to be in halfspace
#define PP_EPS_PROJECTION			PP_EPS_ZERO	// Precision to calculate projection
//--------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 13237
// Number of including inequality hyperplanes: 23
// Objective function: -153.138
//--------------------------------------------------------------------------

/*============================== beaconfd LP problem =======================*
// Number of equations: 140
// Subspace dimension: 122
#define PP_PROBLEM_NAME		"beaconfd"
#define PP_M 173	// Number of constraints in mps-file
#define PP_N 262	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE -33592.4858072
//---------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11		// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	PP_EPS_ZERO	// Precision for point to be in halfspace
#define PP_EPS_PROJECTION			PP_EPS_ZERO	// Precision to calculate projection
//---------------------------------------------------------------------
// Elapsed time: 104
// Number of iterations: 548860
// Number of including inequality hyperplanes: 125
// Objective function: -33854.64
//---------------------------------------------------------------------

/*============================== blend LP problem =====================*
// Number of equations: 43
// Subspace dimension: 40
#define PP_PROBLEM_NAME		"blend"
#define PP_M 74	// Number of constraints in mps-file
#define PP_N 83		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 30.812149845828220173774356124984	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-11			// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_PROJECTION		PP_EPS_ZERO		// Precision to calculate projection
//--------------------------------------------------------------------------

/*============================== fit1d LP problem ==========================*
// Number of equations : 1
// Subspace dimension : 1025
#define PP_PROBLEM_NAME		"fit1d"
#define PP_M 24	// Number of equations (after conversion to standard form)
#define PP_N 1026	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 9146.3780924209269467749025024617	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11		// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	PP_EPS_ZERO	// Precision for point to be in halfspace
#define PP_EPS_PROJECTION			PP_EPS_ZERO	// Precision to calculate projection
//--------------------------------------------------------------------------
// Zero point is feasible
//--------------------------------------------------------------------------

/*============================== israel LP problem =========================*
// Number of equations: 0
#define PP_PROBLEM_NAME		"israel"
#define PP_M 174	// Number of constraints in mps-file
#define PP_N 142	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 896644.82186304572966200464196045	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11		// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	PP_EPS_ZERO	// Precision for point to be in halfspace
#define PP_EPS_PROJECTION			PP_EPS_ZERO	// Precision to calculate projection
//--------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 143
// Number of including inequality hyperplanes: 141
// Objective function: 85567.12
//--------------------------------------------------------------------------

/*============================== kb2 LP problem ============================*/
// Number of equations: 16
// Subspace dimension: 25
#define PP_PROBLEM_NAME		"kb2"
#define PP_M 43	// Number of equations (after conversion to standard form)
#define PP_N 41		// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 1749.9001299062057129526866493726
//--------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-12			// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_ZERO*100)// Precision for point to be in halfspace
#define PP_EPS_PROJECTION		PP_EPS_ZERO		// Precision to calculate projection
//--------------------------------------------------------------------------
// Zero point is feasible
//--------------------------------------------------------------------------

/*============================== recipe LP problem =========================*
// Number of equations: 79
// Subspace dimension: 101
#define PP_PROBLEM_NAME		"recipe"
#define PP_M 91	// Number of constraints in mps-file
#define PP_N 180	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 266.616 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11		// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	PP_EPS_ZERO	// Precision for point to be in halfspace
#define PP_EPS_PROJECTION			PP_EPS_ZERO	// Precision to calculate projection
//--------------------------------------------------------------------------
// Elapsed time: 2
// Number of iterations: 13602
// Number of including inequality hyperplanes: 190
// Objective function: 104.8237
//--------------------------------------------------------------------------

/*============================== sc50a LP problem ==========================*
// Number of equations: 20
// Subspace dimension: 28
#define PP_PROBLEM_NAME		"sc50a"
#define PP_M 49	// Number of constraints in mps-file
#define PP_N 48		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 64.575077058564509026860413914575	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11		// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	PP_EPS_ZERO	// Precision for point to be in halfspace
#define PP_EPS_PROJECTION			PP_EPS_ZERO	// Precision to calculate projection
//----------------------------------------------------------------------------
// Zero point is feasible
//----------------------------------------------------------------------------

/*============================== sc50b LP problem ==============================*
#define PP_PROBLEM_NAME		"sc50b"
#define PP_M 48	// Number of constraints in mps-file
#define PP_N 48		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 70	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		1E-7	// Accuracy of belonging to polytope
//----------------------------------------------------------------------------
// Zero point is feasible
//----------------------------------------------------------------------------

/*============================== sc105 LP problem ==============================*
#define PP_PROBLEM_NAME		"sc105"
#define PP_M 104	// Number of constraints in mps-file
#define PP_N 103	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 52.202061211707248062628010857689 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		1E-7	// Accuracy of belonging to polytope
//----------------------------------------------------------------------------
// Zero point is feasible.
//----------------------------------------------------------------------------

/*============================== share2b LP problem ==========================*
#define PP_PROBLEM_NAME		"share2b"	
#define PP_M 96	// Number of constraints in mps-file
#define PP_N 79	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 415.732240741419486545199108738 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11		// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	PP_EPS_ZERO	// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		PP_EPS_ZERO	// Accuracy of belonging to polytope
//----------------------------------------------------------------------------
// Elapsed time: 1550
// Number of iterations: 79885851
// Number of including inequality hyperplanes: 39
// Objective function: 354.6537
//----------------------------------------------------------------------------

//============================================================================/**/