/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: BIP (Block-lterative Projection) method (No MPI)
Module: _Problems-Miscellaneous.h (Miscellaneous LP problems)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
================================================================================*/
#pragma once

/*
//=========================== problem Parameters ===============================
#define PP_EPS_ZERO					1E-12		// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_ZERO*10)// Precision for point to be in halfspace
#define PP_EPS_PROJECTION		PP_EPS_ZERO		// Accuracy of belonging to hyperplane
*/

/*============================== nguyen5 LP problem ============================*/
#define PP_PROBLEM_NAME	"nguyen5"
#define PP_MPS_FORMAT
#define PP_M 4		// Number of constrains
#define PP_N 5		// Number of variables
#define PP_MAX_OBJ_VALUE 21.4549732313109693393471
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 75
// Objective function: 5.08235294117632818711172
// Distance to polytope: 2.35117716511801969942853e-15
// z on hyperplanes: {2, 4, 5, 7, 8}
// Number of including inequality hyperplanes: 5
//------------------------------------------------------------------------------

/*============================== simpleCube LP problem ==========================*
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"simpleCube"
#ifdef PP_MPS_FORMAT
#define PP_M 3		// Number of constraints
#define PP_N 3		// Number of variables
#else
#define PP_M 3	// Number of rows in *.mtx
#define PP_N 6	// Number of cols in *.mtx
#endif
#define PP_MAX_OBJ_VALUE 		60000
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 4
// Objective function: 60000
// Distance to polytope: 0
// z on hyperplanes: {0, 1, 2}
// Number of including inequality hyperplanes: 3
//------------------------------------------------------------------------------

/*============================== cubeInHyperplane LP problem ===================*
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"cubeInHyperplane"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 4		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		90000
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 127
// Objective function: 42000
// Distance to polytope: 0
// z on hyperplanes: {0, 3}
// Number of including inequality hyperplanes: 1
//------------------------------------------------------------------------------

/*============================== simple1FxVar LP problem ========================*
// Simple LP problem & x_1=150
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"simple1FxVar"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 3		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 52500
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 3
// Objective function: 15000
// Distance to polytope: 0
// z on hyperplanes: {4, 6, 7}
// Number of including inequality hyperplanes: 2
//------------------------------------------------------------------------------

/*============================== simple1 LP problem =============================*
#define PP_PROBLEM_NAME	"simple1"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 7		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		55000
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 5
// Objective function: 10000
// Distance to polytope: 0
// z on hyperplanes: {2, 4, 5}
// Number of including inequality hyperplanes: 3
//------------------------------------------------------------------------------

/*============================== simple1.1 LP problem ===========================*
// Simple LP problem with alternating objective function
#define PP_PROBLEM_NAME	"simple1.1"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 7		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		40000
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 7
// Objective function: 17984.1820987654318741988
// Distance to polytope: 0
// z on hyperplanes: {3}
// Number of including inequality hyperplanes: 1
//------------------------------------------------------------------------------

/*============================== simple2 LP problem =============================*
// Simple LP problem & x_3=200; x_2>=110; x_0<=190
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"simple2"
#ifdef PP_MPS_FORMAT
#define PP_M 5		// Number of constraints
#define PP_N 4		// Number of variables
#else
#define PP_M 5		// Number of rows in *.mtx
#define PP_N 8		// Nnumber of cols in *.mtx
#endif
#define PP_MAX_OBJ_VALUE 	63500
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 55
// Objective function: 10500
// Distance to polytope: 0
// z on hyperplanes: {4, 5, 6, 7}
// Number of including inequality hyperplanes: 3
//------------------------------------------------------------------------------

/*============================== simple3 LP problem =============================*
//#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"simple3"
#ifdef PP_MPS_FORMAT
#define PP_M 5		// Number of constraints
#define PP_N 5		// Number of variables
#else
#define PP_M 5		// Number of rows in *.mtx
#define PP_N 8		// Nnumber of cols in *.mtx
#endif
#define PP_MAX_OBJ_VALUE 		55000
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 621
// Objective function: 24999.9999999999563442543
// Distance to polytope: 1.00485917355761594892983e-13
// z on hyperplanes: {1, 7}
// Number of including inequality hyperplanes: 1
//------------------------------------------------------------------------------

/*============================== simple1min LP problem ==========================*
#define PP_PROBLEM_NAME	"simple1min"
#define PP_M 5		// Number of rows in *.mtx
#define PP_N 8		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		-5000
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 2
// Objective function: -10000
// Distance to polytope: 0
// z on hyperplanes: {4}
// Number of including inequality hyperplanes: 1
//------------------------------------------------------------------------------

/*============================== simple_zcv LP problem ==========================*
#define PP_PROBLEM_NAME	"simple_zcv"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 7		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		50000
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 5
// Objective function: 45000
// Distance to polytope: 0
// z on hyperplanes: {1, 2}
// Number of including inequality hyperplanes: 2
//------------------------------------------------------------------------------

/*============================== simple_lcv LP problem ==========================*
#define PP_PROBLEM_NAME	"simple_lcv"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 7		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		50000.2
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 5
// Objective function: 47777.9777777777781011537
// Distance to polytope: 0
// z on hyperplanes: {0, 2}
// Number of including inequality hyperplanes: 2
//------------------------------------------------------------------------------

/*============================== simple_lcv_neg LP problem ======================*
#define PP_PROBLEM_NAME	"simple_lcv_neg"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 7		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		49998
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 5
// Objective function: 43323.2565157750286743976
// Distance to polytope: 0
// z on hyperplanes: {3}
// Number of including inequality hyperplanes: 1
//------------------------------------------------------------------------------

/*============================== angle03 LP problem =============================*
#define PP_PROBLEM_NAME	"angle03"
#define PP_M 3		// Number of rows in *.mtx
#define PP_N 6		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		3000
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 18
// Objective function: 1100
// Distance to polytope: 0
// z on hyperplanes: {0, 3, 5}
// Number of including inequality hyperplanes: 3
//------------------------------------------------------------------------------

/*============================== angle04 LP problem =============================*
#define PP_PROBLEM_NAME	"angle04"
#define PP_M 3		// Number of rows in *.mtx
#define PP_N 7		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		3300
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 106
// Objective function: 113.99999999999951683094
// Distance to polytope: 1.27162277778984453168361e-14
// z on hyperplanes: {0, 3, 5, 6}
// Number of including inequality hyperplanes: 4
//------------------------------------------------------------------------------

/*============================== cone3-0 LP problem =============================*
#define PP_PROBLEM_NAME	"cone3-0"
#define PP_M 11		// Number of rows in *.mtx
#define PP_N 14		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		132.5
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 128
// Objective function: 51.2430710132439557469297
// Distance to polytope: 1.09893689810250844314941e-14
// z on hyperplanes: {5, 12}
// Number of including inequality hyperplanes: 2
//------------------------------------------------------------------------------

/*============================== rnd3-10 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd3-10"
#define PP_M 13		// Number of rows in *.mtx
#define PP_N 16		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		852.0289179009677
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 3
// Objective function: 783.708298413821353278763
// Distance to polytope: 0
// z on hyperplanes: {3}
// Number of including inequality hyperplanes: 1
//------------------------------------------------------------------------------

/*============================== rnd5-100 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd5-100"
#define PP_M 105		// Number of rows in *.mtx
#define PP_N 110		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE	1848.437080587865
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 4
// Objective function: 1436.47875368017207620142
// Distance to polytope: 0
// z on hyperplanes: {88}
// Number of including inequality hyperplanes: 1
//------------------------------------------------------------------------------

/*==============================================================================*/