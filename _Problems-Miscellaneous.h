/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: BIP (Block-lterative Projection) method (No MPI)
Module: _Problems-Miscellaneous.h (Miscellaneous LP problems)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
==============================================================================*/
#pragma once

//#define PP_COS_MODE

//=========================== problem Parameters ========================
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO	// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO)	// Accuracy of belonging to hyperplane
#define PP_EPS_COS					1E-2	// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO/100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO/100)	// Precision for point to be inside cone

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

/*============================== cubeInHyperplane LP problem ===================*
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"cubeInHyperplane"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 4		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		60000
//-------------------------------------------------------------------------------

/*============================== simple1FxVar LP problem ========================*
// Simple LP problem & x_1=150
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"simple1FxVar"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 3		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 52500
//-------------------------------------------------------------------------------

/*============================== simple1 LP problem =============================*
#define PP_PROBLEM_NAME	"simple1"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 7		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		55000
//-------------------------------------------------------------------------------

/*============================== simple1.1 LP problem ===========================*
// Simple LP problem with alternating objective function
#define PP_PROBLEM_NAME	"simple1.1"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 7		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		40000
//-------------------------------------------------------------------------------

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

/*============================== simple1min LP problem ==========================*
#define PP_PROBLEM_NAME	"simple1min"
#define PP_M 5		// Number of rows in *.mtx
#define PP_N 8		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		-5000
//-------------------------------------------------------------------------------

/*============================== simple_zcv LP problem ==========================*
#define PP_PROBLEM_NAME	"simple_zcv"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 7		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		50000
//-------------------------------------------------------------------------------

/*============================== simple_lcv LP problem ==========================*
#define PP_PROBLEM_NAME	"simple_lcv"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 7		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		50000.2
//-------------------------------------------------------------------------------

/*============================== simple_lcv_neg LP problem ======================*
#define PP_PROBLEM_NAME	"simple_lcv_neg"
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 7		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		49998
//-------------------------------------------------------------------------------

/*============================== angle03 LP problem =============================*
#define PP_PROBLEM_NAME	"angle03"
#define PP_M 3		// Number of rows in *.mtx
#define PP_N 6		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		3000
//-------------------------------------------------------------------------------

/*============================== angle04 LP problem =============================*
#define PP_PROBLEM_NAME	"angle04"
#define PP_M 3		// Number of rows in *.mtx
#define PP_N 7		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		3300
//-------------------------------------------------------------------------------

/*============================== cone3-0 LP problem =============================*
#define PP_PROBLEM_NAME	"cone3-0"
#define PP_M 11		// Number of rows in *.mtx
#define PP_N 14		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		132.5
//-------------------------------------------------------------------------------

/*============================== rnd3-10 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd3-10"
#define PP_M 13		// Number of rows in *.mtx
#define PP_N 16		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE 		852.0289179009677
//------------------------------------------------------------------------------

/*============================== rnd5-100 LP problem ===========================*/
#define PP_PROBLEM_NAME	"rnd5-100"
#define PP_M 105		// Number of rows in *.mtx
#define PP_N 110		// Number of cols in *.mtx
#define PP_MAX_OBJ_VALUE	1848.437080587865
//------------------------------------------------------------------------------

/*==============================================================================*/