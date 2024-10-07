/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: BIP (Block-lterative Projection) method (No MPI)
Module: Problems20-1.h (LP problems of dimension 20 with 1 randome inequality)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
LP problems were obtained using BSF-LPP-Generator.
==============================================================================*/
#pragma once

//=========================== Method Parameters ===============================
#define PP_EPS_ZERO					1E-11		// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO	// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		PP_EPS_ZERO	// Accuracy of belonging to hyperplane
#define PP_EPS_COS					PP_EPS_ZERO	// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO /100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO /100)	// Precision for point to be inside cone
//-------------------------------------------------------------------------------
#define PP_M 21		// Number of equations (number of rows in *.mtx)
#define PP_N 41		// Number of variables (number of cols in *.mtx)

/*============================== rnd20-0 LP problem ==============================*/
#define PP_PROBLEM_NAME	"rnd20-0"
#define PP_MAX_OBJ_VALUE 41900
//------------------------------------------------------------------




/*============================== rnd20-1-1 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd20-1-1"


//------------------------------------------------------------------



//------------------------------------------------------------------/**/

/*============================== rnd20-1-2 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd20-1-2"


//------------------------------------------------------------------



//------------------------------------------------------------------/**/

/*============================== rnd20-1-3 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd20-1-3"


//------------------------------------------------------------------



//------------------------------------------------------------------/**/

/*============================== rnd20-1-4 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd20-1-4"


//------------------------------------------------------------------



//------------------------------------------------------------------/**/

/*============================== rnd20-1-5 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd20-1-5"


//------------------------------------------------------------------



//------------------------------------------------------------------/**/