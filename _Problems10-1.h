/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: FIP (Feasible point Iterative Projection) method (No MPI)
Module: Problems10-1.h (LP problems of dimension 10 with 1 randome inequality)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
LP problems were obtained using BSF-LPP-Generator.
================================================================================*/
#pragma once

//=========================== Method Parameters ================================
#define PP_M 11		// Number of equations (number of rows in *.mtx)
#define PP_N 21		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------------------

/*============================== rnd10-0 LP problem ============================*/
#define PP_PROBLEM_NAME	"rnd10-0"
#define PP_MAX_OBJ_VALUE 10900			
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 2
// Objective function: 1000
// Distance to polytope: 0
// z on hyperplanes: {4, 11, 12, 13, 14, 16, 17, 18, 19, 20}
// Number of including inequality hyperplanes: 10
//-------------------------------------------------------------------------------

/*============================== rnd10-1-1 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd10-1-1"


//-------------------------------------------------------------------------------

/*============================== rnd10-1-2 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd10-1-2"


//-------------------------------------------------------------------------------

/*============================== rnd10-1-3 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd10-1-3"


//-------------------------------------------------------------------------------

/*============================== rnd10-1-4 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd10-1-4"


//-------------------------------------------------------------------------------

/*============================== rnd10-1-5 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd10-1-5"


//-------------------------------------------------------------------------------

/*============================== rnd10-1-6 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd10-1-6"


//-------------------------------------------------------------------------------
/*===============================================================================*/