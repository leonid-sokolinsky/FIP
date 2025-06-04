/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: FIP (Feasible Iterative Projection) method (No MPI)
Module: Problems05-1.h (LP problems of dimension 5 with 1 randome inequality)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
LP problems were obtained using BSF-LPP-Generator.
=================================================================================*/
#pragma once

//=========================== Problem Parameters ===============================
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 11		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------------------

/*============================== rnd5-0 LP problem =============================*/
// Exact solution:	100   200   200   200   200
#define PP_PROBLEM_NAME	"rnd5-0"
#define PP_MAX_OBJ_VALUE 2900
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 119
// Objective function: 2346.64317630873756570509
// Distance to polytope: 1.42108547152020037174225e-13
// z on hyperplanes: {0, 1, 2, 3}
// Number of including inequality hyperplanes: 4
//-------------------------------------------------------------------------------

/*============================== rnd5-1-1 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd5-1-1"
#define PP_MAX_OBJ_VALUE 2584.349489711176
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 31
// Objective function: 2179.55445184212794629275
// Distance to polytope: 0
// z on hyperplanes: {0}
// Number of including inequality hyperplanes: 1
//-------------------------------------------------------------------------------

/*============================== rnd5-1-2 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd5-1-2"
#define PP_MAX_OBJ_VALUE 2657.525612546482
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 60
// Objective function: 2087.12102390933250717353
// Distance to polytope: 0
// z on hyperplanes: {1, 5}
// Number of including inequality hyperplanes: 2
//-------------------------------------------------------------------------------

/*============================== rnd5-1-3 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd5-1-3"
#define PP_MAX_OBJ_VALUE 2424.919153814418
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

/*============================== rnd5-1-4 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd5-1-4"
#define PP_MAX_OBJ_VALUE 2300.112758702127
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

/*============================== rnd5-1-5 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd5-1-5"
#define PP_MAX_OBJ_VALUE 2626.473236211187
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

/*============================== rnd5-1-6 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd5-1-6"
#define PP_MAX_OBJ_VALUE 2675.351994188833
#define PP_EPS_ZERO				1E-15					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

/*===============================================================================*/