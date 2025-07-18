/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: FIP (Feasible point Iterative Projection) method (No MPI)
Module: _Problems-NetLib-LP.h (Problems from the NETLIB LP Test Problem Set)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
================================================================================*/
#pragma once

#define PP_MPS_FORMAT

/*============================== adlittle LP problem ===========================*/
// Number of equations : 15
// Subspace dimension : 82
#define PP_PROBLEM_NAME		"adlittle"
#define PP_M 56	// Number of constraints in mps-file
#define PP_N 97	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 		-225494.96316238038228101176621492
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 13202
// Objective function: -580172.408514236100018024
// Distance to polytope: 9.92761536223621320227147e-13
// Number of including inequality hyperplanes: 40
//------------------------------------------------------------------------------

/*============================== afiro LP problem ==============================*
// Number of equations : 8
// Subspace dimension : 24
#define PP_PROBLEM_NAME	"afiro"
#define PP_M 27		// Number of constraints in mps-file
#define PP_N 32		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 464.75314285714285714285714285714
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
/// Elapsed time: 0
// Number of iterations: 4678
// Objective function: -153.078793414766835212504
// Distance to polytope: 9.95279125405900389970587e-13
// Number of including inequality hyperplanes: 24
//------------------------------------------------------------------------------

/*============================== beaconfd LP problem ===========================*
// Number of equations: 140
// Subspace dimension: 122
#define PP_PROBLEM_NAME		"beaconfd"
#define PP_M 173	// Number of constraints in mps-file
#define PP_N 262	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE -33592.4858072
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 10
// Number of iterations: 66253
// Objective function: -33856.9000641296297544613
// Distance to polytope: 9.97561425517768269030821e-13
// Number of including inequality hyperplanes: 135
//------------------------------------------------------------------------------

/*============================== blend LP problem ==============================*
// Number of equations: 43
// Subspace dimension: 40
#define PP_PROBLEM_NAME		"blend"
#define PP_M 74	// Number of constraints in mps-file
#define PP_N 83		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 30.812149845828220173774356124984	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Zero point is feasible
//------------------------------------------------------------------------------

/*============================== fit1d LP problem ==============================*
// Number of equations : 1
// Subspace dimension : 1025
#define PP_PROBLEM_NAME		"fit1d"
#define PP_M 24	// Number of equations (after conversion to standard form)
#define PP_N 1026	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 9146.3780924209269467749025024617	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Zero point is feasible
//------------------------------------------------------------------------------

/*============================== grow7 LP problem ==============================*
// Number of equations: 140
// Subspace dimension: 161
#define PP_PROBLEM_NAME		"grow7"
#define PP_M 140	// Number of equations (after conversion to standard form)
#define PP_N 301	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 47787811.814711502616766956242865	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Zero point is feasible
//------------------------------------------------------------------------------

/*============================== israel LP problem =============================*
// Number of equations: 0
// Subspace dimension: 142
#define PP_PROBLEM_NAME		"israel"
#define PP_M 174	// Number of constraints in mps-file
#define PP_N 142	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 896644.82186304572966200464196045	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	1E-12	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 9
// Objective function: 84416.6666666666715173051
// Distance to polytope: 0
// Number of including inequality hyperplanes: 141
//------------------------------------------------------------------------------

/*============================== kb2 LP problem ================================*
// Number of equations: 16
// Subspace dimension: 25
#define PP_PROBLEM_NAME		"kb2"
#define PP_M 43	// Number of equations (after conversion to standard form)
#define PP_N 41		// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 1749.9001299062057129526866493726
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Zero point is feasible
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

/*============================== recipe LP problem =============================*
// Number of equations: 79
// Subspace dimension: 101
#define PP_PROBLEM_NAME		"recipe"
#define PP_M 91	// Number of constraints in mps-file
#define PP_N 180	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 266.616 // Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 1
// Number of iterations: 4616
// Objective function: 104.823166666645093414445
// Distance to polytope: 9.34363697524531744420528e-13
// Number of including inequality hyperplanes: 190
//------------------------------------------------------------------------------

/*============================== sc105 LP problem ==============================*
#define PP_PROBLEM_NAME		"sc105"
#define PP_M 104	// Number of constraints in mps-file
#define PP_N 103	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 52.202061211707248062628010857689 // Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Zero point is feasible
//------------------------------------------------------------------------------

/*============================== sc50a LP problem ==============================*
// Number of equations: 20
// Subspace dimension: 28
#define PP_PROBLEM_NAME		"sc50a"
#define PP_M 49	// Number of constraints in mps-file
#define PP_N 48		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 64.575077058564509026860413914575	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Zero point is feasible
//------------------------------------------------------------------------------

/*============================== sc50b LP problem ==============================*
#define PP_PROBLEM_NAME		"sc50b"
#define PP_M 48	// Number of constraints in mps-file
#define PP_N 48		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 70	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Zero point is feasible
//------------------------------------------------------------------------------

/*============================== scagr7 LP problem ==============================*
// Number of equations : 84
// Subspace dimension : 56
#define PP_PROBLEM_NAME	"scagr7"
#define PP_M 129		// Number of constraints in mps-file
#define PP_N 140		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 2331389.824330984	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 10
// Number of iterations: 253434
// Objective function: 1632971.02128593600355089
// Distance to polytope: 9.97917797778629684133842e-13
// Number of including inequality hyperplanes: 29
//------------------------------------------------------------------------------

/*============================== scsd1 LP problem ==============================*
// Number of equations : 77
// Subspace dimension : 683
#define PP_PROBLEM_NAME	"scsd1"
#define PP_M 77			// Number of constraints in mps-file
#define PP_N 760		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE -8.6666666743333647292533502995263	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 137
// Number of iterations: 170896
// Objective function: -14.7251531794086805149391
// Distance to polytope: 9.91047368421529525762281e-13
// Number of including inequality hyperplanes: 200
//------------------------------------------------------------------------------

/*============================== share2b LP problem ============================*
#define PP_PROBLEM_NAME		"share2b"	
#define PP_M 96	// Number of constraints in mps-file
#define PP_N 79	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 415.732240741419486545199108738 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 344
// Number of iterations: 20119607
// Objective function: 301.263653950768571121444
// Distance to polytope: 9.96934394263653700685132e-13
// Number of including inequality hyperplanes: 68
//------------------------------------------------------------------------------

/*============================== stocfor1 LP problem ============================*
// Number of equations: 63
// Subspace dimension: 48
#define PP_PROBLEM_NAME		"stocfor1"	
#define PP_M 117	// Number of constraints in mps-file
#define PP_N 111	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 41131.976219436406065682760731514 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION		(PP_EPS_ZERO*10)		// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE	(PP_EPS_PROJECTION*10)	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 1377
// Number of iterations: 42123668
// Objective function: 25426.8061500511248596013
// Distance to polytope: 9.96497898019737736204466e-13
// Number of including inequality hyperplanes: 33
//------------------------------------------------------------------------------

//============================================================================/**/