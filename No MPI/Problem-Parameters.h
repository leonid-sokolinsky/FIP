/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: FIP (Feasible point Iterative Projection) method (No MPI)
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
================================================================================*/

/*----------------------------- Problems ---------------------------------------*/

//#define PP_PATH "Problems/"

/**
#include "../_Problems05-1.h"
#define PP_PATH "../../Set-of-LP-Problems/Rnd-LP/"

/**
#include "../_Problems10-1.h"
#define PP_PATH "../../Set-of-LP-Problems/Rnd-LP/"

/**
#include "../_Problems-Miscellaneous.h"
#define PP_PATH "../../Set-of-LP-Problems/Miscellaneous-LP/"

/**/
#include "../_Problems-NetLib-LP.h"
#define PP_PATH "../../Set-of-LP-Problems/NetLib-LP/"

/*------------------------------------------------------------------------------*/

#define PP_METHOD_NAME "FIP (no MPI)"

//-------------------------- Compilation Modes ---------------------------------
#define PP_DEBUG
#define PP_MAXPROJECTION
//#define PP_MATRIX_OUTPUT
//#define PP_SAVE_ITER_IN_X0
#define PP_SAVE_RESULT

//================================ Problem Parameters ==========================
#define PP_MM				(PP_M+2*PP_N)		// Maximal number of constraints
#define PP_MAX_ITER_COUNT	10000000000		// Maximal count of iterations
//-------------------------- Input/Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	4 // Number of Elements to output
#define PP_SETW 24
//------------------------- Matrix format ----------------
#define PP_INFINITY			1E+308		// Highest bound in *_hi.mtx
#define PP_MTX_PREFIX		"lp_"
#define PP_MTX_POSTFIX_A	".mtx"
#define PP_MTX_POSTFIX_B	"_b.mtx"
#define PP_MTX_POSTFIX_C	"_c.mtx"
#define PP_MTX_POSTFIX_HI	"_hi.mtx"
#define PP_MTX_POSTFIX_LO	"_lo.mtx"
#define PP_MTX_POSTFIX_X0	"_x0.mtx" // Starting point
#define PP_MTX_POSTFIX_Z0	"_z0.mtx" // Boundary point
//------------------------- MPS format ----------------
#define PP_MPS_NAME_LENGTH	9
#define PP_MPS_PREFIX		"lp_"
#define PP_MPS_EXTENSION	".mps"
#define PP_MAX_NUMBER_OF_ROWS	(PP_M+3)
#define PP_MAX_NUMBER_OF_COLS	(PP_N*(PP_M+1)+1)
//-------------------------- Jobs  -----------------------
// Not used.
//------------- Locations -------------
#define PP_ON_HYPERPLANE			1
#define PP_INSIDE_HALFSPACE			2
#define PP_OUTSIDE_HALFSPACE		3