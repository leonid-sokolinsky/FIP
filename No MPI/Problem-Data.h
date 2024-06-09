/*==============================================================================
Project: LiFe
Theme: BIP (Block-lterative Projection) method (No MPI)
Module: Problem-Data.h (Problem Data)
Prefix: PD
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Types.h"		// Problem Parameters 
using namespace std;
//========================== Problem variables ====================================
static int		PD_m;				// Current number of inequalities
static int		PD_n;				// Space dimension
static int		PD_iterNo;			// Number of iterations
static double	PD_cos;				// Cosine between previous and next projecting vectors
static double	PD_normPrevious_r;	// Norm of previous projecting vector
static double	PD_shiftLength;		// Shift length
static bool		PD_noFile_x0;		// Input file *_x0.mtx does not exist
//========================== Problem structures ====================================
static PT_matrix_T PD_A;			// Matrix of coefficients of inequalities
static PT_column_T PD_b;			// Column of the constant terms of the system Ax <= PD_b
static PT_vector_T PD_c;			// Objective Function Coefficients
static PT_vector_T PD_x;			// Current approximation
static PT_vector_T PD_hi;			// Higher bound
static PT_vector_T PD_lo;			// Lower bound
static PT_column_T PD_norm_a;		// Column of norms of matrix rows
static PT_vector_T PD_previous_r;	// Previous projecting vector
#ifdef PP_SCALE_RATIO
static PT_vector_T PD_minus_e;		// minus_e = (-1,...,-1)
#endif // PP_SCALE_RATIO
//========================== Input/Output ====================================
static string PD_problemName;