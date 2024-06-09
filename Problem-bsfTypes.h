/*==============================================================================
Project: LiFe
Theme: BIP (Block-lterative Projection) method (MPI)
Module: Problem-bsfTypes.h (Predefined BSF Problem Types)
Prefix: PT_bsf
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton
==============================================================================*/
#pragma once
#include "Problem-Types.h"		// Problem Types 
//=========================== BSF Types =========================
struct PT_bsf_parameter_T {
	PT_vector_T x;			// Current approximation
};

struct PT_bsf_mapElem_T {		// Type of map-list elements
	double* a;				// Pointer to row: a[0],...,a[n-1]
	double* b;				// Pointer to constant term b: a_0x_0+...+a_{n-1}x_{n-1} \leq b
};

struct PT_bsf_reduceElem_T {	// Type of reduce-list elements for Job 0 (default)	
	PT_vector_T projectingVector;	// Orthogonal projecting vector
	int nonZeroCounter;				// Counter of nonzero elements
};

struct PT_bsf_reduceElem_T_1 {	// Type of reduce-list elements for Job 1
	// Not used
};

struct PT_bsf_reduceElem_T_2 {	// Type of reduce-list elements for Job 2
	// Not used
};

struct PT_bsf_reduceElem_T_3 {	// Type of reduce-list elements for Job 3
// Not used
};