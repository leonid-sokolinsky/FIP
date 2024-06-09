/*==============================================================================
Project: LiFe
Theme: BIP (Block-lterative Projection) method (No MPI)
Module: Problems22-1.h (LP problems of dimension 22 with 1 randome inequality)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
LP problems were obtained using BSF-LPP-Generator.
=================================================================================*/
#pragma once

//=========================== Method Parameters =================================
#define PP_MAX_B_NO_CORRECT				200		// Maximum b that does not require correction
#define PP_EPS_POINT_IN_RND_HALFSPACE	1E-2	// Precision for random inequality in PointInPolytope()

/*============================== rnd22-1-1 LP problem ==============================*/
#define PP_PROBLEM_NAME	"rnd22-1-1"
#define PP_M 23		// Number of equations (number of rows in *.mtx)
#define PP_N 45		// Number of variables (number of cols in *.mtx)
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_COS					1E-8	// Precision for cos == 1
#define PP_EPS_TINY_PPROJ_VEC		1E-12	// Tiny pseudoprojecting vector
#define PP_EPS_POINT_INSIDE_CONE	1E-12	// Precision for point to be inside cone
#define PP_EPS_SHIFT				1E-9	// If shift length < PP_EPS_SHIFT then turn off cosine mode
//-------------------------------------------------------------------------------

//------------------------------------------------------------------/**/