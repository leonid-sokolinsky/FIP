/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: BIP (Block-lterative Projection) method (No MPI)
Module: Problems24-1.h (LP problems of dimension 24 with 1 randome inequality)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
LP problems were obtained using BSF-LPP-Generator.
==============================================================================*/
#pragma once

//=========================== Method Parameters =========================
#define PP_ETA_TO_APEX				1000000			// Distance from apex base to apex point

/*============================== rnd24-0 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd24-0"
#define PP_M 25		// Number of equations (number of rows in *.mtx)
#define PP_N 49		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_INFINITY					1E+308			// Highest bound in *_hi.mtx
#define PP_ETA_TO_APEX				10000000		// Distance from apex base to apex point
//------------------------------------------------------------------/**/

/*============================== rnd24-1-1 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd24-1-1"
#define PP_M 25		// Number of equations (number of rows in *.mtx)
#define PP_N 49		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_INFINITY					1E+308			// Highest bound in *_hi.mtx
#define PP_ETA_TO_APEX				10000000		// Distance from apex base to apex point
//------------------------------------------------------------------/**/

/*============================== rnd24-1-2 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd24-1-2"
#define PP_M 25		// Number of equations (number of rows in *.mtx)
#define PP_N 49		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_INFINITY					1E+308			// Highest bound in *_hi.mtx
#define PP_ETA_TO_APEX				1000000000		// Distance from apex base to apex point
//------------------------------------------------------------------/**/

/*============================== rnd24-1-3 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd24-1-3"
#define PP_M 25		// Number of equations (number of rows in *.mtx)
#define PP_N 49		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_INFINITY					1E+308			// Highest bound in *_hi.mtx
#define PP_ETA_TO_APEX				10000000		// Distance from apex base to apex point
//------------------------------------------------------------------/**/

/*============================== rnd24-1-4 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd24-1-4"
#define PP_M 25		// Number of equations (number of rows in *.mtx)
#define PP_N 49		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_INFINITY					1E+308			// Highest bound in *_hi.mtx
#define PP_ETA_TO_APEX				10000000		// Distance from apex base to apex point
//------------------------------------------------------------------/**/

/*============================== rnd24-1-5 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd24-1-5"
#define PP_M 25		// Number of equations (number of rows in *.mtx)
#define PP_N 49		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_INFINITY					1E+308			// Highest bound in *_hi.mtx
#define PP_ETA_TO_APEX				10000000		// Distance from apex base to apex point
//------------------------------------------------------------------/**/