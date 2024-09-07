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

//=========================== Method Parameters =========================
#define PP_ETA_TO_APEX				1000000			// Distance from apex base to apex point

/*============================== rnd20-0 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd20-0"
#define PP_M 21		// Number of equations (number of rows in *.mtx)
#define PP_N 41		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_INFINITY					1E+308			// Highest bound in *_hi.mtx
#define PP_ETA_TO_APEX				100000			// Distance from apex base to apex point
//------------------------------------------------------------------/**/

/*============================== rnd20-1-1 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd20-1-1"
#define PP_M 21		// Number of equations (number of rows in *.mtx)
#define PP_N 41		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_INFINITY					1E+308			// Highest bound in *_hi.mtx
#define PP_ETA_TO_APEX				1000000			// Distance from apex base to apex point
//------------------------------------------------------------------/**/

/*============================== rnd20-1-2 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd20-1-2"
#define PP_M 21		// Number of equations (number of rows in *.mtx)
#define PP_N 41		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_INFINITY					1E+308			// Highest bound in *_hi.mtx
#define PP_ETA_TO_APEX				1000000			// Distance from apex base to apex point
//------------------------------------------------------------------/**/

/*============================== rnd20-1-3 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd20-1-3"
#define PP_M 21		// Number of equations (number of rows in *.mtx)
#define PP_N 41		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_INFINITY					1E+308			// Highest bound in *_hi.mtx
#define PP_ETA_TO_APEX				1000000			// Distance from apex base to apex point
//------------------------------------------------------------------/**/

/*============================== rnd20-1-4 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd20-1-4"
#define PP_M 21		// Number of equations (number of rows in *.mtx)
#define PP_N 41		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_INFINITY					1E+308			// Highest bound in *_hi.mtx
#define PP_ETA_TO_APEX				1000000			// Distance from apex base to apex point
//------------------------------------------------------------------/**/

/*============================== rnd20-1-5 LP problem ==============================*
#define PP_PROBLEM_NAME	"rnd20-1-5"
#define PP_M 21		// Number of equations (number of rows in *.mtx)
#define PP_N 41		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_INFINITY					1E+308			// Highest bound in *_hi.mtx
#define PP_ETA_TO_APEX				1000000			// Distance from apex base to apex point
//------------------------------------------------------------------/**/