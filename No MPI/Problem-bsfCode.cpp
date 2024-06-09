/*==============================================================================
Project: LiFe
Theme:	BIP (Block-lterative Projection) method (No MPI)
Module: Problem-bsfCode.cpp (Implementation of Problem Code)
Prefix:	PC_bsf	- BSF Predefined Problem Functions
		CF		- Common Functionc
		PF		- Private Functions
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	for (int j = 0; j < PD_n; j++)
		parameterOutP->x[j] = parameterIn.x[j];
}

void PC_bsf_Init(bool* success) {
	PD_problemName = PP_PROBLEM_NAME;

	*success = SF_MTX_Load__Problem();
	if (*success == false)
		return;

	*success = SF_MTX_LoadPoint(PD_x, PP_MTX_POSTFIX_X0);
	if (*success == false)
		return;

	SF_MakeColumnOfNorms(PD_A, PD_norm_a);

	if (SF_PointBelongsPolytope(PD_x, PP_EPS_ZERO)) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster) {
			if (PD_noFile_x0)
				cout << "Zero point is feasible.";
			else
				cout << "Starting point is feasible.";
			SF_MTX_SavePoint(PD_x, PP_MTX_POSTFIX_Z0);
			cout << " It is saved in file *_z0.mtx\n";
		}
		*success = false;
		return;
	}

#ifdef PP_SCALE_RATIO
	SF_Vector_MakeMinus_e(PD_minus_e);
	SF_PolytopeHomothety(PD_minus_e, PP_SCALE_RATIO);
	SF_PointHomothety(PD_x, PD_minus_e, PP_SCALE_RATIO);
#endif // PP_SCALE_RATIO

	PD_iterNo = 0;
#ifdef PP_COS_MODE
	PD_cos = 0;
	PD_normPrevious_r = 0;
#endif
}

void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {

	cout << "# " << BSF_sv_iterCounter << "\tTime " << round(elapsedTime);
	cout << "\tx =";
	SF_Print_Vector(parameter.x);
	cout << "\tF(t) = " << setw(PP_SETW) << SF_ObjF(parameter.x) << "\tDistance to M: " << setprecision(PP_BSF_PRECISION + 14) << SF_Distance_PointToPolytope(parameter.x) << endl;
	cout << setprecision(PP_BSF_PRECISION);
}

void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
};

void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
}

void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
}

void PC_bsf_JobDispatcher(PT_bsf_parameter_T* parameter, int* job, bool* exit, double t) {
	// Not used
}

void PC_bsf_MainArguments(int argc, char* argv[]) {
	// Not used
}

void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success) {
	int exitCode;
	PT_vector_T r; // Orthogonal projecting vector

	/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
	cout << PF_MAP_LIST_INDEX << ") ";
#endif // PP_DEBUG /**/

	SF_OrthogonalProjectingVectorOntoHalfspace_i(BSF_sv_parameter.x, PF_MAP_LIST_INDEX, r, PP_EPS_ON_HYPERPLANE, &exitCode);

	SF_Vector_Copy(r, reduceElem->projectingVector);

	/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
	cout << "||r|| = " << SF_Vector_Norm(r) << "\t\t";
#endif // PP_DEBUG /**/

	switch (exitCode) {
	case PP_EXITCODE_DEGENERATE_INEQUALITY:
		reduceElem->nonZeroCounter = 0;
		/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "DEGENERATE INEQUALITY:\treduceElem->nonZeroCounter = " << reduceElem->nonZeroCounter << endl;
#endif // PP_DEBUG /**/
		break;
	case PP_EXITCODE_INSIDE_HALFSPACE:
		reduceElem->nonZeroCounter = 0;
		/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "INSIDE HALFSPACE:\t\treduceElem->nonZeroCounter = " << reduceElem->nonZeroCounter << endl;
#endif // PP_DEBUG /**/
		break;
	case PP_EXITCODE_ON_HYPERPLANE:
		reduceElem->nonZeroCounter = 0;
		/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "ON HYPERPLANE:\t\treduceElem->nonZeroCounter = " << reduceElem->nonZeroCounter << endl;
#endif // PP_DEBUG /**/
		break;
	case PP_EXITCODE_NONDEGENERATE_PROJECTING:
		reduceElem->nonZeroCounter = 1;
		/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "NONDEGENERATE PROJECTION:\treduceElem->nonZeroCounter = " << reduceElem->nonZeroCounter << endl;
#endif // PP_DEBUG /**/
		break;
	default:
		cout << "PC_bsf_MapF: Unexpected exit code from function SF_OrthogonalProjectingVectorOntoHalfspace_i()!" << endl;
		*success = false;
		return;
	}
}

void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem, int* success) {
	// Not used
}

void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem, int* success) {
	// Not used
}

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {
	// Not used
}

void PC_bsf_MapInit(PT_bsf_parameter_T parameter) {
	// Not used
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "=================================================== BIP ====================================================" << endl;
	cout << "Problem name: " << PD_problemName << endl;
#ifdef PP_USE_LEASTPROJECTION
	cout << "Mode: Least projections " << endl;
#else
	cout << "Mode: Pseudoprojections " << endl;
#endif // PP_USE_LEASTPROJECTION
	cout << "No MPI" << endl;
#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
#else
	cout << "Number of Threads: " << omp_get_num_procs() << endl;
#endif // PP_BSF_NUM_THREADS
#else
	cout << "OpenMP is turned off!" << endl;
#endif // PP_BSF_OMP

#ifdef PP_BSF_FRAGMENTED_MAP_LIST
	cout << "Map List is Fragmented." << endl;
#else
	cout << "Map List is not Fragmented." << endl;
#endif
	cout << "Before conversion: m =\t" << PP_M << "\tn = " << PP_N << endl;
	cout << "After conversion:  m =\t" << PD_m << "\tn = " << PD_n << endl;
#ifdef PP_COS_MODE
	cout << "Cos mode is on." << endl;
#else
	cout << "Cos mode is off." << endl;
#endif
#ifdef PP_SCALE_RATIO
	cout << "Affine scaling: " << PP_SCALE_RATIO << endl;
#endif
	cout << "PP_EPS_ZERO = " << PP_EPS_ZERO << endl;
	cout << "PP_EPS_ZERO = " << PP_EPS_ZERO << endl;
	cout << "PP_EPS_ON_HYPERPLANE = " << PP_EPS_ON_HYPERPLANE << endl;
	cout << "PP_EPS_COS = " << PP_EPS_COS << endl;
	cout << "PP_EPS_MOVING = " << PP_EPS_MOVING << endl;
	cout << "PP_EPS_POINT_INSIDE_CONE = " << PP_EPS_POINT_INSIDE_CONE << endl;
	cout << "--------------- Data ---------------\n";
#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix PD_A & Column PD_b -------" << endl;
	SF_Print_Inequalities();
#endif // PP_MATRIX_OUTPUT
	cout << "Obj Function:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
		cout << setw(PP_SETW) << PD_c[j];
	if (PP_OUTPUT_LIMIT < PD_n)
		cout << " ...";
	cout << endl;
	cout << "x0 =\t\t";
	SF_Print_Vector(PD_x); cout << "\tF(x) = " << setw(PP_SETW) << SF_ObjF(PD_x) << endl;

#ifdef PP_DEBUG
	cout << "x0 on hyperplanes: "; SF_Print_Numbers_of_HyperplanesIncludingPoint(PD_x, PP_EPS_ZERO);
#endif // PP_DEBUG

	cout << "-------------------------------------------" << endl;
}

void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
#ifdef PP_SCALE_RATIO
	SF_PolytopeHomothety(PD_minus_e, 1 / (double)PP_SCALE_RATIO);
	SF_PointHomothety(parameter.x, PD_minus_e, 1 / (double)PP_SCALE_RATIO);
#endif // PP_SCALE_RATIO
	cout << "=============================================" << endl;
	cout << "Elapsed time: " << t << endl;
	cout << "Number of iterations: " << PD_iterNo << endl;
	cout << "=============================================" << endl;
	cout << "Boundary point z:\t";
	SF_Print_Vector(parameter.x); cout << endl;
	if (SF_PointBelongsPolytope(parameter.x, PP_EPS_ON_HYPERPLANE))
		cout << "z on polytope.\n";
	else
		cout << "z NOT in polytope!!!\n";
	cout << "Distance to polytope: " << setprecision(24) << SF_Distance_PointToPolytope(parameter.x) << endl;
	cout << "z on hyperplanes: "; SF_Print_Numbers_of_HyperplanesIncludingPoint(parameter.x, PP_EPS_ON_HYPERPLANE);

#ifdef PP_SAVE_RESULT
	if (SF_MTX_SavePoint(parameter.x, PP_MTX_POSTFIX_Z0))
		cout << "\nBoundary point is saved into the file *.z0" << endl;
#endif // PP_SAVE_RESULT
}

void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// Not used
}

void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// Not used
}

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// Not used
}

void PC_bsf_ProcessResults(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* exit) {
	PD_iterNo++;

#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "-------------> PC_bsf_ProcessResults: Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = "
			<< PP_MAX_ITER_COUNT << endl;
		*exit = true;
		return;
	};
#endif // PP_MAX_ITER_COUNT

	if (reduceResult->nonZeroCounter == 0) {
		*exit = true;
		return;
	}

	SF_Vector_DivideEquals(reduceResult->projectingVector, (double)(reduceResult->nonZeroCounter));

	/*DEBUG PC_bsf_ProcessResults**
#ifdef PP_DEBUG
	cout << "_______________________________ " << PD_iterNo << " ___________________________________\n";
	cout << "||r|| = " << SF_Vector_Norm(reduceResult->projectingVector) << endl;
#endif // PP_DEBUG /**/

	if (SF_Vector_Is_Tiny(reduceResult->projectingVector, PP_EPS_ZERO)) {
		*exit = true;
		return;
	}

	SF_Vector_Addition(parameter->x, reduceResult->projectingVector, parameter->x);

	double norm_r = SF_Vector_Norm(reduceResult->projectingVector);

	if (PD_normPrevious_r != 0) {
		PD_cos = SF_Vector_DotProduct(PD_previous_r, reduceResult->projectingVector) / (PD_normPrevious_r * norm_r);
	}

	SF_Vector_Copy(reduceResult->projectingVector, PD_previous_r);
	PD_normPrevious_r = norm_r;

	/*DEBUG PC_bsf_ProcessResults**
#ifdef PP_DEBUG
	//cout << "x = "; SF_Print_Vector(parameter->x); cout << endl;
	cout << "x on hyperplanes: "; SF_Print_Numbers_of_HyperplanesIncludingPoint(parameter->x, PP_EPS_ZERO);
	cout << "cos = " << PD_cos << endl;
	cout << "Distance to polytope: " << SF_Distance_PointToPolytope(parameter->x) << endl;
#endif // PP_DEBUG /**/

#ifdef PP_COS_MODE
	if (1 - PD_cos >= PP_EPS_COS)
		return;

	SF_MovingTowardsPolytope(parameter->x, reduceResult->projectingVector, PP_EPS_MOVING);

	/*DEBUG PC_bsf_ProcessResults**
#ifdef PP_DEBUG
	cout << "--------------------------------- SF_MovingTowardsPolytope -----------------------------------------\n";
	cout << "SF_Shift = " << PD_shiftLength << endl;
	//cout << "x = "; SF_Print_Vector(parameter->x); cout << endl;
	cout << "x on hyperplanes: "; SF_Print_Numbers_of_HyperplanesIncludingPoint(parameter->x, PP_EPS_ZERO);
	cout << "Distance to polytope: " << SF_Distance_PointToPolytope(parameter->x) << endl;
#endif // PP_DEBUG /**/
#endif
}

void PC_bsf_ProcessResults_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* exit) {
	// Not used
}

void PC_bsf_ProcessResults_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* exit) {
	// Not used
}

void PC_bsf_ProcessResults_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* exit) {
	// Not used
}

void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	SF_Vector_Addition(x->projectingVector, y->projectingVector, z->projectingVector);
	z->nonZeroCounter = x->nonZeroCounter + y->nonZeroCounter;
}

void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {
	// Not used
}

void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {
	// Not used
}

void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {
	// Not used
}

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	SF_Vector_Copy(PD_x, parameter->x);
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PD_m;
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->a = PD_A[i];
	elem->b = &(PD_b[i]);
}

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; };
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; };
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; };
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; };
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; };
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; };
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; };
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; };

//---------------------------------- Shared Functions -------------------------

static void SF_Conversion() { // Transformation to inequalities & dimensionality reduction
	int m_equation = PD_m;
	int m_inequality;
	int m_lowerBound;
	int m_higherBound;

	for (int i = 0; i < m_equation; i++) { // Conversion to inequalities
		for (int j = 0; j < PD_n; j++)
			PD_A[PD_m][j] = -PD_A[i][j];
		PD_b[PD_m] = -PD_b[i];
		PD_m++; assert(PD_m <= PP_MM);
	}

	for (int i = 0; i < PD_m; i++) // Remove negative sign for zero value
		for (int j = 0; j < PD_n; j++)
			if (PD_A[i][j] == 0)
				PD_A[i][j] = 0;
	m_inequality = PD_m;

	for (int i = 0; i < PD_n; i++) { // Adding lower bound conditions
		for (int j = 0; j < PD_n; j++)
			PD_A[i + PD_m][j] = 0;
		PD_A[i + PD_m][i] = -1;
		if (PD_lo[i] == 0)
			PD_b[i + PD_m] = 0;
		else
			PD_b[i + PD_m] = -PD_lo[i];
	}
	PD_m += PD_n; assert(PD_m <= PP_MM);
	m_lowerBound = PD_m;

	for (int i = 0; i < PD_n; i++) { // Adding higher bound conditions
		if (PD_hi[i] != PP_INFINITY) {
			for (int j = 0; j < PD_n; j++)
				PD_A[PD_m][j] = 0;
			PD_A[PD_m][i] = 1;
			PD_b[PD_m] = PD_hi[i];
			PD_m++; assert(PD_m <= PP_MM);
		}
	}
	m_higherBound = PD_m;

	/**
	cout << "-----------------------------------------------------\n";
	SF_Print_Inequalities();
	cout << "-----------------------------------------------------\n";
	cout << "PD_c: "; SF_Print_Vector(PD_c); cout << endl;/**/

	SF_RemoveFreeVariables(m_equation, m_inequality, m_lowerBound, m_higherBound);

	/**
	cout << "-----------------------------------------------------\n";
	SF_Print_Inequalities();
	cout << "-----------------------------------------------------\n";/**/
}

inline void SF_ConversionSimple() { // Transformation to inequalities
	int m_equation = PD_m;

	for (int i = 0; i < m_equation; i++) { // Conversion to inequalities
		for (int j = 0; j < PD_n; j++)
			PD_A[PD_m][j] = -PD_A[i][j];
		PD_b[PD_m] = -PD_b[i];
		PD_m++; assert(PD_m <= PP_MM);
	}

	for (int i = 0; i < PD_m; i++) // Remove negative sign for zero value
		for (int j = 0; j < PD_n; j++)
			if (PD_A[i][j] == 0)
				PD_A[i][j] = 0;

	for (int i = 0; i < PD_n; i++) { // Adding lower bound conditions
		for (int j = 0; j < PD_n; j++)
			PD_A[i + PD_m][j] = 0;
		PD_A[i + PD_m][i] = -1;
		if (PD_lo[i] == 0)
			PD_b[i + PD_m] = 0;
		else
			PD_b[i + PD_m] = -PD_lo[i];
	}
	PD_m += PD_n; assert(PD_m <= PP_MM);

	for (int i = 0; i < PD_n; i++) { // Adding higher bound conditions
		if (PD_hi[i] != PP_INFINITY) {
			for (int j = 0; j < PD_n; j++)
				PD_A[PD_m][j] = 0;
			PD_A[PD_m][i] = 1;
			PD_b[PD_m] = PD_hi[i];
			PD_m++; assert(PD_m <= PP_MM);
		}
	}
}

inline double SF_Distance_PointToHalfspace_i(PT_vector_T x, int i) {
	double a_dot_z_minus_b = SF_Vector_DotProduct(PD_A[i], x) - PD_b[i];

	if (PD_norm_a[i] < PP_EPS_ZERO) //Degenerate equation
		return 0;

	if (a_dot_z_minus_b < 0) // Point belongs to halfspace
		return 0;

	return a_dot_z_minus_b / PD_norm_a[i];
}

inline double SF_Distance_PointToHyperplane_i(PT_vector_T x, int i) {
	if (PD_norm_a[i] < PP_EPS_ZERO) //Degenerate equation
		return 0;
	else
		return fabs(SF_Vector_DotProduct(PD_A[i], x) - PD_b[i]) / PD_norm_a[i];
}

inline double SF_Distance_PointToPoint(PT_vector_T x, PT_vector_T y) {
	PT_vector_T z;
	SF_Vector_Subtraction(x, y, z);
	return SF_Vector_Norm(z);
}

inline double SF_Distance_PointToPolytope(PT_vector_T x) { // Measure of distance from point to polytope
	double maxDistance = 0;
	double distance;

	for (int i = 0; i < PD_m; i++) {
		distance = SF_Distance_PointToHalfspace_i(x, i);
		if (distance > 0)
			maxDistance = PF_MAX(maxDistance, distance);
	}
	return maxDistance;
}

inline void SF_MakeColumnOfNorms(PT_matrix_T A, PT_column_T norm_a) {
	for (int i = 0; i < PD_m; i++)
		norm_a[i] = SF_Vector_Norm(A[i]);
}

inline void SF_MovingOnPolytope(PT_vector_T startPoint, PT_vector_T directionVector, PT_vector_T finishPoint, double eps) {
	double leftBound = 0;
	double rightBound = PP_DBL_MAX;
	double factor;
	PT_vector_T unitDirectionVector;

	double vectorLength = SF_Vector_Norm(directionVector);
	assert(vectorLength >= eps);

	SF_Vector_DivideByNumber(directionVector, vectorLength, unitDirectionVector);


	PD_shiftLength = PP_START_SHIFT_LENGTH;
	factor = PD_shiftLength;

	while (rightBound - leftBound >= eps) {
		SF_Shift(startPoint, unitDirectionVector, PD_shiftLength, finishPoint);
		if (SF_PointBelongsPolytope(finishPoint, eps / 10)) {
			leftBound = PD_shiftLength;
			PD_shiftLength += factor;
		}
		else {
			rightBound = PD_shiftLength;
			factor /= 2;
			PD_shiftLength -= factor;
		}
	}

	PD_shiftLength -= factor;
	SF_Shift(startPoint, unitDirectionVector, PD_shiftLength, finishPoint);
	return;
}

inline void SF_MovingTowardsPolytope(PT_vector_T point, PT_vector_T directionVector, double eps) {
	double leftBound = 0;
	double rightBound = PP_DBL_MAX;
	double factor;
	double vectorLength;
	PT_vector_T shiftedPoint;
	PT_vector_T unitDirectionVector;
	static int outerHalspace_i[PP_MM];	// Index of out half-spaces
	int mo;								// Number of out half-spaces
	bool pointInsideCone;

	mo = 0;
	for (int i = 0; i < PD_m; i++)
		if (!SF_PointBelongsHalfspace_i(point, i, eps)) {
			outerHalspace_i[mo] = i;
			mo++;
		}

	vectorLength = SF_Vector_Norm(directionVector);
	if (vectorLength < eps)
		return;

	SF_Vector_DivideByNumber(directionVector, SF_Vector_Norm(directionVector), unitDirectionVector);

	PD_shiftLength = PP_START_SHIFT_LENGTH;
	factor = PD_shiftLength;

	while (rightBound - leftBound >= eps) {
		SF_Shift(point, unitDirectionVector, PD_shiftLength, shiftedPoint);

		pointInsideCone = true;
		for (int i = 0; i < mo; i++)
			if (SF_PointBelongsHalfspace_i(shiftedPoint, outerHalspace_i[i], eps)) {
				pointInsideCone = false;
				break;
			}

		if (pointInsideCone) {
			leftBound = PD_shiftLength;
			PD_shiftLength += factor;
		}
		else {
			rightBound = PD_shiftLength;
			factor /= 2;
			PD_shiftLength -= factor;
		}
	}

	PD_shiftLength -= factor;
	SF_Shift(point, unitDirectionVector, PD_shiftLength, point);
}

static bool SF_MTX_Load__Problem() {

	//--------------- Reading A ------------------
	if (!SF_MTX_Load_A())
		return false;

	//--------------- Reading b ------------------
	if (!SF_MTX_Load_b())
		return false;

	//--------------- Reading lo ------------------
	if (!SF_MTX_Load_lo())
		return false;

	//--------------- Reading c ------------------
	if (!SF_MTX_Load_c())
		return false;

	//--------------- Reading hi ------------------
	if (!SF_MTX_Load_hi())
		return false;

	//---------- Conversion to inequalities -----------
#ifdef PP_SIMPLE_CONVERSION
	SF_ConversionSimple();
#else
	SF_Conversion();
#endif

	/**
	cout << "-----------------------------------------------------\n";
	SF_Print_Inequalities();
	cout << "-----------------------------------------------------\n";
	cout << "PD_c: "; SF_Print_Vector(PD_c); cout << endl;/**/

	return true;
}

inline bool SF_MTX_Load_A() {
	int nor;	// Number of matrix rows
	int noc;	// Number of matrix columns
	int non;	// Number of non-zero elements
	const char* mtxFile;
	FILE* stream;// Input stream
	char str[80] = { '\0' };
	char* chr = str;
	string PD_MTX_File;

	PD_MTX_File = PP_PATH;
	PD_MTX_File += PP_MTX_PREFIX;
	PD_MTX_File += PD_problemName;
	PD_MTX_File += PP_MTX_POSTFIX_A;
	mtxFile = PD_MTX_File.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SF_MTX_SkipComments(stream);
	if (fscanf(stream, "%d%d%d", &nor, &noc, &non) < 3) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Unexpected end of file " << mtxFile << endl;
		return false;
	}

	if (nor >= noc) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Number of rows m = " << nor << " must be < " << "Number of columns n = " << noc << "\n";
		return false;
	}

	if (noc != PP_N) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Invalid input data: PP_N must be = " << noc << "\n";
		return false;
	}

	if (nor != PP_M) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Invalid input data: PP_M must be = " << nor << "\n";
		return false;
	}

	PD_m = nor;
	PD_n = noc;

	if (2 * nor + noc > PP_MM) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Invalid input data: number of inequalities m = " << 2 * nor + noc
			<< " must be < PP_MM + 1 =" << PP_MM + 1 << "\n";
		return false;
	}

	for (int k = 0; k < non; k++) {
		int i, j;

		if (fscanf(stream, "%d%d%s", &i, &j, str) < 3) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}

		i -= 1;
		j -= 1;
		if (i < 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Negative row index in'" << mtxFile << "'.\n" << endl;
			return false;
		}
		if (j < 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Negative column index in'" << mtxFile << "'.\n" << endl;
			return false;
		}
		PD_A[i][j] = strtod(str, &chr);
	}

	fclose(stream);

	return true;
}

inline bool SF_MTX_Load_b() {
	int nor;	// Number of matrix rows
	int noc;	// Number of matrix columns
	const char* mtxFile;
	FILE* stream;// Input stream
	char str[80] = { '\0' };
	char* chr = str;
	string PD_MTX_File;

	PD_MTX_File = PP_PATH;
	PD_MTX_File += PP_MTX_PREFIX;
	PD_MTX_File += PD_problemName;
	PD_MTX_File += PP_MTX_POSTFIX_B;
	mtxFile = PD_MTX_File.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SF_MTX_SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (PD_m != nor) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int i = 0; i < PD_m; i++) {
		if (fscanf(stream, "%s", str) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Unexpected end of file '" << mtxFile << "'." << endl;
			return false;
		}
		PD_b[i] = strtod(str, &chr);
	}
	fclose(stream);

	return true;
}

inline bool SF_MTX_Load_c() {
	int nor;	// Number of matrix rows
	int noc;	// Number of matrix columns
	const char* mtxFile;
	FILE* stream;// Input stream
	char str[80] = { '\0' };
	char* chr = str;
	string PD_MTX_File;

	PD_MTX_File = PP_PATH;
	PD_MTX_File += PP_MTX_PREFIX;
	PD_MTX_File += PD_problemName;
	PD_MTX_File += PP_MTX_POSTFIX_C;
	mtxFile = PD_MTX_File.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SF_MTX_SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%s", str) < 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Unexpected end of file" << endl;
			return false;
		}
		PD_c[j] = -strtod(str, &chr);
	}
	fclose(stream);

	return true;
}

inline bool SF_MTX_Load_hi() {
	int nor;	// Number of matrix rows
	int noc;	// Number of matrix columns
	const char* mtxFile;
	FILE* stream;// Input stream
	char str[80] = { '\0' };
	char* chr = str;
	string PD_MTX_File;

	PD_MTX_File = PP_PATH;
	PD_MTX_File += PP_MTX_PREFIX;
	PD_MTX_File += PD_problemName;
	PD_MTX_File += PP_MTX_POSTFIX_HI;
	mtxFile = PD_MTX_File.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SF_MTX_SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%s", str) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file '" << mtxFile << "'." << endl;
			return false;
		}
		PD_hi[j] = strtod(str, &chr);
	}
	fclose(stream);
	return true;
}

inline bool SF_MTX_Load_lo() {
	int nor;	// Number of matrix rows
	int noc;	// Number of matrix columns
	const char* mtxFile;
	FILE* stream;// Input stream
	char str[80] = { '\0' };
	char* chr = str;
	string PD_MTX_File;

	PD_MTX_File = PP_PATH;
	PD_MTX_File += PP_MTX_PREFIX;
	PD_MTX_File += PD_problemName;
	PD_MTX_File += PP_MTX_POSTFIX_LO;
	mtxFile = PD_MTX_File.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SF_MTX_SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%s", str) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Unexpected end of file '" << mtxFile << "'." << endl;
			return false;
		}
		PD_lo[j] = strtod(str, &chr);
	}

	fclose(stream);

	return true;
}

inline bool SF_MTX_LoadPoint(PT_vector_T x, string postfix) {
	int nor;	// Number of matrix rows
	int noc;	// Number of matrix columns
	const char* mtxFile;
	FILE* stream;// Input stream
	char str[80] = { '\0' };
	char* chr = str;
	string PD_MTX_File;

	PD_MTX_File = PP_PATH;
	PD_MTX_File += PP_MTX_PREFIX;
	PD_MTX_File += PD_problemName;
	PD_MTX_File += postfix;
	mtxFile = PD_MTX_File.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SF_MTX_SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%s", str) < 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Unexpected end of file" << endl;
			return false;
		}
		x[j] = strtod(str, &chr);
	}
	fclose(stream);

	return true;
}

static bool SF_MTX_SavePoint(PT_vector_T x, string postfix) {
	const char* mtxFile;
	FILE* stream;// Input stream
	string PD_MTX_File;

	PD_MTX_File = PP_PATH;
	PD_MTX_File += PP_MTX_PREFIX;
	PD_MTX_File += PD_problemName;
	PD_MTX_File += postfix;
	mtxFile = PD_MTX_File.c_str();
	stream = fopen(mtxFile, "w");
	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	fprintf(stream, "%d %d\n", PD_n, 1);

	for (int j = 0; j < PD_n; j++)
		fprintf(stream, "%.16f\n", x[j]);

	fclose(stream);
	return true;
}

inline void SF_MTX_SkipComments(FILE* stream) {
	fpos_t pos;	// Position in the input stream
	int res;
	res = fscanf(stream, "\n");
	fgetpos(stream, &pos);
	while (getc(stream) == '%') {
		while (getc(stream) != 10);
		res = fscanf(stream, "\n");
		fgetpos(stream, &pos);
	};
	fsetpos(stream, &pos);
}

inline double SF_ObjF(PT_vector_T x) {
	double s = 0;
	for (int j = 0; j < PD_n; j++)
		s += PD_c[j] * x[j];
	return s;
}

inline void	SF_ObliqueProjectingVectorOntoHalfspace_i(PT_vector_T z, int i, PT_vector_T g, PT_vector_T o, double eps, int* exitCode) {
	// Oblique projecting vector o of point z onto Half-space H_i with respect to vector g
	double a_dot_g;	// <a,g>
	double a_dot_z_minus_b;	// <a,z> - b
	double factor;	// (b - <a,z>) / <a,g>

	if (PD_norm_a[i] < PP_EPS_ZERO) {
		SF_Vector_Zeroing(o);
		*exitCode = 0; // #define PP_EXITCODE_DEGENERATE_INEQUALITY		0
		return;
	}

	a_dot_z_minus_b = SF_Vector_DotProduct(PD_A[i], z) - PD_b[i]; // <a,z> - b

	if (fabs(a_dot_z_minus_b) < PP_EPS_ZERO) { // <a,z> - b < 0
		*exitCode = 1; // #define PP_EXITCODE_ON_HYPERPLANE				1
		SF_Vector_Zeroing(o);
		return;
	}

	if (a_dot_z_minus_b <= -PP_EPS_ZERO) { // <a,z> - b < 0
		*exitCode = 2; // #define PP_EXITCODE_INSIDE_HALFSPACE			2
		SF_Vector_Zeroing(o);
		return;
	}

	a_dot_g = SF_Vector_DotProduct(PD_A[i], g); // <a,g>


	if (fabs(a_dot_g) < PP_EPS_ZERO) {
		*exitCode = 4; // #define PP_EXITCODE_PARALLEL				4
		SF_Vector_Zeroing(o);
		return;
	}

	if (a_dot_g >= PP_EPS_ZERO) {
		*exitCode = 5; // #define PP_EXITCODE_RECESSIVE				5
		SF_Vector_Zeroing(o);
		return;
	}

	factor = a_dot_z_minus_b / a_dot_g; // (<a,z> - b) / <a,g>

	// Oblique projection vector: o = -(<a,z> - b)g/<a, g> = -factor * g
	SF_Vector_MultiplyByNumber(g, -factor, o);

	*exitCode = 9; // #define PP_EXITCODE_NONDEGENERATE_PROJECTING	9
	return;
}

inline void SF_OrthogonalProjectingVectorOntoHalfspace_i(PT_vector_T z, int i, PT_vector_T r, double eps, int* exitCode) {
	double factor;
	double a_dot_z_minus_b = SF_Vector_DotProduct(PD_A[i], z) - PD_b[i]; // <a,z> - b
	double distance = fabs(a_dot_z_minus_b) / PD_norm_a[i];

	if (PD_norm_a[i] < PP_EPS_ZERO) {
		SF_Vector_Zeroing(r);
		*exitCode = 0; // #define PP_EXITCODE_DEGENERATE_INEQUALITY		0
		return;
	}

	if (distance < eps) {
		SF_Vector_Zeroing(r);
		*exitCode = 1; // #define PP_EXITCODE_ON_HYPERPLANE				1
		return;
	}

	if (a_dot_z_minus_b < 0) { // <a,z> - b < 0
		SF_Vector_Zeroing(r);
		*exitCode = 2; // #define PP_EXITCODE_INSIDE_HALFSPACE			2
		return;
	}

	factor = -a_dot_z_minus_b / (PD_norm_a[i] * PD_norm_a[i]); // (b - <z,a>) / ||a||^2
	SF_Vector_MultiplyByNumber(PD_A[i], factor, r); // r = a(b - <z,a>) / ||a||^2
	*exitCode = 9; // #define PP_EXITCODE_NONDEGENERATE_PROJECTING	9
}

inline void SF_OrthogonalProjectingVectorOntoHyperplane_i(PT_vector_T x, int i, PT_vector_T p) {
	assert(SF_Vector_NormSquare(PD_A[i]));
	SF_Vector_MultiplyByNumber(PD_A[i], -(SF_Vector_DotProduct(PD_A[i], x) - PD_b[i]) / SF_Vector_NormSquare(PD_A[i]), p);
}

inline bool SF_PointBelongsHalfspace_i(PT_vector_T x, int i, double eps) {
	if (PD_norm_a[i] < eps) //Degenerate equation
		return true;
	double a_dot_x_minus_b = SF_Vector_DotProduct(PD_A[i], x) - PD_b[i];
	double distanceToHyperplane = fabs(a_dot_x_minus_b) / PD_norm_a[i];
	if (distanceToHyperplane < eps)
		return true;
	if (a_dot_x_minus_b < 0)
		return true;
	return false;
}

inline bool SF_PointBelongsHyperplane_i(PT_vector_T x, int i, double eps) {
	if (SF_Distance_PointToHyperplane_i(x, i) < eps)
		return true;
	else
		return false;
}

inline bool SF_PointBelongsPolytope(PT_vector_T x, double eps) { // If the point belongs to the polytope with prescigion of eps
	for (int i = 0; i < PD_m; i++)
		if (!SF_PointBelongsHalfspace_i(x, i, eps))
			return false;
	return true;
}

inline void SF_PointHomothety(PT_vector_T x, PT_vector_T center, double ratio) { // https://en.wikipedia.org/wiki/Homothety
	if (ratio == 1)
		return;
	assert(ratio > 0);
	for (int j = 0; j < PD_n; j++)
		x[j] = ratio * x[j] - (ratio - 1) * center[j];
}

inline bool SF_PointInsideHalfspace_i(PT_vector_T x, int i, double eps) {
	if (PD_norm_a[i] < PP_EPS_ZERO) //Degenerate equation
		return true;
	double a_dot_x_minus_b = SF_Vector_DotProduct(PD_A[i], x) - PD_b[i];
	double distanceToHyperplane = fabs(a_dot_x_minus_b) / PD_norm_a[i];
	if (distanceToHyperplane < eps)
		return false;
	if (a_dot_x_minus_b < 0)
		return true;
	return false;
}

inline void SF_PolytopeHomothety(PT_vector_T center, double ratio) { // https://en.wikipedia.org/wiki/Homothety
	if (ratio == 1)
		return;
	assert(ratio > 0);

	for (int i = 0; i < PD_m; i++) {
		PD_b[i] = ratio * PD_b[i] - (ratio - 1) * SF_Vector_DotProduct(PD_A[i], center);
	}
}

inline void SF_Print_Inequalities() {
	for (int i = 0; i < PD_m; i++) {
		cout << i << ")";
		for (int j = 0; j < PD_n; j++)
			cout << setw(PP_SETW) << PD_A[i][j];
		cout << "\t<=" << setw(PP_SETW) << PD_b[i] << endl;
	}
}

inline void SF_Print_Numbers_of_HyperplanesIncludingPoint(PT_vector_T x, double eps) {
	for (int i = 0; i < PD_m; i++)
		if (SF_PointBelongsHyperplane_i(x, i, eps))
			cout << i << " ";
	cout << endl;
}

inline void SF_Print_Vector(PT_vector_T x) {
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << x[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
}

inline double SF_RelativeError(double trueValue, double calculatedValue) {
	if (fabs(trueValue) >= PP_EPS_ZERO)
		return fabs(calculatedValue - trueValue) / fabs(trueValue);
	else
		return fabs(calculatedValue - trueValue);
}

inline void SF_RemoveFreeVariables(int m_equation, int m_inequality, int m_lowerBound, int m_higherBound) {
	int freeVariable_j_col[PP_M];
	int m_freeVariable_j_col = 0;
	int rowToDelete_i[PP_MM];
	int m_rowToDelete = 0;
	bool unique;
	bool ok;

	// Find free variables
	for (int j_col = 0; j_col < PD_n; j_col++) { // Find zero element in PD_c
		if (PD_c[j_col] != 0)
			continue;
		for (int i_row = 0; i_row < m_equation; i_row++) { // Find PD_A i_row with nonzero element in column j_col
			if (PD_A[i_row][j_col] == 0)
				continue;

			// Check uniqueness in column j_col
			unique = true;
			for (int i = 0; i < m_equation; i++) {
				if (i == i_row)
					continue;
				if (PD_A[i][j_col] != 0) {
					unique = false;
					break;
				}
			}
			if (!unique)
				continue;

			// Check lower bound
			ok = true;
			for (int i_lowerBound = m_inequality; i_lowerBound < m_lowerBound; i_lowerBound++) {
				if (PD_A[i_lowerBound][j_col] == 0)
					continue;
				if (PD_A[i_lowerBound][j_col] != -1) {
					ok = false;
					break;
				}
				if (PD_b[i_lowerBound] != 0) {
					ok = false;
					break;
				}
				rowToDelete_i[m_rowToDelete] = i_lowerBound;
				m_rowToDelete++; assert(m_rowToDelete <= PP_MM);
				break;
			}
			if (!ok)
				continue;

			// Check higher bound
			ok = true;
			for (int i_higherBound = m_lowerBound; i_higherBound < m_higherBound; i_higherBound++) {
				if (PD_A[i_higherBound][j_col] != 0) {
					ok = false;
					m_rowToDelete--;
					break;
				}
			}
			if (!ok)
				continue;

			freeVariable_j_col[m_freeVariable_j_col] = j_col;
			m_freeVariable_j_col++; assert(m_freeVariable_j_col <= PP_M);
			rowToDelete_i[m_rowToDelete] = m_equation + i_row;
			m_rowToDelete++; assert(m_rowToDelete <= PP_MM);
			// Check sign of free variable
			if (PD_A[i_row][j_col] < 0) {
				// Change sign of inequality
				for (int j = 0; j < PD_n; j++)
					PD_A[i_row][j] = -PD_A[i_row][j];
				PD_b[i_row] = -PD_b[i_row];
			}
			break;
		}
	}

	{// Eliminate columns with free variables
		static bool colToDeleteLable[PP_N];
		for (int j = 0; j < m_freeVariable_j_col; j++) {
			assert(freeVariable_j_col[j] < PP_N);
			colToDeleteLable[freeVariable_j_col[j]] = true;
		}

		for (int j = 0; j < PD_n; j++) {
			if (colToDeleteLable[j]) {
				for (int i = 0; i < PD_m; i++)
					PD_A[i][j] = PD_A[i][PD_n - 1];
				PD_c[j] = PD_c[PD_n - 1];
				colToDeleteLable[j] = colToDeleteLable[PD_n - 1];
				j--; assert(j >= 0);
				PD_n--; assert(PD_n >= 0);

				/**
				cout << "-----------------------------------------------------\n";
				SF_Print_Inequalities();
				cout << "-----------------------------------------------------\n";/**/

			}
		}
	}

	{// Eliminate rows corresponding to free variables
		static bool rowToDeleteLable[PP_MM];
		for (int i = 0; i < m_rowToDelete; i++) {
			rowToDeleteLable[rowToDelete_i[i]] = true;
		}
		for (int i = 0; i < PD_m; i++) {
			if (rowToDeleteLable[i]) {
				for (int j = 0; j < PD_n; j++)
					PD_A[i][j] = PD_A[PD_m - 1][j];
				PD_b[i] = PD_b[PD_m - 1];
				rowToDeleteLable[i] = rowToDeleteLable[PD_m - 1];
				i--; assert(i >= 0);
				PD_m--; assert(PD_m >= 0);

				/**
				cout << "-----------------------------------------------------\n";
				SF_Print_Inequalities();
				cout << "-----------------------------------------------------\n";/**/
			}
		}
	}
}

inline void SF_Shift(PT_vector_T point, PT_vector_T unitVector, double shiftLength, PT_vector_T shiftedPoint) {
	for (int j = 0; j < PD_n; j++)
		shiftedPoint[j] = point[j] + unitVector[j] * shiftLength;
}

inline void SF_Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x + y
	for (int j = 0; j < PD_n; j++)
		z[j] = x[j] + y[j];
}

inline bool SF_Vector_ApproximateEquality(PT_vector_T x, PT_vector_T y, double eps) { // x = y
	for (int j = 0; j < PD_n; j++)
		if (fabs(x[j] - y[j]) >= eps)
			return false;
	return true;
}

inline void SF_Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
	for (int j = 0; j < PD_n; j++)
		toPoint[j] = fromPoint[j];
}

inline void SF_Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = x/r
	for (int j = 0; j < PD_n; j++)
		y[j] = x[j] / r;
}

inline void SF_Vector_DivideEquals(PT_vector_T x, double r) {  // x = x/r
	for (int j = 0; j < PD_n; j++)
		x[j] /= r;
}

inline double SF_Vector_DotProduct(PT_vector_T x, PT_vector_T y) {
	double sum = 0;
	for (int j = 0; j < PD_n; j++)
		sum += x[j] * y[j];
	return sum;
}

inline bool SF_Vector_Is_Tiny(PT_vector_T x, double eps) {
	return SF_Vector_Norm(x) < eps;
}

inline void SF_Vector_MakeLike(PT_vector_T x, double lengthOfLikeVector, PT_vector_T likeVector) {
	double norm_x = SF_Vector_Norm(x);
	if (norm_x == 0)
		SF_Vector_Zeroing(likeVector);
	else
		SF_Vector_MultiplyByNumber(x, lengthOfLikeVector / norm_x, likeVector);
}

void SF_Vector_MakeMinus_e(PT_vector_T minus_e) {
	for (int j = 0; j < PD_n; j++)
		minus_e[j] = -1;
}

inline void SF_Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector) { // equalPoint += minusVector
	for (int j = 0; j < PD_n; j++)
		equalPoint[j] -= minusVector[j];
}

inline void SF_Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = r*x
	for (int j = 0; j < PD_n; j++)
		y[j] = x[j] * r;
}

inline void SF_Vector_MultiplyEquals(PT_vector_T x, double r) {  // x = r*x
	for (int j = 0; j < PD_n; j++)
		x[j] *= r;
}

inline double SF_Vector_Norm(PT_vector_T x) {
	return sqrt(SF_Vector_NormSquare(x));
}

inline double SF_Vector_NormSquare(PT_vector_T x) {
	double sum = 0;

	for (int j = 0; j < PD_n; j++) {
		sum += x[j] * x[j];
	}
	return sum;
}

inline void SF_Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector) { // equalVector += plusVector
	for (int j = 0; j < PD_n; j++)
		equalVector[j] += plusVector[j];
}

inline void SF_Vector_Round(PT_vector_T x, double eps) {
	double floorValue;
	double fractionalPart;
	double sign;
	double absValue;

	for (int j = 0; j < PD_n; j++) {
		if (fabs(x[j]) < eps) {
			x[j] = 0;
			continue;
		}
		absValue = fabs(x[j]);
		sign = x[j] > 0 ? 1 : -1;
		floorValue = floor(absValue);
		fractionalPart = absValue - floorValue;
		if (1 - fractionalPart < eps) {
			x[j] = sign * (floorValue + 1);
			x[j] = sign * (floorValue + 1);
			x[j] = sign * (floorValue + 1);
			x[j] = sign * (floorValue + 1);
			continue;
		}
		if (fractionalPart < eps)
			x[j] = sign * floorValue;
	}
}

inline void SF_Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
	for (int j = 0; j < PD_n; j++)
		z[j] = x[j] - y[j];
}

inline void SF_Vector_Zeroing(PT_vector_T x) {  // x = 0
	for (int j = 0; j < PD_n; j++) x[j] = 0;
}

//---------------------------------- Private Functions -------------------------

// Not used