/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: BIP (Block-lterative Projection) method (MPI)
Module: Problem-bsfCode.cpp (Implementation of Problem Code)
Prefix:	PC_bsf	- BSF Predefined Problem Functions
		SF		- Shared Functionc
		PF		- Private Functions
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;
using namespace SF;

//---------------------------------- BSF Predefined Problem Functions -----------------------------

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	Vector_Copy(parameterIn.x, parameterOutP->x);
}

void PC_bsf_Init(bool* success) {
	PD_problemName = PP_PROBLEM_NAME;
	PD_m = 0;
	PD_n = 0;

	if (!EpsilonsAreOK(PP_EPS_ZERO, PP_EPS_PROJECTION, PP_EPS_ON_HYPERPLANE)) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "PC_bsf_Init error: The following condition must be satisfied:\n"
			<< "PP_EPS_ON_HYPERPLANE >  PP_EPS_PROJECTION > PP_EPS_ZERO > DBL_EPSILON = "
			<< DBL_EPSILON << endl;
		*success = false;
		return;
	}

#ifdef PP_MPS_FORMAT
	*success = MPS___Load_Problem();
#else
	*success = MTX__Load_Problem();
#endif // PP_MPS_FORMAT

	if (*success == false)
		return;

	MakeColumnOfNorms(PD_A, PD_norm_a);

	for (int i = 0; i < PD_m; i++)
		if (PD_norm_a[i] < PP_EPS_ZERO) { //Degenerate equation!
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "PC_bsf_Init error: Equation/inequality " << i << " is degenerate with pecision of PP_EPS_ZERO = "
				<< PP_EPS_ZERO << "!!!\n";
			*success = false;
			return;
		}

	// Number of constraints being equations
	PD_meq = 0;
	for (int i = 0; i < PD_m; i++)
		if (PD_isEquation[i])
			PD_meq++;
	assert(PD_meq < PD_n);

	PD_neq = PD_n - PD_meq;	// Dimension of the subspace of intersection of equation-hyperplanes

	*success = MTX_LoadPoint(PD_x, PP_MTX_POSTFIX_X0);
	if (*success == false)
		return;

	MakeColumnOfNorms(PD_A, PD_norm_a);

	for (int i = 0; i < PD_m; i++)
		if (PD_norm_a[i] < PP_EPS_ZERO) { //Degenerate equation!
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Equation/inequality " << i << " is degenerate with pecision of PP_EPS_ZERO = "
				<< PP_EPS_ZERO << "!!!\n";
			*success = false;
			return;
		}

	if (PointBelongsToPolytope(PD_x, PP_EPS_ON_HYPERPLANE)) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster) {
			cout << "Starting point is feasible.";
#ifdef PP_SAVE_RESULT
			MTX_SavePoint(PD_x, PP_MTX_POSTFIX_Z0);
			cout << " It is saved in file *_z0.mtx\n";
#			endif // PP_SAVE_RESULT
		}
		*success = false;
		return;
	}
	PD_iterNo = 0;
}

void PC_bsf_IterInit(PT_bsf_parameter_T parameter) {
	// Not used
}

void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {

	cout << "# " << BSF_sv_iterCounter << "\t Time " << round(elapsedTime);
	//cout << "\tx ="; Print_Vector(parameter.x);
	cout << "\tF(t) = " << setprecision(24) << ObjF(parameter.x) << "\tDistance: " << Distance_PointToPolytope(parameter.x)
		<< setprecision(PP_SETW / 2) << endl;
#ifdef PP_SAVE_ITER_IN_X0
	if (MTX_SavePoint(parameter.x, PP_MTX_POSTFIX_X0))
		cout << ". x is saved into x0." << endl;
#endif // PP_SAVE_ITER_IN_X0
}

void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
}

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
	// PF_MAP_LIST_INDEX
	OrthogonalProjectingVectorOntoHalfspace_i(BSF_sv_parameter.x, mapElem->constraint_i, reduceElem->projectingVector, success);
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

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "=================================================== " << PP_METHOD_NAME << " ====================================================" << endl;
	cout << "Problem name: " << PD_problemName << endl;

#ifdef PP_MPS_FORMAT
	cout << "Input format: MPS" << endl;
	cout << "m = " << PD_m << "\tn = " << PD_n << " (after conversion into standard form)" << endl;
#else
	cout << "Input format: MTX (with elimination of free variables)" << endl;
	cout << "Before elimination: m =\t" << PP_M << "\tn = " << PP_N << endl;
	cout << "After elimination:  m =\t" << PD_m << "\tn = " << PD_n << endl;
#endif // PP_MPS_FORMAT

	cout << "Number of equations: " << PD_meq << endl;
	if (PD_meq > 0)
		cout << "Subspace dimension: " << PD_neq << endl;

	cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;

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
	cout << "Map List is Fragmented" << endl;
#else
	cout << "Map List is not Fragmented" << endl;
#endif

	cout << "PP_EPS_ZERO\t\t" << PP_EPS_ZERO << endl;
	cout << "PP_EPS_PROJECTION\t" << PP_EPS_PROJECTION << endl;
	cout << "PP_EPS_ON_HYPERPLANE\t" << PP_EPS_ON_HYPERPLANE << endl;
	cout << "--------------- Data ---------------\n";
#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix PD_A & Column PD_b -------" << endl;
	Print_Constraints();
#endif // PP_MATRIX_OUTPUT
	cout << "Obj Function:\t";
	Print_Vector(PD_c);
	cout << endl;
	cout << "x0 =\t\t";
	Print_Vector(PD_x); cout << "\tF(x) = " << setw(PP_SETW) << ObjF(PD_x) << endl;

#ifdef PP_DEBUG
	cout << "x0 on hyperplanes: "; Print_HyperplanesIncludingPoint(PD_x, PP_EPS_ON_HYPERPLANE);
#endif // PP_DEBUG

	cout << "\n-------------------------------------------" << endl;
}

void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	cout << setprecision(24);
	cout << "================================================" << endl;
	cout << "// Elapsed time: " << t << endl;
	cout << "// Number of iterations: " << PD_iterNo << endl;
	cout << "// Objective function: " << ObjF(parameter.x) << endl;
	cout << "// Distance to polytope: " << Distance_PointToPolytope(parameter.x) << endl;
	cout << "// z on hyperplanes: "; Print_HyperplanesIncludingPoint(parameter.x, PP_EPS_ON_HYPERPLANE); cout << endl;
	cout << "// Number of including inequality hyperplanes: " << Number_IncludingNeHyperplanes(parameter.x, PP_EPS_ON_HYPERPLANE) << endl;
	cout << "================================================" << endl;

#ifdef PP_SAVE_RESULT
	if (MTX_SavePoint(parameter.x, PP_MTX_POSTFIX_Z0))
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

	if (reduceCounter == 0) {
		*exit = true;
		return;
	}

	/*DEBUG PC_bsf_ProcessResults**
	#ifdef PP_DEBUG
	cout << "z = "; Print_Vector(parameter->x); cout << endl;
	#endif // PP_DEBUG /**/

#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "-------------> PC_bsf_ProcessResults: Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = "
			<< PP_MAX_ITER_COUNT << endl;
		*exit = true;
		return;
	};
#endif // PP_MAX_ITER_COUNT

	Vector_DivideEquals(reduceResult->projectingVector, (double)(reduceCounter));

#ifdef PP_DEBUG
	PT_vector_T x_prev;
	Vector_Copy(parameter->x, x_prev);
#endif

	Vector_Addition(parameter->x, reduceResult->projectingVector, parameter->x);

#ifdef PP_DEBUG
	double dist = Distance_PointToPoint(parameter->x, x_prev);
	if (dist > 0)
		if (dist < DBL_EPSILON) {
			cout << "PC_bsf_ProcessResults error: The distance between the approximations is less than the machine epsilon! You should increase PP_EPS_PROJECTION." << endl;
			*exit = true;
			return;
		}
#endif

	/*DEBUG PC_bsf_ProcessResults**
#ifdef PP_DEBUG
	cout << "_______________________________ " << PD_iterNo << " ___________________________________\n";
	cout << "r = "; Print_Vector(reduceResult->projectingVector);
	cout << "\t||r|| = " << Vector_Norm(reduceResult->projectingVector) << endl;
	cout << "x = "; Print_Vector(parameter->x); cout << endl;
	cout << "x on hyperplanes: "; Print_HyperplanesIncludingPoint(parameter->x, PP_EPS_ON_HYPERPLANE); cout << endl;
#endif // PP_DEBUG /**/

	double length = Vector_Norm(reduceResult->projectingVector);

	if (length < PP_EPS_PROJECTION) {
		*exit = true;
		return;
	}

#ifdef PP_DEBUG
	static bool notFirsIteration;
	static PT_vector_T projectingVector_prev;

	if (notFirsIteration) {
		if (Vector_Equal(reduceResult->projectingVector, projectingVector_prev)) {
			cout << "PC_bsf_ProcessResults error: The significand bit depth has been exceeded! You should increase PP_EPS_PROJECTION.\n";
			*exit = true;
			return;
		}
	}
	else
		notFirsIteration = true;

	Vector_Copy(reduceResult->projectingVector, projectingVector_prev);
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
	Vector_Addition(x->projectingVector, y->projectingVector, z->projectingVector);
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
	Vector_Copy(PD_x, parameter->x);
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PD_m;
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->constraint_i = i;
}

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; }
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; }
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; }
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; }
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; }
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; }
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; }
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; }

//---------------------------------- Shared Functions -------------------------
namespace SF {

	static inline unsigned long long BinomialCoefficient(int n, int k) { // https://habr.com/ru/articles/274689/
		unsigned long long res = 1;
		if (k > n / 2) k = n - k;
		if (k == 1)  return n;
		if (k == 0)  return 1;
		for (int i = 1; i <= k; ++i) {
			res *= n - k + i;
			res /= i;
		}
		return res;
	}

	static inline void Bitscale_Create(bool* bitscale, int m, int* hyperplanes, int mh) {
		for (int i = 0; i < m; i++)
			bitscale[i] = false;
		for (int ih = 0; ih < mh; ih++)
			bitscale[hyperplanes[ih]] = true;
	}

	static inline double Distance_PointToHalfspace_i(PT_vector_T x, int i) {
		double a_DoT_z_MinuS_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];

		if (a_DoT_z_MinuS_b < 0) // Point belongs to halfspace
			return 0;

		return a_DoT_z_MinuS_b / PD_norm_a[i];
	}

	static inline double Distance_PointToHyperplane_i(PT_vector_T x, int i) {
		return fabs(Vector_DotProduct(PD_A[i], x) - PD_b[i]) / PD_norm_a[i];
	}

	static inline double Distance_PointToPoint(PT_vector_T x, PT_vector_T y) {
		PT_vector_T z;
		Vector_Subtraction(x, y, z);
		return Vector_Norm(z);
	}

	static inline double Distance_PointToPolytope(PT_vector_T x) { // Measure of distance from point to polytope
		double maxDistance = 0;
		double distance;

		for (int i = 0; i < PD_m; i++) {
			if (PD_isEquation[i])
				distance = Distance_PointToHyperplane_i(x, i);
			else
				distance = Distance_PointToHalfspace_i(x, i);
			if (distance > 0)
				maxDistance = PF_MAX(maxDistance, distance);
		}
		return maxDistance;
	}

	static inline double DistanceSQR_PointToPoint(PT_vector_T x, PT_vector_T y) {
		PT_vector_T z;
		Vector_Subtraction(x, y, z);
		return Vector_NormSquare(z);
	}

	static inline bool EpsilonsAreOK(double eps_zero, double eps_projection, double eps_on_hyperplane) {
		return (eps_zero > DBL_EPSILON && eps_projection > eps_zero && eps_on_hyperplane > eps_projection);
	}

	static inline void Flat_BipProjection(int* flatHyperplanes, int m_flat, PT_vector_T v, double eps_projection, int maxProjectingIter, PT_vector_T w, int* success) {
		PT_vector_T p;
		PT_vector_T r;
		int iterCount = 0;
		double length_r;

		Vector_Copy(v, w);
		*success = true;

		do {
			Vector_Zeroing(r);
			for (int i = 0; i < m_flat; i++) {
				OrthogonalProjectingVectorOntoHyperplane_i(w, flatHyperplanes[i], p);
				Vector_PlusEquals(r, p);
			}

			Vector_DivideEquals(r, m_flat);

#ifdef PP_DEBUG
			PT_vector_T w_prev;
			Vector_Copy(w, w_prev);
#endif // PP_DEBUG

			Vector_PlusEquals(w, r);

#ifdef PP_DEBUG
			double dist = Distance_PointToPoint(w, w_prev);
			if (dist < DBL_EPSILON * 10) { // Significand bit depth is exceeded!
				*success = -1;
				return;
			}
#endif // PP_DEBUG

			if (iterCount > maxProjectingIter) {
				*success = -2;
				break;
			}

			length_r = Vector_Norm(r);

			/*DEBUG Flat_BipProjection**
			#ifdef PP_DEBUG
			if (iterCount % PP_PROJECTION_COUNT == 0)
				cout << "Worker " << BSF_sv_mpiRank << ": \t sublist_index = " << BSF_sv_numberInSublist << ": \t Length of r = " << length_r << endl;
			#endif // PP_DEBUG /**/

		} while (length_r >= eps_projection);

		/*DEBUG PC_bsf_MapF**
		#ifdef PP_DEBUG
		cout << "Flat_BipProjection: iterCount = " << iterCount << endl;
		#endif // PP_DEBUG /**/
	}

	static inline void Flat_MaxProjection(int* flatHyperplanes, int m_flat, PT_vector_T v, double eps_projection, int maxProjectingIter, PT_vector_T w, int* success) {
		PT_vector_T p;
		PT_vector_T p_max;
		double max_length;
		int iterCount = 0;

		Vector_Copy(v, w);
		*success = true;

		do {
			max_length = 0;
			Vector_Zeroing(p_max);
			for (int i = 0; i < m_flat; i++) {
				OrthogonalProjectingVectorOntoHyperplane_i(w, flatHyperplanes[i], p);
				double norm_p = Vector_Norm(p);
				if (norm_p > max_length) {
					max_length = norm_p;
					Vector_Copy(p, p_max);
				}
			}

			Vector_PlusEquals(w, p_max);

			iterCount++;
			if (iterCount > maxProjectingIter) {
				*success = -2;
				break;
			}

			/*DEBUG Flat_MaxProjection**
			#ifdef PP_DEBUG
			if (iterCount % PP_PROJECTION_COUNT == 0)
				cout << "Worker " << BSF_sv_mpiRank << ": \tsublist_index = " << BSF_sv_numberInSublist << ": \t max_length = " << max_length << endl;
			#endif // PP_DEBUG /**/

		} while (max_length >= eps_projection);

		/*DEBUG Flat_MaxProjection**
		#ifdef PP_DEBUG
		cout << "Flat_MaxProjection: iterCount = " << iterCount << endl;
		#endif // PP_DEBUG /**/
	}

	static inline void JumpingOnPolytope(PT_vector_T startPoint, PT_vector_T direcionVector, PT_vector_T finishPoint, double eps_on_hyperplane, double eps_zero, bool* parallelHPlanes, int* success) {
		PT_vector_T o;		// Oblique projection vector
		PT_vector_T o_min;	// Oblique projection vector with minimum length
		double length_o;
		double a_DoT_d;
		double norm_d;
		double norm_a_DoT_norm_d;
		int location_z;
		double a_DoT_z_MinuS_b;
		double* d = direcionVector;		// Direction vector
		double* z = startPoint;
		double minLength_o = PP_INFINITY;

		*success = true;

		/*DEBUG JumpingOnPolytope**
		#ifdef PP_DEBUG
		cout << "d =\t"; Print_Vector(d); cout << endl;
		#endif // PP_DEBUG /**/

		norm_d = Vector_Norm(d);
		if (norm_d < eps_zero) {
			/*DEBUG JumpingOnPolytope**
			#ifdef PP_DEBUG
			cout << "Worker " << BSF_sv_mpiRank << ": JumpingOnPolytope: norm_d < eps_zero => return \n";
			#endif // PP_DEBUG /**/
			Vector_Copy(startPoint, finishPoint);
			*success = false;
			return;
		}

		Vector_Zeroing(o_min);

		for (int i = 0; i < PD_m; i++) {
			if (parallelHPlanes[i])
				continue;

			a_DoT_d = Vector_DotProduct(PD_A[i], d);
			norm_a_DoT_norm_d = a_DoT_d / (PD_norm_a[i] * norm_d);

			location_z = PointLocation_i(z, i, eps_on_hyperplane, &a_DoT_z_MinuS_b);
			switch (location_z) {
			case PP_ON_HYPERPLANE:

				if (fabs(norm_a_DoT_norm_d) < eps_zero) { // Vector d is parallel to hyperplane
					/*DEBUG JumpingOnPolytope**
					#ifdef PP_DEBUG
					cout << "Worker " << BSF_sv_mpiRank << ": JumpingOnPolytope: " << i
						<< ") \tStart and finish points belong to hyperplane. => continue\n";
					#endif // PP_DEBUG /**/
					continue;
				}

				if (norm_a_DoT_norm_d >= eps_zero) {
					/*DEBUG JumpingOnPolytope**
					#ifdef PP_DEBUG
					cout << "Worker " << BSF_sv_mpiRank << ": JumpingOnPolytope: " << i
						<< ") \tStart point belong to hyperplane, finish point is outside half-space. => return\n";
					#endif // PP_DEBUG /**/
					Vector_Copy(startPoint, finishPoint);
					*success = false;
					return;
				}

			case PP_INSIDE_HALFSPACE:

				if (fabs(norm_a_DoT_norm_d) < eps_zero) { // Vector d is parallel to hyperplane
					/*DEBUG JumpingOnPolytope**
				#ifdef PP_DEBUG
				cout << "Worker " << BSF_sv_mpiRank << ": JumpingOnPolytope: " << i << ") \tVector d is parallel to hyperplane. => continue\n";
				#endif // PP_DEBUG /**/
					continue;
				}

				if (norm_a_DoT_norm_d < 0) { // Vector looks inside of half-space
					/*DEBUG JumpingOnPolytope**
					#ifdef PP_DEBUG
					cout << "Worker " << BSF_sv_mpiRank << ": JumpingOnPolytope: " << i << ") \tVector looks inside of half-space. => continue\n";
					#endif // PP_DEBUG /**/
					continue;
				}

				// norm_a_DoT_norm_d > 0	// Vector looks out from half-space
				// Oblique projection vector: o = -(<a,z> - b)d/<a, d>
				Vector_MultiplyByNumber(d, -a_DoT_z_MinuS_b / a_DoT_d, o);
				length_o = Vector_Norm(o);
				if (minLength_o > length_o) {
					minLength_o = length_o;
					Vector_Copy(o, o_min);
				}
				/*DEBUG JumpingOnPolytope**
				#ifdef PP_DEBUG
				PT_vector_T u_next;
				Vector_Addition(startPoint, o, u_next);
				cout << "Worker " << BSF_sv_mpiRank << ": " << i << ") \tJump: length_o = " << length_o << ":\tu_next = "; Print_Vector(u_next);
				cout << "\tObjF(u_next) = " << ObjF(u_next) << endl;
				#endif // PP_DEBUG /**/
				break;
			case PP_OUTSIDE_HALFSPACE:
				cout << "JumpingOnPolytope error: Point is outside halfspace!\n";
				assert(false);
			default:
				assert(false);
			}
		}
		Vector_Addition(startPoint, o_min, finishPoint);
	}

	static inline void MakeColumnOfNorms(PT_matrix_T A, PT_column_T norm_a) {
		for (int i = 0; i < PD_m; i++)
			norm_a[i] = Vector_Norm(A[i]);
	}

	static inline void MakeListOfNotIncludingHalfspaces(PT_vector_T x, int* notIncludingHalfspacesList, double eps_on_hyperplane) {
		int mo = 0;
		for (int i = 0; i < PD_m; i++)
			if (!PointBelongsToHalfspace_i(x, i, eps_on_hyperplane)) {
				notIncludingHalfspacesList[mo] = i;
				mo++;
			}
		if (mo < PD_m)
			notIncludingHalfspacesList[mo] = -1;
	}

	static inline void MakeNeHyperplaneList(PT_vector_T u, int* neHyperplanes_u, int* mneh_u, double eps_on_hyperplane) {
		// List of hyperplanes that are not equations and include point u.
		*mneh_u = 0;
		for (int i = 0; i < PD_m; i++) {
			if (PD_isEquation[i])
				continue;
			if (PointBelongsToHyperplane_i(u, i, eps_on_hyperplane)) {
				neHyperplanes_u[*mneh_u] = i;
				(*mneh_u)++;
			}
		}
	}

	static bool MPS___Load_Problem() {
		const char* mpsFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MPS_file;
		string word;
		int n_row = 0;
		int n_col = 0;
		int n_up = 0;
		int n_fx = 0;

		PT_MPS_row_T* rows = (PT_MPS_row_T*)calloc(PP_MAX_NUMBER_OF_ROWS, sizeof(PT_MPS_row_T));
		if (rows == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem error: Can't allocate memory for array 'rows'.\n";
			return false;
		}

		PT_MPS_column_T* columns = (PT_MPS_column_T*)calloc(PP_MAX_NUMBER_OF_COLS, sizeof(PT_MPS_column_T));
		if (rows == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem error: Can't allocate memory for array 'columns'.\n";
			return false;
		}

		double* loBounds = (double*)calloc(PP_N, sizeof(double));
		if (rows == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem error: Can't allocate memory for array 'loBounds'.\n";
			return false;
		}

		PT_MPS_upBound_T* upBounds = (PT_MPS_upBound_T*)calloc(PP_N, sizeof(PT_MPS_upBound_T));
		if (rows == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem error: Can't allocate memory for array 'upBounds'.\n";
			return false;
		}

		PT_MPS_fxVariable_T* fxVariables = (PT_MPS_fxVariable_T*)calloc(PP_N, sizeof(PT_MPS_fxVariable_T));
		if (rows == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem error: Can't allocate memory for array 'fxVariables'.\n";
			return false;
		}

		MPS_file = PP_PATH;
		MPS_file += PP_MPS_PREFIX;
		MPS_file += PD_problemName;
		MPS_file += PP_MPS_EXTENSION;
		mpsFile = MPS_file.c_str();
		stream = fopen(mpsFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mpsFile << "'.\n";
			return false;
		}

		MPS__SkipComment(stream);

		if (!MPS__ReadKeyWord(stream, &word, "NAME")) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem: Syntax error '" << word << "', expected 'NAME'.\n";
			return false;
		}

		MPS__SkipComment(stream);

		if (!MPS__ReadKeyWord(stream, &word, "ROWS")) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem: Syntax error '" << word << "', expected 'ROWS'.\n";
			return false;
		}

		if (!MPS__ReadRows(stream, rows, &n_row))
			return false;

		MPS__SkipComment(stream);

		if (!MPS__ReadKeyWord(stream, &word, "COLUMNS")) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem: Syntax error '" << word << "', expected 'COLUMNS'.\n";
			return false;
		}

		if (!MPS__ReadColumns(stream, columns, &n_col))
			return false;

		MPS__SetColumnIndexes(columns, n_col, &PD_n);

		MPS__SkipComment(stream);

		if (!MPS__ReadKeyWord(stream, &word, "RHS")) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem: Syntax error '" << word << "', expected 'RHS'.\n";
			return false;
		}

		if (!MPS__ReadRHS(stream, rows, n_row))
			return false;

		MPS__SkipComment(stream);

		if (MPS__ReadKeyWord(stream, &word, "BOUNDS"))
			if (!MPS__ReadBounds(stream, columns, n_col, loBounds, upBounds, &n_up, fxVariables, &n_fx))
				return false;

		MPS__SkipComment(stream);

		if (!MPS__ReadKeyWord(stream, &word, "ENDATA")) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem error: '" << word << "', expected 'ENDATA'.\n";
			return false;
		}

		if (!MPS__MakeProblem(rows, n_row, columns, n_col, loBounds, upBounds, n_up, fxVariables, n_fx))
			return false;

		free(rows);
		free(columns);
		free(loBounds);
		free(upBounds);
		return true;
	}

	static inline bool MPS__MakeProblem(PT_MPS_row_T* row, int n_row, PT_MPS_column_T* column, int n_col, double* loBound,
		PT_MPS_upBound_T* upBound, int n_up, PT_MPS_fxVariable_T* fxVariable, int n_fx)
	{
		PD_m = 0;

		for (int i = 0; i < PP_MM; i++) {
			for (int j = 0; j < PP_N; j++)
				PD_A[i][j] = 0;
			PD_b[i] = 0;
			PD_isEquation[i] = false;
		}

		for (int j = 0; j < PP_N; j++)
			PD_c[j] = 0;

		for (int i_row = 0; i_row < n_row; i_row++) {
			switch (row[i_row].type) {
			case 'E':
				if (!MPS_AddEquation(row[i_row].name, row[i_row].RHS_value, column, n_col))
					return false;
				break;
			case 'G':
				if (!MPS_AddInequality_G(row[i_row].name, row[i_row].RHS_value, column, n_col))
					return false;
				break;
			case 'L':
				if (!MPS_AddInequality_L(row[i_row].name, row[i_row].RHS_value, column, n_col))
					return false;
				break;
			case 'N':
				if (!MPS_AddObjectiveFunction(row[i_row].name, column, n_col))
					return false;
				break;
			default:
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__MakeProblem error:Unpredictable row type " << row[i_row].type << endl;
				return false;
			}
		}

		if (PD_m != PP_M) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS__MakeProblem error: Number of constraints in mps-file = " << PD_m << " not equal to PP_M = " << PP_M << ".\n";
			return false;
		}

		if (PD_n != PP_N) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS__MakeProblem error: Number of variables in mps-file = " << PD_n << " not equal to PP_M = " << PP_N << ".\n";
			return false;
		}

		for (int j = 0; j < PD_n; j++) { // Adding lower bounds
			PD_A[PD_m + j][j] = -1;
			if (loBound[j] > 0)
				PD_b[PD_m + j] = -loBound[j];
		}
		PD_m += PD_n;

		for (int j_up = 0; j_up < n_up; j_up++) { // Adding upper bounds
			PD_A[PD_m][upBound[j_up].varIndex] = 1;
			PD_b[PD_m] = upBound[j_up].value;
			PD_m++;
			assert(PD_m <= PP_MM);
		}

		for (int j_fx = 0; j_fx < n_fx; j_fx++) { // Adding fixed variables
			PD_A[PD_m][fxVariable[j_fx].varIndex] = 1;
			PD_b[PD_m] = fxVariable[j_fx].value;
			PD_isEquation[PD_m] = true;
			PD_m++;
			assert(PD_m <= PP_MM);
		}

		return true;
	}

	static inline bool MPS__ReadBounds(FILE* stream, PT_MPS_column_T* column, int n_col, double* loBound,
		PT_MPS_upBound_T* upBound, int* n_up, PT_MPS_fxVariable_T* fxVariable, int* n_fx)
	{
		char ch, typeCh1, typeCh2;
		PT_MPS_name_T baundName;
		PT_MPS_name_T varName;
		int j_var;
		double boundValue;
		fpos_t pos;	// Position in the input stream

		assert(*n_up == 0);
		assert(*n_fx == 0);

		for (int j_lo = 0; j_lo < PD_n; j_lo++)
			loBound[j_lo] = 0;

		ch = getc(stream);
		if (ch != ' ') {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS__ReadBounds error: Expected 'space'!\n";
			return false;
		}

		while (ch == ' ') {
			typeCh1 = getc(stream);
			typeCh2 = getc(stream);
			ch = getc(stream);

			MPS_ReadName(stream, baundName);
			MPS_SkipSpaces(stream);
			MPS_ReadName(stream, varName);

			j_var = -1;
			for (int j_col = 0; j_col < n_col; j_col++)
				if (MPS_SameNames(column[j_col].varName, varName)) {
					j_var = column[j_col].j;
					break;
				}
			if (j_var == -1) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__ReadBounds Error: Variable " << varName << " was not found.\n";
				return false;
			}

			MPS_SkipSpaces(stream);

			MPS_ReadValue(stream, &boundValue);

			switch (typeCh1) {
			case 'L':
				loBound[j_var] = boundValue;
				break;
			case 'U':
				upBound[*n_up].value = boundValue;
				upBound[*n_up].varIndex = j_var;
				(*n_up)++;
				break;
			case 'F':
				fxVariable[*n_fx].value = boundValue;
				fxVariable[*n_fx].varIndex = j_var;
				(*n_fx)++;
				break;
			default:
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__ReadBounds Error: unpredictable bound type " << typeCh1 << endl;
				return false;
			}

			MPS_NewLine(stream);
			fgetpos(stream, &pos);
			ch = getc(stream);
		}

		fsetpos(stream, &pos);
		return true;
	}

	static inline bool MPS__ReadColumns(FILE* stream, PT_MPS_column_T* column, int* n_col) {
		char ch;
		fpos_t pos;	// Position in the input stream

		*n_col = 0;

		fgetpos(stream, &pos);
		ch = getc(stream);
		fsetpos(stream, &pos);
		if (ch != ' ') {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS__ReadColumns error: Expected 'space'!\n";
			return false;
		}

		while (ch == ' ') {
			fgetpos(stream, &pos);
			ch = getc(stream);
			fsetpos(stream, &pos);
			if (ch != ' ')
				return true;

			if (!MPS_ReadColumnLine(stream, column, n_col))
				return false;
		}

		return true;
	}

	static inline bool MPS__ReadKeyWord(FILE* stream, string* word, string pattern) {
		char str[80] = { '\0' };
		fpos_t pos;	// Position in the input stream

		fgetpos(stream, &pos);

		if (fscanf(stream, "%s", str) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "MPS__ReadKeyWord error: Failure to read string!\n";
			return false;
		}

		*word = str;

		if (*word != pattern) {
			fsetpos(stream, &pos);
			return false;
		}

		MPS_NewLine(stream);

		return true;
	}

	static inline bool MPS__ReadRHS(FILE* stream, PT_MPS_row_T* row, int n_row) {
		PT_MPS_name_T RHS_name = { '\0' };
		fpos_t pos;	// Position in the input stream
		char ch;

		for (int i_row = 0; i_row < n_row; i_row++)
			row[i_row].RHS_value = 0;

		fgetpos(stream, &pos);
		ch = getc(stream);
		if (ch == 'B') {
			fsetpos(stream, &pos);
			return true;
		}
		while (ch == ' ') {
			fsetpos(stream, &pos);
			MPS_ReadRHS_line(stream, row, n_row, RHS_name);

			fgetpos(stream, &pos);
			ch = getc(stream);
		}
		fsetpos(stream, &pos);
		return true;
	}

	static inline bool MPS__ReadRows(FILE* stream, PT_MPS_row_T* row, int* n_row) {
		char ch;
		int j;
		fpos_t pos;	// Position in the input stream

		*n_row = 0;

		fgetpos(stream, &pos);
		ch = getc(stream);
		while (ch == ' ') {
			j = 1;

			ch = getc(stream);
			j++;

			if (ch == ' ') {
				ch = getc(stream);
				j++;
			}

			if (ch <= ' ') {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__ReadRows: Syntax error in row " << *n_row + 1 << endl;
				return false;
			}

			switch (ch) {
			case 'E':
			case 'G':
			case 'L':
			case 'N':
				row[*n_row].type = ch;
				break;
			default:
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__ReadRows: Syntax error - unpredictable row type " << ch << endl;
				return false;
			}

			while (j < 4) {
				ch = getc(stream);
				if (ch < ' ') {
					if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
						cout << "MPS__ReadRows: Syntax error - not ASCII symbol.\n";
					return false;
				}
				j++;
			}

			MPS_SkipSpaces(stream);

			MPS_ReadName(stream, row[*n_row].name);

			(*n_row)++;
			if (*n_row > PP_MAX_NUMBER_OF_ROWS) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__ReadRows: The number of row is greater than PP_MAX_NUMBER_OF_ROWS = " << PP_MAX_NUMBER_OF_ROWS << endl;
				return false;
			}

			if (!MPS_UniqueRowName(row, *n_row, row[*n_row - 1].name)) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__ReadRows:Error: The name '" << row[*n_row - 1].name << "' is not unique.";
				return false;
			}

			MPS_NewLine(stream);
			fgetpos(stream, &pos);
			ch = getc(stream);
		}

		fsetpos(stream, &pos);
		return true;
	}

	static inline void MPS__SkipComment(FILE* stream) {
		fpos_t pos;	// Position in the input stream
		int res;

		fgetpos(stream, &pos);

		while (getc(stream) == '*') {
			while (getc(stream) != 10);
			res = fscanf(stream, "\n");
			fgetpos(stream, &pos);
		}

		fsetpos(stream, &pos);
	}

	static inline bool MPS_AddEquation(PT_MPS_name_T rowName, double RHS_value, PT_MPS_column_T* column, int n_col) {
		bool empty = true;

		for (int i_col = 0; i_col < n_col; i_col++)
			if (MPS_SameNames(column[i_col].rowName, rowName)) {
				if (!PD_A[PD_m][column[i_col].j] == 0) {
					if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
						cout << "MPS_AddEquation error: Coefficient redefinition of the variable " << column[i_col].varName
						<< " in row " << column[i_col].rowName << ".\n";
					return false;
				}
				PD_A[PD_m][column[i_col].j] = column[i_col].value;
				empty = false;
			}
		if (empty) {
			if (RHS_value != 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_AddEquation error: Zero row " << rowName << " has non-zero RHS value.\n";
				return false;
			}
			return true;
		}
		PD_b[PD_m] = RHS_value;
		PD_isEquation[PD_m] = true;
		PD_m++;
		assert(PD_m < PP_MM);
		return true;
	}

	static inline bool MPS_AddInequality_G(PT_MPS_name_T rowName, double RHS_value, PT_MPS_column_T* column, int n_col) {
		bool empty = true;

		for (int i_col = 0; i_col < n_col; i_col++)
			if (MPS_SameNames(column[i_col].rowName, rowName)) {
				if (!PD_A[PD_m][column[i_col].j] == 0) {
					if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
						cout << "MPS_AddInequality_G error: Coefficient redefinition of the variable " << column[i_col].varName
						<< " in row " << column[i_col].rowName << ".\n";
					return false;
				}
				PD_A[PD_m][column[i_col].j] = -column[i_col].value;
				empty = false;
			}
		if (empty) {
			if (RHS_value != 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_AddInequality_G error: Zero row " << rowName << " has non-zero RHS value.\n";
				return false;
			}
			return true;
		}
		PD_b[PD_m] = -RHS_value;
		PD_m++;
		assert(PD_m < PP_MM);
		return true;
	}

	static inline bool MPS_AddInequality_L(PT_MPS_name_T rowName, double RHS_value, PT_MPS_column_T* column, int n_col) {
		bool empty = true;

		for (int i_col = 0; i_col < n_col; i_col++)
			if (MPS_SameNames(column[i_col].rowName, rowName)) {
				if (!PD_A[PD_m][column[i_col].j] == 0) {
					if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
						cout << "MPS_AddInequality_L error: Coefficient redefinition of the variable " << column[i_col].varName
						<< " in row " << column[i_col].rowName << ".\n";
					return false;
				}
				PD_A[PD_m][column[i_col].j] = column[i_col].value;
				empty = false;
			}
		if (empty) {
			if (RHS_value != 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_AddInequality_L error: Zero row " << rowName << " has non-zero RHS value.\n";
				return false;
			}
			return true;
		}
		PD_b[PD_m] = RHS_value;
		PD_m++;
		assert(PD_m <= PP_MM);
		return true;
	}

	static inline bool MPS_AddObjectiveFunction(PT_MPS_name_T rowName, PT_MPS_column_T* column, int n_col) {
		for (int i_col = 0; i_col < n_col; i_col++)
			if (MPS_SameNames(column[i_col].rowName, rowName)) {
				if (!PD_c[column[i_col].j] == 0) {
					if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
						cout << "MPS_AddObjectiveFunction error: Coefficient redefinition of the objective function for " << column[i_col].varName << ".\n";
					return false;
				}
#ifdef MPS_MIN_OF_OBJECTIVE_FUNCTION
				PD_c[column[i_col].j] = column[i_col].value;
#else
				PD_c[column[i_col].j] = -column[i_col].value;
#endif
			}
		return true;
	}

	static inline void MPS_CopyName(char* name_x, char* name_y) {
		for (int p = 0; p < PP_MPS_NAME_LENGTH; p++)
			name_y[p] = name_x[p];
	}

	static inline void MPS_NewLine(FILE* stream) {
		char ch;
		fpos_t pos;	// Position in the input stream

		fgetpos(stream, &pos);
		if (getc(stream) == -1) // EOF
			return;
		fsetpos(stream, &pos);

		do {
			fgetpos(stream, &pos);
		} while (getc(stream) != '\n');

		fsetpos(stream, &pos);

		do {
			fgetpos(stream, &pos);
			ch = getc(stream);
		} while (ch == '\n' || ch == '\r');

		fsetpos(stream, &pos);
	}

	static inline bool MPS_ReadColumnLine(FILE* stream, PT_MPS_column_T* column, int* n_col) {
		fpos_t pos;	// Position in the input stream
		double value;

		MPS_SkipSpaces(stream);

		if (!MPS_ReadName(stream, column[*n_col].varName))
			return false;

		MPS_SkipSpaces(stream);

		if (!MPS_ReadName(stream, column[*n_col].rowName))
			return false;

		// Checking for duplicates
		for (int j_col = 0; j_col < *n_col; j_col++)
			if (MPS_SameNames(column[j_col].varName, column[*n_col].varName) && MPS_SameNames(column[j_col].rowName, column[*n_col].rowName)) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_ReadColumnLine error: Redefinition of item varName = " << column[j_col].varName
					<< ", rowName = " << column[j_col].rowName << ".\n";
				return false;
			}

		MPS_SkipSpaces(stream);

		if (!MPS_ReadValue(stream, &value))
			return false;

		column[*n_col].value = value;

		(*n_col)++;
		if (*n_col > PP_MAX_NUMBER_OF_COLS) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadColumnLine error: The number of column is greater than PP_MAX_NUMBER_OF_COLS = " << PP_MAX_NUMBER_OF_COLS << endl;
			return false;
		}

		MPS_SkipSpaces(stream);

		fgetpos(stream, &pos);
		if (getc(stream) < ' ') {
			fsetpos(stream, &pos);
			MPS_NewLine(stream);
			return true;
		};
		fsetpos(stream, &pos);

		MPS_CopyName(column[*n_col - 1].varName, column[*n_col].varName);

		if (!MPS_ReadName(stream, column[*n_col].rowName))
			return false;

		// Checking for duplicates
		for (int j_col = 0; j_col < *n_col; j_col++)
			if (column[j_col].varName == column[*n_col].varName && column[j_col].rowName == column[*n_col].rowName) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_ReadColumnLine error: Redefinition of item varName = " << column[j_col].varName
					<< ", rowName = " << column[j_col].rowName << ".\n";
				return false;
			}

		MPS_SkipSpaces(stream);

		if (!MPS_ReadValue(stream, &value))
			return false;

		column[*n_col].value = value;

		(*n_col)++;
		if (*n_col > PP_MAX_NUMBER_OF_COLS) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS__ReadColumns: The number of column is greater than PP_MAX_NUMBER_OF_COLS = " << PP_MAX_NUMBER_OF_COLS << endl;
			return false;
		}

		MPS_NewLine(stream);

		return true;
	}

	static inline bool MPS_ReadName(FILE* stream, char* name) {
		char ch;
		fpos_t pos;	// Position in the input stream

		for (int p = 0; p < PP_MPS_NAME_LENGTH; p++)
			name[p] = '\0';

		fgetpos(stream, &pos);
		if (getc(stream) < ' ') {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadName: Syntax error - not ASCII symbol." << endl;
			return false;
		}
		fsetpos(stream, &pos);

		MPS_SkipSpaces(stream);

		for (int p = 0; p < 8; p++) {
			fgetpos(stream, &pos);
			ch = getc(stream);
			if (ch <= ' ') {
				fsetpos(stream, &pos);
				break;
			}
			name[p] = ch;
		}

		return true;
	}

	static inline bool MPS_ReadRHS_line(FILE* stream, PT_MPS_row_T* row, int n_row, PT_MPS_name_T RHS_name)
	{ // RHS_name must be equal to {'\0'} when calling MPS_ReadRHS_line at the first time
		fpos_t pos;	// Position in the input stream
		char ch;
		PT_MPS_name_T next_RHS_name;
		PT_MPS_name_T rowName;
		double RHS_value;
		int rowIndex;

		for (int p = 0; p < PP_MPS_NAME_LENGTH; p++)
			next_RHS_name[p] = '\0';

		for (int p = 0; p < 4; p++) {
			ch = getc(stream);
			if (ch != ' ') {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_ReadRHS_line: Syntax error, expected ' '\n";
				return false;
			}
		}

		int p = 0;
		fgetpos(stream, &pos);
		while (ch == ' ') {
			fgetpos(stream, &pos);
			ch = getc(stream);
			p++;
		}
		fsetpos(stream, &pos);

		if (p > 8)
			next_RHS_name[0] = ' ';
		else
			if (!MPS_ReadName(stream, next_RHS_name))
				return false;

		if (RHS_name[0] != '\0') {
			if (!MPS_SameNames(RHS_name, next_RHS_name)) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_ReadRHS_line: Syntax error, Multiple RHS names.\n";
				return false;
			}
		}
		else
			MPS_CopyName(next_RHS_name, RHS_name);

		MPS_SkipSpaces(stream);

		if (!MPS_ReadName(stream, rowName))
			return false;

		rowIndex = MPS_SearchRowByName(row, n_row, rowName);

		if (rowIndex < 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Syntax error, non-existent row name '" << rowName << "'.\n";
			return false;
		}

		if (row[rowIndex].type == 'N') {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Invalid row type 'N'.\n";
			return false;
		}

		MPS_SkipSpaces(stream);

		if (fscanf(stream, "%lf", &RHS_value) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Unexpected end of line!\n";
			return false;
		}
		row[rowIndex].RHS_value = RHS_value;

		MPS_SkipSpaces(stream);

		fgetpos(stream, &pos);
		if (getc(stream) < ' ') {
			fsetpos(stream, &pos);
			MPS_NewLine(stream);
			fgetpos(stream, &pos);
			ch = getc(stream);
			fsetpos(stream, &pos);
			return true;
		};
		fsetpos(stream, &pos);

		if (!MPS_ReadName(stream, rowName))
			return false;

		rowIndex = MPS_SearchRowByName(row, n_row, rowName);

		if (rowIndex < 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Syntax error, non-existent row name '" << rowName << "'.\n";
			return false;
		}

		if (row[rowIndex].type == 'N') {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Invalid row type 'N'.\n";
			return false;
		}

		MPS_SkipSpaces(stream);

		if (fscanf(stream, "%lf", &RHS_value) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Unexpected end of line!\n";
			return false;
		}
		row[rowIndex].RHS_value = RHS_value;

		MPS_SkipSpaces(stream);

		fgetpos(stream, &pos);
		if (char ch = getc(stream) > ' ') {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Illegal ASCII character '" << ch << "' at end of line.\n";
			return false;
		}
		fsetpos(stream, &pos);

		MPS_NewLine(stream);

		return true;
	}

	static inline bool MPS_ReadValue(FILE* stream, double* value) {
		if (fscanf(stream, "%lf", value) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadValue: Error: Non-ASCII character.\n";
			return false;
		}
		return true;
	}

	static inline bool MPS_SameNames(PT_MPS_name_T name_x, PT_MPS_name_T name_y) {
		for (int p = 0; p < PP_MPS_NAME_LENGTH; p++) {
			if (name_x[p] == '\0' && name_y[p] == '\0')
				return true;
			if (name_x[p] != name_y[p])
				return false;
		}
		return true;
	}

	static inline void MPS__SetColumnIndexes(PT_MPS_column_T* column, int n_col, int* n) {

		*n = 0;

		for (int j_col = 0; j_col < n_col; j_col++)
			column[j_col].j = -1;

		for (int j_col = 0; j_col < n_col; j_col++) {
			if (column[j_col].j >= 0)
				continue;
			column[j_col].j = *n;
			(*n)++;
			assert(*n <= PP_N);
			for (int j = j_col + 1; j < n_col; j++) {
				if (column[j].j >= 0)
					continue;
				if (MPS_SameNames(column[j_col].varName, column[j].varName))
					column[j].j = column[j_col].j;
			}
		}
	}

	static inline int MPS_SearchColByName(PT_MPS_column_T* column, int j_col, PT_MPS_name_T name) {
		for (int j = 0; j < j_col; j++)
			if (MPS_SameNames(column[j].varName, name))
				return j;
		return -1;
	}

	static inline int MPS_SearchRowByName(PT_MPS_row_T* row, int n_row, PT_MPS_name_T name) {
		for (int i = 0; i < n_row; i++)
			if (MPS_SameNames(row[i].name, name))
				return i;
		return -1;
	}

	static inline void MPS_SkipSpaces(FILE* stream) {
		fpos_t pos;	// Position in the input stream

		do {
			fgetpos(stream, &pos);
		} while (getc(stream) == ' ');

		fsetpos(stream, &pos);
	}

	static inline bool MPS_UniqueRowName(PT_MPS_row_T* rows, int n_row, PT_MPS_name_T name) {
		int count = 0;

		for (int i = 0; i < n_row; i++)
			if (MPS_SameNames(rows[i].name, name)) {
				if (count > 0)
					return false;
				count = 1;
			}
		return true;
	}

	static void MTX_Conversion() { // Removing free variables

		for (int i = 0; i < PD_m; i++)
			PD_isEquation[i] = true;

		MTX_RemoveFreeVariables();

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

		/**
		cout << "-----------------------------------------------------\n";
		Print_Constraints();
		cout << "-----------------------------------------------------\n";
		cout << "PD_c: "; Print_Vector(PD_c); cout << endl;/**/
	}

	static bool MTX__Load_Problem() {

		//--------------- Reading A ------------------
		if (!MTX_Load_A())
			return false;

		//--------------- Reading b ------------------
		if (!MTX_Load_b())
			return false;

		//--------------- Reading lo ------------------
		if (!MTX_Load_lo())
			return false;

		//--------------- Reading c ------------------
		if (!MTX_Load_c())
			return false;

		//--------------- Reading hi ------------------
		if (!MTX_Load_hi())
			return false;

		//---------- Conversion to inequalities -----------
		MTX_Conversion();

		/**
		cout << "-----------------------------------------------------\n";
		Print_Constraints();
		cout << "-----------------------------------------------------\n";
		cout << "PD_c: "; Print_Vector(PD_c); cout << endl;/**/

		return true;
	}

	static inline bool MTX_Load_A() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		int non;	// Number of non-zero elements
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += PP_MTX_POSTFIX_A;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		MTX_SkipComments(stream);
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
				cout << "MTX_Load_A error: PP_N must be = " << noc << "\n";
			return false;
		}

		if (nor != PP_M) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MTX_Load_A error:  PP_M must be = " << nor << "\n";
			return false;
		}

		PD_m = nor;
		PD_n = noc;

		if (nor + noc > PP_MM) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Invalid input data: number of constraints m = " << nor + noc
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

	static inline bool MTX_Load_b() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += PP_MTX_POSTFIX_B;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		MTX_SkipComments(stream);
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

	static inline bool MTX_Load_c() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += PP_MTX_POSTFIX_C;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		MTX_SkipComments(stream);
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

	static inline bool MTX_Load_hi() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += PP_MTX_POSTFIX_HI;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		MTX_SkipComments(stream);
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

	static inline bool MTX_Load_lo() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += PP_MTX_POSTFIX_LO;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		MTX_SkipComments(stream);
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

	static inline bool MTX_LoadPoint(PT_vector_T x, string postfix) {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += postfix;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'.\n";
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}
		if (nor != PD_n) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of rows in'" << mtxFile << "'. Must be " << PD_n << "\n";
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

	static inline void MTX_RemoveFreeVariables(void) {
		int freeVariable_j_col[PP_M];
		int m_freeVariable_j_col = 0;
		bool unique;

		// Find free variables
		for (int j_col = 0; j_col < PD_n; j_col++) { // Find zero element in PD_c
			if (PD_c[j_col] != 0)
				continue;
			if (PD_lo[j_col] != 0)
				continue;
			if (PD_hi[j_col] != PP_INFINITY)
				continue;

			for (int i_row = 0; i_row < PD_m; i_row++) { // Find PD_A i_row with nonzero element in column j_col
				if (PD_A[i_row][j_col] == 0)
					continue;

				// Check uniqueness in column j_col
				unique = true;
				for (int i = 0; i < PD_m; i++) {
					if (i == i_row)
						continue;
					if (PD_A[i][j_col] != 0) {
						unique = false;
						break;
					}
				}
				if (!unique)
					continue;

				freeVariable_j_col[m_freeVariable_j_col] = j_col;
				m_freeVariable_j_col++; assert(m_freeVariable_j_col <= PP_M);

				// Check sign of free variable
				if (PD_A[i_row][j_col] < 0) {
					// Change sign of inequality
					for (int j = 0; j < PD_n; j++)
						PD_A[i_row][j] = -PD_A[i_row][j];
					PD_b[i_row] = -PD_b[i_row];
				}

				PD_isEquation[i_row] = false;

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
					Print_Constraints();
					cout << "-----------------------------------------------------\n";/**/
				}
			}
		}
	}

	static bool MTX_SavePoint(PT_vector_T x, string postfix) {
		const char* mtxFile;
		FILE* stream;// Input stream
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += postfix;
		mtxFile = MTX_file.c_str();
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

	static inline void MTX_SkipComments(FILE* stream) {
		fpos_t pos;	// Position in the input stream
		int res;
		res = fscanf(stream, "\n");
		fgetpos(stream, &pos);
		while (getc(stream) == '%') {
			while (getc(stream) != 10);
			res = fscanf(stream, "\n");
			fgetpos(stream, &pos);
		}
		fsetpos(stream, &pos);
	}

	static inline int Number_IncludingNeHyperplanes(PT_vector_T x, double eps_on_hyperplane) {
		int number = 0;

		for (int i = 0; i < PD_m; i++) {
			if (PD_isEquation[i])
				continue;
			if (PointBelongsToHyperplane_i(x, i, eps_on_hyperplane))
				number++;
		}
		return number;
	}

	static inline double ObjF(PT_vector_T x) {
		double s = 0;
		for (int j = 0; j < PD_n; j++)
			s += PD_c[j] * x[j];
		return s;
	}

	static inline void OrthogonalProjectingVectorOntoHalfspace_i(PT_vector_T z, int i, PT_vector_T r, int* success) {
		double factor;
		double a_DoT_z_MinuS_b = Vector_DotProduct(PD_A[i], z) - PD_b[i]; // <a,z> - b

		if (!PD_isEquation[i])
			if (a_DoT_z_MinuS_b <= 0) { // <a,z> - b < 0
				*success = false;
				return;
			}

		factor = -a_DoT_z_MinuS_b / (PD_norm_a[i] * PD_norm_a[i]); // (b - <z,a>) / ||a||^2
		Vector_MultiplyByNumber(PD_A[i], factor, r); // r = a(b - <z,a>) / ||a||^2
		*success = true;
	}

	static inline void OrthogonalProjectingVectorOntoHyperplane_i(PT_vector_T x, int i, PT_vector_T p) {
		double ns = Vector_NormSquare(PD_A[i]);
		Vector_MultiplyByNumber(PD_A[i], -(Vector_DotProduct(PD_A[i], x) - PD_b[i]) / ns, p);
	}

	static inline bool PointBelongsToFlat(PT_vector_T x, int* hyperplaneList, int hyperplaneCount, double eps_on_hyperplane) { // If the point belongs to the flat
		for (int i = 0; i < hyperplaneCount; i++)
			if (!PointBelongsToHyperplane_i(x, hyperplaneList[i], eps_on_hyperplane))
				return false;
		return true;
	}

	static inline bool PointBelongsToHalfspace_i(PT_vector_T x, int i, double eps_on_hyperplane) {
		double a_DoT_x_MinuS_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];
		double distanceToHyperplane = fabs(a_DoT_x_MinuS_b) / PD_norm_a[i];
		if (distanceToHyperplane < eps_on_hyperplane)
			return true;
		if (PD_isEquation[i])
			return false;
		if (a_DoT_x_MinuS_b < 0)
			return true;
		return false;
	}

	static inline bool PointBelongsToHyperplane_i(PT_vector_T x, int i, double eps_on_hyperplane) {
		double dist = Distance_PointToHyperplane_i(x, i);

		/*DEBUG PointBelongsToHyperplane_i**
#ifdef PP_DEBUG
		if (dist > eps_on_hyperplane && dist < eps_on_hyperplane * 10) {
			cout << "Distance from testing point is less than " << PD_eps_on_hyperplane*10 << ", but greater than " << PD_eps_on_hyperplane << "!\n";
			//system("pause");
		}
#endif // PP_DEBUG /**/

		if (dist < eps_on_hyperplane)
			return true;
		else
			return false;
	}

	static inline bool PointBelongsToPolytope(PT_vector_T x, double eps_on_hyperplane) { // If the point belongs to the polytope with prescigion of eps_on_hyperplane
		for (int i = 0; i < PD_m; i++)
			if (!PointBelongsToHalfspace_i(x, i, eps_on_hyperplane))
				return false;
		return true;
	}

	static inline bool PointInsideHalfspace_i(PT_vector_T x, int i, double eps_on_hyperplane) {
		double a_DoT_x_MinuS_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];
		double distanceToHyperplane = fabs(a_DoT_x_MinuS_b) / PD_norm_a[i];
		if (distanceToHyperplane < eps_on_hyperplane)
			return false;
		if (a_DoT_x_MinuS_b < 0)
			return true;
		return false;
	}

	static inline bool PointIsBoundary(PT_vector_T x, double eps_on_hyperplane) {
		if (!PointBelongsToPolytope(x, eps_on_hyperplane))
			return false;

		for (int i = 0; i < PD_m; i++) {
			if (PD_isEquation[i])
				continue;
			if (PointBelongsToHyperplane_i(x, i, eps_on_hyperplane))
				return true;
		}
		return false;
	}

	static inline bool PointIsVertex(PT_vector_T x, double eps_on_hyperplane) {
		int count_x = 0; // The number of inequality hyperplanes including the point x

		for (int i = 0; i < PD_m; i++) {
			if (PD_isEquation[i])
				continue;
			if (PointBelongsToHyperplane_i(x, i, eps_on_hyperplane))
				count_x++;
		}

		if (count_x >= PD_neq)
			return true;
		else
			return false;
	}

	static inline int PointLocation_i(PT_vector_T x, int i, double eps_on_hyperplane, double* a_DoT_x_MinuS_b) {
		*a_DoT_x_MinuS_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];
		double dist = fabs(*a_DoT_x_MinuS_b) / PD_norm_a[i];

		if (dist < eps_on_hyperplane)// <a,x> = b
			return PP_ON_HYPERPLANE;

		if (*a_DoT_x_MinuS_b < 0)								// <a,x> < b
			return PP_INSIDE_HALFSPACE;

		return PP_OUTSIDE_HALFSPACE;							// <a,x> > b
	}

	static inline void PolytopeHomothety(PT_vector_T center, double ratio) { // https://en.wikipedia.org/wiki/Homothety
		if (ratio == 1)
			return;
		assert(ratio > 0);

		for (int i = 0; i < PD_m; i++) {
			PD_b[i] = ratio * PD_b[i] - (ratio - 1) * Vector_DotProduct(PD_A[i], center);
		}
	}

	static inline void Print_Constraints() {
		for (int i = 0; i < PD_m; i++) {
			cout << i << ")";
			for (int j = 0; j < PD_n; j++)
				cout << setw(PP_SETW) << PD_A[i][j];
			cout << "\t" << (PD_isEquation[i] ? "==" : "<=") << setw(PP_SETW) << PD_b[i] << endl;
		}
	}

	static inline void Print_HalfspacesIncludingPoint(PT_vector_T x, double eps_on_hyperplane) {
		bool comma = false;

		cout << "{";

		for (int i = 0; i < PD_m; i++) {
			if (PointBelongsToHalfspace_i(x, i, eps_on_hyperplane)) {
				if (comma)
					cout << ", ";
				else
					comma = true;
				cout << i;
			}
		}

		cout << "}";
	}

	static inline void Print_HalfspacesOutOfPoint(PT_vector_T x, double eps_on_hyperplane) {
		bool comma = false;

		cout << "{";

		for (int i = 0; i < PD_m; i++) {
			if (!PointBelongsToHalfspace_i(x, i, eps_on_hyperplane)) {
				if (comma)
					cout << ", ";
				else
					comma = true;
				cout << i;
			}
		}

		cout << "}";
	}

	static inline void Print_HyperplanesIncludingPoint(PT_vector_T x, double eps_on_hyperplane) {
		bool comma = false;

		cout << "{";

		for (int i = 0; i < PD_m; i++) {
			if (PointBelongsToHyperplane_i(x, i, eps_on_hyperplane)) {
				if (comma)
					cout << ", ";
				else
					comma = true;
				cout << i;
			}
		}

		cout << "}";
	}

	static inline void Print_Number_of_edges(PT_vector_T x, double eps_on_hyperplane) {
		int mne;
		unsigned long long ull_mne;

		mne = 0;
		for (int i = 0; i < PD_m; i++) {
			if (PD_isEquation[i])
				continue;
			if (PointBelongsToHyperplane_i(x, i, eps_on_hyperplane))
				mne++;
		}

		if (mne == PD_neq)
			ull_mne = (unsigned long long) mne;
		else {
			if (mne > 62) {
				cout << "Warning: Can't calculate binomial coefficient for number of including hyperplanes mne = "
					<< mne << " > 62" << endl;
				return;
			}
			ull_mne = BinomialCoefficient(mne, PD_neq - 1);
		}
		cout << ull_mne << endl;
	}

	static inline void Print_Vector(PT_vector_T x) {
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << x[j];
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	}

	static inline double RelativeError(double trueValue, double calculatedValue) {
		if (trueValue == 0)
			return fabs(calculatedValue - trueValue);
		else
			return fabs(calculatedValue - trueValue) / fabs(trueValue);
	}

	static inline void Tuning_Eps_PointBelongsToFlat(PT_vector_T x, int* hyperplaneList, int hyperplaneCount, double* eps) {
		while (!PointBelongsToFlat(x, hyperplaneList, hyperplaneCount, *eps))
			(*eps) *= 2;
	}

	static inline void Tuning_Eps_PointBelongsToPolytope(PT_vector_T x, double* eps) {
		while (!PointBelongsToPolytope(x, *eps))
			(*eps) *= 2;
	}

	static inline void TWIDDLE // https://doi.org/10.1145/362384.362502
	(int* x, int* y, int* z, int* p, bool* done) {
		int i, j, k;
		j = 0;
		*done = false;

		do {
			j++;
		} while (p[j] <= 0);

		if (p[j - 1] == 0) {
			i = j - 1;
			while (i != 1) {
				p[i] = -1;
				i -= 1;
			}
			p[j] = 0;
			p[1] = *x = *z = 1;
			*y = j;
			return;
		}

		if (j > 1)
			p[j - 1] = 0;

		do {
			j++;
		} while (p[j] > 0);

		i = k = j - 1;

		i++;
		while (p[i] == 0) {
			p[i] = -1;
			i++;
		}

		if (p[i] == -1) {
			p[i] = *z = p[k];
			*x = i;
			*y = k;
			p[k] = -1;
			return;
		}

		if (i == p[0]) {
			*done = true;
			return;
		}

		*z = p[j] = p[i];
		p[i] = 0;
		*x = j;
		*y = i;
	}

	static inline void TWIDDLE_CodeToSubset(int code, int* a, int* c, int n, int m, int* x, int* y, int* z, int* p, bool* done, int* nextI) {
		if (*nextI == 0) {
			for (int k = 0; k < m; k++)
				c[k] = a[n - m + k];
			if (code == 0) {
				(*nextI)++;
				return;
			}
		}

		do {
			TWIDDLE(x, y, z, p, done);
			assert(!*done);
			c[*z - 1] = a[*x - 1];
			(*nextI)++;
		} while (*nextI < code);
	}

	static inline void TWIDDLE_Make_p(int* p, int n, int m) {
		// p - auxiliary integer array for generating all combinations of m out of n objects.
		assert(n >= m && m > 0);
		p[0] = n + 1;
		p[n + 1] = -2;
		for (int j = 1; j <= n - m; j++)
			p[j] = 0;
		for (int j = n - m + 1; j <= n; j++)
			p[j] = j - n + m;
	}

	static inline void Shift(PT_vector_T point, PT_vector_T shiftVector, double factor, PT_vector_T shiftedPoint) {
		for (int j = 0; j < PD_n; j++)
			shiftedPoint[j] = point[j] + shiftVector[j] * factor;
	}

	static inline void Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x + y
		for (int j = 0; j < PD_n; j++)
			z[j] = x[j] + y[j];
	}

	static inline void Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
		for (int j = 0; j < PD_n; j++)
			toPoint[j] = fromPoint[j];
	}

	static inline void Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = x/r
		for (int j = 0; j < PD_n; j++)
			y[j] = x[j] / r;
	}

	static inline void Vector_DivideEquals(PT_vector_T x, double r) {  // x = x/r
		for (int j = 0; j < PD_n; j++)
			x[j] /= r;
	}

	static inline double Vector_DotProduct(PT_vector_T x, PT_vector_T y) {
		double sum = 0;
		for (int j = 0; j < PD_n; j++)
			sum += x[j] * y[j];
		return sum;
	}

	static inline bool Vector_Equal(PT_vector_T x, PT_vector_T y) { // x = y
		for (int j = 0; j < PD_n; j++)
			if (x[j] != y[j])
				return false;
		return true;
	}

	static inline void Vector_MakeLike(PT_vector_T x, double lengthOfLikeVector, PT_vector_T likeVector) {
		double norm_x = Vector_Norm(x);
		if (norm_x == 0)
			Vector_Zeroing(likeVector);
		else
			Vector_MultiplyByNumber(x, lengthOfLikeVector / norm_x, likeVector);
	}

	static inline void Vector_MakeMinus_e(PT_vector_T minus_e) {
		for (int j = 0; j < PD_n; j++)
			minus_e[j] = -1;
	}

	static inline void Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector) { // equalPoint += minusVector
		for (int j = 0; j < PD_n; j++)
			equalPoint[j] -= minusVector[j];
	}

	static inline void Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = r*x
		for (int j = 0; j < PD_n; j++)
			y[j] = x[j] * r;
	}

	static inline void Vector_MultiplyEquals(PT_vector_T x, double r) {  // x = r*x
		for (int j = 0; j < PD_n; j++)
			x[j] *= r;
	}

	static inline double Vector_Norm(PT_vector_T x) {
		double norm_x = sqrt(Vector_NormSquare(x));
		return norm_x;
	}

	static inline double Vector_NormSquare(PT_vector_T x) {
		double sum = 0;

		for (int j = 0; j < PD_n; j++) {
			sum += x[j] * x[j];
		}
		return sum;
	}

	static inline void Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector) { // equalVector += plusVector
		for (int j = 0; j < PD_n; j++)
			equalVector[j] += plusVector[j];
	}

	static inline void Vector_Round(PT_vector_T x, double eps_round) {
		double floorValue;
		double fractionalPart;
		double sign;
		double absValue;

		if (eps_round == 0)
			return;

		for (int j = 0; j < PD_n; j++) {
			if (fabs(x[j]) < eps_round) {
				x[j] = 0;
				continue;
			}
			absValue = fabs(x[j]);
			sign = x[j] > 0 ? 1 : -1;
			floorValue = floor(absValue);
			fractionalPart = absValue - floorValue;
			if (1 - fractionalPart < eps_round) {
				x[j] = sign * (floorValue + 1);
				continue;
			}
			if (fractionalPart < eps_round)
				x[j] = sign * floorValue;
		}
	}

	static inline void Vector_SetValue(PT_vector_T x, double v) {  // x = (v,...,v)
		for (int j = 0; j < PD_n; j++) x[j] = v;
	}

	static inline void Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
		for (int j = 0; j < PD_n; j++)
			z[j] = x[j] - y[j];
	}

	static inline void Vector_Zeroing(PT_vector_T x) {  // x = 0
		for (int j = 0; j < PD_n; j++) x[j] = 0;
	}

}

//---------------------------------- Private Functions -------------------------
// Not used