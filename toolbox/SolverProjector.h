#ifndef SOLVER_PROJECTOR_H
#define SOLVER_PROJECTOR_H

class SolverProjector{

public:
	//Problem
	MatrixXd T;
	MatrixXd W;
	MatrixXd eqLHS;
	MatrixXd eqRHS;
	VectorXd x0;

	//Solver
	int mode;
	bool usePreFactorization = true;
	int nVars;
	int nEq;
	int x;
	int pTx;
	int tanNormal;
	int tanLHS;
	int tanRHS;
	int preFactorization;
	MatrixXd TWW;
	MatrixXd TWWT;

	//Temporaries
	int verbose = 4;
	int t_iter;
	int t_projectD;
	int t_projectLinear;
	int t_factorization;

	//Aux
	int lambdaMultiTangent = 10;

	SolverProjector();
	~SolverProjector();
}

SolverProjector::SolverProjector(){}
SolverProjector::~SolverProjector(){}


#endif