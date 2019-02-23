#ifndef SOLVER_PROJECTOR_BD_H
#define SOLVER_PROJECTOR_BD_H

class SolverProjectorBD{

public:
	//Problem
	MatrixXd K;
	MatrixXd lb;
	MatrixXd ub;
	int dim;
	int distortions;
	bool flips;
	int minsv;
	int maxsv;

	SolverProjectorBD();
	~SolverProjectorBD();
}

SolverProjectorBD::SolverProjectorBD(){
	
}

SolverProjectorBD::~SolverProjectorBD(){}


#endif