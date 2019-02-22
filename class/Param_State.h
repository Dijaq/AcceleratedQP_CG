#ifndef PARAM_STATE_H
#define PARAM_STATE_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

class Param_State{
public:

	Eigen::MatrixXd V;
	Eigen::MatrixXd F;
	Eigen::MatrixXd eq_lhs;
	Eigen::MatrixXd eq_rhs;

	Param_State();
};

Param_State::Param_State()
{}

#endif