#pragma once
#include <Eigen/Core>
#include <Eigen/LU>
#include "GetState.h"
EIGEN_NO_DEBUG

#define STATE_DIM 4 //Dimension of state vector
#define INPUT_DOF 2 //Dimension of input vector

class MPC
{
public:
	MPC();
	~MPC();

	Eigen::Matrix<double, STATE_DIM, 1> xf_end;//Target delivery position

	void setDestination(const Eigen::Matrix<double, STATE_DIM, 1>& x_end);//Set the target position
	void calcStateFunc(Eigen::Matrix<double, STATE_DIM, 1>& f, const Eigen::Matrix<double, STATE_DIM, 1>& cur_x, const Eigen::Matrix<double, INPUT_DOF, 1>& cur_u);//Calculate state equation
	int calcInputbyGradientMethod(Eigen::Matrix<double, INPUT_DOF, 1>& input, const Eigen::Matrix<double, STATE_DIM, 1>& cur_x, double T, double dt, pos3d_t LShould, pos3d_t RShould, pos3d_t Body, double length_arm, double* BodyValues);//Calculate the input by gradient method
	void calcGradient(Eigen::VectorXd& F, Eigen::VectorXd& x, const Eigen::VectorXd& cur_x, const Eigen::VectorXd& cur_u, const int N, const double dt, pos3d_t LShould, pos3d_t RShould, pos3d_t Body, double length_arm, double*BodyValues);//Calculate the gradient vector
	void JointAngle2EndPosition(Eigen::Vector2d& pos, const Eigen::Vector2d& angle);                         //Forward Kinematics to find position of the EndEffector
	void JointVel2EndVel(Eigen::Vector2d& vel, const Eigen::Vector2d& angle, const Eigen::Vector2d& anvel);  //Use Jacobian to get EndEffector velocity

	double calcVisibCost(const Eigen::Vector2d& angle, double* BodyValues); 
	double calcSafeCost(const Eigen::Vector2d& angle, double* BodyValues);
	double calcArmComfortR(pos3d_t ShouldR, double Arm, double Bangle, pos3d_t EndEff);
	pos4d_t InvKineArmR(pos3d_t EndEff_SR, pos3d_t ShouldR, double Arm, double phi);
	pos4d_t JointDispR(pos4d_t Joint_value, pos4d_t Max, pos4d_t Min);
	double calcSmoothCost(const Eigen::Matrix<double, STATE_DIM, 1>& cur_x,const Eigen::Matrix<double, INPUT_DOF, 1>& cur_u, const int N);

private:

	static const double EPS;//Small value
	static const int CYCLE_NEWTON;//Calculation cycle of gradient method

	static const double R_C1;//Weight coefficients for Barrier function
	static const double R_C2;
	static const double R_C3;
	static const double R_C4;


	Eigen::Matrix<double, STATE_DIM, STATE_DIM> S;//Weight coefficients matrix

	Eigen::Matrix<double, 2, 2> M;  //Inertia Matrix
	Eigen::Matrix<double, 2, 1> C;  //Coriolis term

	pos4d_t Max_R, Min_R, Min_L, Max_L;
	
	void calcInertiaMat(const Eigen::Vector2d& q); //Calculate Inertia Matrix
	void calcCoriolisTerm(const Eigen::Vector2d& q, const Eigen::Vector2d& dq);//Calculate the Coriolis term
	double calcObjective(const Eigen::VectorXd& x_seq, const Eigen::VectorXd& u_seq, const int N, pos3d_t LShould, pos3d_t RShould, pos3d_t Body, double length_arm, double* BodyValues);                          //Calculate the total cost function
	double calcObjectivePhi(const Eigen::Matrix<double, STATE_DIM, 1>& goal_x, const int N);//Calculate the cost function Phi
	double calcObjectiveB(const Eigen::Matrix<double, STATE_DIM, 1>& cur_x, const Eigen::Matrix<double, INPUT_DOF, 1>& cur_u);//Calculate the cost function B
	double calcHamiltonian(const Eigen::Matrix<double, STATE_DIM, 1>& cur_x, const Eigen::Matrix<double, INPUT_DOF, 1>& cur_u, const Eigen::Matrix<double, STATE_DIM, 1>& cur_lmd, pos3d_t LShould, pos3d_t RShould, pos3d_t Body, double length_arm, double* BodyValues, const int N);//Calculate the Hamiltonian function
	void calcdHdu(Eigen::Matrix<double, 1, INPUT_DOF>& dHdu, const Eigen::Matrix<double, STATE_DIM, 1>& cur_x, const Eigen::Matrix<double, STATE_DIM, 1>& cur_lmd);//Partial differentiate the Hamiltonian function by u
	double calcObjectivePhiDash(const Eigen::Matrix<double, STATE_DIM, 1>& goal_x, const int N);//Calculate the cost function PhiDash, Velocity instead of position
};

