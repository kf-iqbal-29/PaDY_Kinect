#include "ArmParams.h"
#include "MPC.h"
#include<sstream>
#include<fstream>
#include <iostream>

//#define DEBUG //For debug
#define SAMPLING_TIME 30 //(ms) 

//////* Definition of static member variables *//////
const double MPC::EPS = 0.001;//Small value
const int MPC::CYCLE_NEWTON = 200;//Calculation cycle of gradient method
const double MPC::R_C1 = 1000.0;     //Weight coefficients for barrier function
const double MPC::R_C2 = 1000.0;     
const double MPC::R_C3 = 1000.0;     
const double MPC::R_C4 = 1000.0;  




MPC::MPC()
{   /* Initialize weighting matrix */ 
	S.setZero();
	S(0, 0) = 0.0;
	S(1, 1) = 0.0;
	S(2, 2) = 5.0;
	S(3, 3) = 5.0;
}


MPC::~MPC()
{
}


//Set the target position
void MPC::setDestination(const Eigen::Matrix<double, STATE_DIM, 1>& x_end)
{
	xf_end = x_end;
}

//Calculate the end position from joint angles
void MPC::JointAngle2EndPosition(Eigen::Vector2d& pos, const Eigen::Vector2d& angle)
{
	double c1, c12, s1, s12;

	if (fabs(cos(angle(0))) < 0.001) c1 = 0;
	else c1 = cos(angle(0));
	if (fabs(cos((angle(0) + angle(1)))) < 0.001) c12 = 0;
	else c12 = cos((angle(0) + angle(1)));
	if (fabs(sin(angle(0))) < 0.001) s1 = 0;
	else s1 = sin(angle(0));
	if (fabs(sin((angle(0) + angle(1)))) < 0.001) s12 = 0;
	else s12 = sin((angle(0) + angle(1)));

	pos(0) = L1*c1 + L2*c12;
	pos(1) = L1*s1 + L2*s12;
}

//Calculate the velocity of endpoint from the joint velocities
void MPC::JointVel2EndVel(Eigen::Vector2d& vel, const Eigen::Vector2d& angle, const Eigen::Vector2d& anvel)
{
	double c1, c12, s1, s12;

	if (fabs(cos(angle(0))) < 0.001) c1 = 0;
	else c1 = cos(angle(0));
	if (fabs(cos((angle(0) + angle(1)))) < 0.001) c12 = 0;
	else c12 = cos((angle(0) + angle(1)));
	if (fabs(sin(angle(0))) < 0.001) s1 = 0;
	else s1 = sin(angle(0));
	if (fabs(sin((angle(0) + angle(1)))) < 0.001) s12 = 0;
	else s12 = sin((angle(0) + angle(1)));

	vel(0) = ((-L1*s1) + (-L2*s12)) * anvel(0) + (-L2*s12) * anvel(1);
	vel(1) = ((L1*c1) + (L2*c12)) * anvel(0) + (L2*c12) * anvel(1);
}

//Calculate Inertia matrix
void MPC::calcInertiaMat(const Eigen::Vector2d& q)
{
	M(0, 0) = M1*LG1*LG1 + M2*L1*L1 + M2*LG2*LG2 + INERTIA1 + INERTIA2 + 2 * M2*L1*LG2*cos(q(1));
	M(0, 1) = M2*LG2*LG2 + INERTIA2 + M2*L1*LG2*cos(q(1));
	M(1, 0) = M(0, 1);
	M(1, 1) = M2*LG2*LG2 + INERTIA2;
}

//Calculate the Coriolis term
void MPC::calcCoriolisTerm(const Eigen::Vector2d& q, const Eigen::Vector2d& dq)
{
	C(0) = -M2*L1*LG2*(2 * dq(0) + dq(1))*dq(1)*sin(q(1));
	C(1) = M2*L1*LG2*dq(0)*dq(0)*sin(q(1));
}


//Calculate state equation
void MPC::calcStateFunc(Eigen::Matrix<double, STATE_DIM, 1>& f, const Eigen::Matrix<double, STATE_DIM, 1>& cur_x, const Eigen::Matrix<double, INPUT_DOF, 1>& cur_u)
{
	//Calculate matrices for Lagrangian Dynamics
	calcInertiaMat(cur_x.block<2, 1>(0, 0));                              
	calcCoriolisTerm(cur_x.block<2, 1>(0, 0), cur_x.block<2, 1>(2, 0));  

	//Calculate X_dot vector using system equations and lagrangian dynamics 
	f.block<2, 1>(0, 0) = cur_x.block<2, 1>(2, 0); 
	f.block<2, 1>(2, 0) = -M.inverse()*C + M.inverse()*cur_u;
}

//Calculate the input by gradient method
int MPC::calcInputbyGradientMethod(Eigen::Matrix<double, INPUT_DOF, 1>& input, const Eigen::Matrix<double, STATE_DIM, 1>& cur_x, double T, double dt, pos3d_t LShould, pos3d_t RShould, pos3d_t Body, double length_arm, double* BodyValues)
{
	int STEP_N;//Step size
	const double DEV_R = 0.618;      //Golden Ratio
	const double EXP_STEP = 0.05;    //Step size for linear search
	double J, pre_J, temp_J, J1, J2; //Cost values
	double scale, beta, alpha_l, alpha_u, alpha1, alpha2;
	Eigen::VectorXd pre_dir, dir;    //Descending direction of gradient		// VectorXd is dynamic vector. You can set its size later.
	Eigen::VectorXd x;              //Group of state vector
	Eigen::VectorXd u;              //Group of input vector
	Eigen::VectorXd temp_u;
	Eigen::VectorXd pre_F, F;      //Gradient of Hamiltonian function
	Eigen::Matrix<double, STATE_DIM, 1> f; //State equation
	Eigen::Vector2d temp_log; 
	double BAngle = BodyValues[0]; 

	//std::ofstream ofs2("TrajData.csv",std::ios::app);

	//Calculate the step size
	STEP_N = static_cast<int>(T / dt);
	//std::cout<<"Steps "<<STEP_N<<std::endl; 
//#ifdef DEBUG
	//std::cout << "STEP N : " << STEP_N << std::endl;
//#endif //DEBUG
	if (STEP_N < 1)
		return 0;

	//Resize vectors. At every step, we need to add new values
	pre_dir.resize(INPUT_DOF * STEP_N);
	dir.resize(INPUT_DOF * STEP_N);
	x.resize(STATE_DIM * (STEP_N + 1));
	u.resize(INPUT_DOF * STEP_N);
	temp_u.resize(INPUT_DOF * STEP_N);
	pre_F.resize(INPUT_DOF * STEP_N);
	F.resize(INPUT_DOF * STEP_N);

	//Initialization
	u.setConstant(0.0);
	pre_J = 10000000000; // assume a high cost in the beginning

	/* Start gradient method for calculating the input */
	for (int cycle = 0; cycle < CYCLE_NEWTON; ++cycle)
	{
		#ifdef DEBUG
				std::cout << "cycle : " << cycle << std::endl;
		#endif //DEBUG

		/* Calculate the gradient vector */
		calcGradient(F, x, cur_x, u, STEP_N, dt, LShould, RShould, Body, length_arm, BodyValues); // F is the output
		
		#ifdef DEBUG
				std::cout << "F jacobians are calculated" << std::endl;
		#endif //DEBUG

		/* Decide the step direction */
		if (cycle == 0)                        // Initialize variable dir in (+)Gradient direction  
			dir = F;
		else{
			//Use the conjugate gradient direction
			beta = F.norm() / pre_F.norm(); //Fletcher|Reeves method
			dir = F + beta * pre_dir;
			}

		/* Decide the step size */
		temp_J = calcObjective(x, u, STEP_N, LShould, RShould, Body, length_arm, BodyValues);
		for (int i = 0; i < static_cast<int>(1.0 / EXP_STEP); ++i){
			//Set the closed interval
			alpha_l = i * EXP_STEP;
			alpha_u = (i + 1) * EXP_STEP;
			//Update the input u
			temp_u = u - alpha_u * dir;

			//Evaluate again state vectors and put them in x vector
			x.block<STATE_DIM, 1>(0, 0) = cur_x;//Set initial state
			for (int n = 0; n < STEP_N; ++n){
				calcStateFunc(f, x.block<STATE_DIM, 1>((n * STATE_DIM), 0), temp_u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0));
				x.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0) = x.block<STATE_DIM, 1>((n * STATE_DIM), 0) + f * dt;
			}
			//Calculate again cost function
			J = calcObjective(x, temp_u, STEP_N, LShould, RShould, Body, length_arm, BodyValues);

			if (J - temp_J > 0.01){
				//std::cout<<J<<" VS "<<temp_J<<std::endl;
				//Calculate new step size using Golden Section Search method
				alpha1 = alpha_l + (1.0 - DEV_R) * (alpha_u - alpha_l);
				alpha2 = alpha_l + DEV_R * (alpha_u - alpha_l);
				for (int n = 0; n < 100; ++n){
					scale = alpha_l;
					//Stop if the size of the interval is sufficiently small
					if ((alpha_u - alpha_l) < 0.001)
						break;
					// Calculate the objective function for each step size
					/* J1 */
					temp_u = u - alpha1 * dir;//Update the input u
					//Calculte the group of state vectors by state equation
					x.block<STATE_DIM, 1>(0, 0) = cur_x;//Set the initial state
					for (int n = 0; n < STEP_N; ++n){
						calcStateFunc(f, x.block<STATE_DIM, 1>((n * STATE_DIM), 0), temp_u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0));
						x.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0) = x.block<STATE_DIM, 1>((n * STATE_DIM), 0) + f * dt;
					}
					J1 = calcObjective(x, temp_u, STEP_N, LShould, RShould, Body, length_arm, BodyValues);// Calculate the objective function
					/* J2 */
					temp_u = u - alpha2 * dir;//Update the input u
					//Calculte the group of state vectors by state equation
					x.block<STATE_DIM, 1>(0, 0) = cur_x;//Set the initial state
					for (int n = 0; n < STEP_N; ++n){
						calcStateFunc(f, x.block<STATE_DIM, 1>((n * STATE_DIM), 0), temp_u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0));
						x.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0) = x.block<STATE_DIM, 1>((n * STATE_DIM), 0) + f * dt;
					}
					J2 = calcObjective(x, temp_u, STEP_N, LShould, RShould, Body, length_arm, BodyValues);// Calculate the objective function

					//if J(alpha1) < J(alpha2)
					if (J1 < J2){
						alpha_u = alpha2;
						alpha2 = alpha1;
						alpha1 = alpha_l + (1.0 - DEV_R) * (alpha_u - alpha_l);
					}
					//if J(alpha1) >= J(alpha2)
					else{
						alpha_l = alpha1;
						alpha1 = alpha2;
						alpha2 = alpha_l + DEV_R * (alpha_u - alpha_l);
					}
				}
				break;
			}
			
			scale = alpha_l;
			temp_J = J;
		}

		/* Calculate the new input using the step size and step direction calculated by conjugate gradient direction*/
		u = u - scale * dir;
		
		/* Calculate the final result of objective function */
		x.block<STATE_DIM, 1>(0, 0) = cur_x;//Set the initial state
		for (int n = 0; n < STEP_N; ++n){
			calcStateFunc(f, x.block<STATE_DIM, 1>((n * STATE_DIM), 0), temp_u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0));
			x.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0) = x.block<STATE_DIM, 1>((n * STATE_DIM), 0) + f * dt;
		    //temp_log = x.block<2, 1>((n * STATE_DIM), 0); 
			//ofs2<< n<<","<<temp_log(0)<<","<<temp_log(1)<<std::endl; 
		}
		J = calcObjective(x, temp_u, STEP_N, LShould, RShould, Body, length_arm, BodyValues);
		
		//ofs2<<std::endl; 
		//ofs2.close(); 

#ifdef DEBUG
		//if (STEP_N == 100){
			std::cout << "scale : " << scale << "  J : " << J << std::endl;
			getchar();
		//}
#endif //DEBUG

		/* Check the convergence */
		if (((J - pre_J) > -0.001) && (cycle > 2)){
			break;
		}

		/* Store the parametes */
		pre_J = J;
		pre_dir = dir;
		pre_F = F;
	}

	//Calculte the group of state vectors by state equation
	x.block<STATE_DIM, 1>(0, 0) = cur_x;//Set the initial state
	for (int n = 0; n < STEP_N; ++n){
		calcStateFunc(f, x.block<STATE_DIM, 1>((n * STATE_DIM), 0), u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0));
		x.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0) = x.block<STATE_DIM, 1>((n * STATE_DIM), 0) + f * dt;
	}

	//Store the input 
	input = u.block<INPUT_DOF, 1>(0, 0); 

#ifdef DEBUG
		for (int n = 0; n < STEP_N; ++n){
			std::cout << "u1 : " << u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0)(0)
				<< "  u2 : " << u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0)(1) << std::endl;
		}
		for (int n = 0; n < (STEP_N + 1); ++n){
			std::cout << "q1 : " << x.block<STATE_DIM, 1>((n * STATE_DIM), 0)(0)*RAD2DEG
				<< "  q2 : " << x.block<STATE_DIM, 1>((n * STATE_DIM), 0)(1)*RAD2DEG << std::endl;
			std::cout << "qq1 : " << x.block<STATE_DIM, 1>((n * STATE_DIM), 0)(2)*RAD2DEG
				<< "  qq2 : " << x.block<STATE_DIM, 1>((n * STATE_DIM), 0)(3)*RAD2DEG << std::endl;
		}
		std::cout << std::endl;
		//getchar();
#endif //DEBUG

	return 1;
}


//Calculate the gradient vector
void MPC::calcGradient(Eigen::VectorXd& F, Eigen::VectorXd& x, const Eigen::VectorXd& cur_x, const Eigen::VectorXd& cur_u, const int N, const double dt, pos3d_t LShould, pos3d_t RShould, pos3d_t Body, double length_arm, double* BodyValues)
{
	Eigen::VectorXd lmd;                        //Adjoint variables vector lambda
	Eigen::Matrix<double, STATE_DIM, 1> f; 
	Eigen::Matrix<double, STATE_DIM, 1> dPhi;   //Partial differentiation of phi
	Eigen::Matrix<double, STATE_DIM, 1> temp_x;
	Eigen::Matrix<double, STATE_DIM, 1> dHdx;   //Partial differentiation of H by x
	Eigen::Matrix<double, INPUT_DOF, 1> temp_u;
	Eigen::Matrix<double, 1, INPUT_DOF> dHdu;   //Partial differentiation of H by u
	
	double temp;

	//Calculte the group of state vectors by state equation
	x.block<STATE_DIM, 1>(0, 0) = cur_x;    //Set the initial state as the current x
	for (int n = 0; n < N; ++n){
		calcStateFunc(f, x.block<STATE_DIM, 1>((n * STATE_DIM), 0), cur_u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0));
		x.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0) = x.block<STATE_DIM, 1>((n * STATE_DIM), 0) + f * dt;            // Discrete time update: x(n+1) = x(n) + f*dt
	}
#ifdef DEBUG
	if (N == 100){
		std::cout << "State vectors are calculated" << std::endl;
		std::cout << x << std::endl;
		getchar();
	}
#endif //DEBUG

	
	lmd.resize(STATE_DIM * (N + 1));
	//Calculate the partial differential of Phi w.r.t. x by approximation 
	temp = calcObjectivePhiDash(x.block<STATE_DIM, 1>((STATE_DIM * N), 0),N);
	for (int i = 0; i < STATE_DIM; ++i){
		temp_x = x.block<STATE_DIM, 1>((STATE_DIM * N), 0);             //Take estimated last state vector
		temp_x(i) += EPS;                                               //Incremental ratio 
		dPhi(i) = (calcObjectivePhiDash(temp_x,N) - temp) / EPS;
	}
	//Calculate costate vectors imposing first the optimal constraint: lambda(N)=dPhi/dx
	lmd.block<STATE_DIM, 1>((STATE_DIM * N), 0) = dPhi;
	for (int n = (N - 1); n > 0; --n){
		//Partial differentiate the Hamiltonian
		temp = calcHamiltonian(x.block<STATE_DIM, 1>((n * STATE_DIM), 0), cur_u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0), lmd.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0), LShould, RShould, Body, length_arm, BodyValues,N);
		for (int i = 0; i < STATE_DIM; ++i){
			temp_x = x.block<STATE_DIM, 1>((n * STATE_DIM), 0);
			temp_x(i) += EPS;
			dHdx(i) = (calcHamiltonian(temp_x, cur_u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0), lmd.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0), LShould, RShould, Body, length_arm, BodyValues,N) - temp) / EPS;
		}

		//Impose the optimal constraint: lambda(n+1) = lambda(n) - dH/dx * dt
		lmd.block<STATE_DIM, 1>((n * STATE_DIM), 0) = lmd.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0) + dHdx * dt;
	}
#ifdef DEBUG
	if (N == 100){
		std::cout << "Lambdas are calculated" << std::endl;
		std::cout << lmd << std::endl;
		getchar();
	}
#endif //DEBUG

	//Calculate the gradient vectors F
	for (int n = 0; n < N; ++n){
		//Partial differentiate the Hamiltonian
		calcdHdu(dHdu, x.block<STATE_DIM, 1>((n * STATE_DIM), 0), lmd.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0));
		//Substitute the result in F
		F.block<INPUT_DOF, 1>((n * INPUT_DOF), 0) = dHdu.transpose();
	}

#ifdef DEBUG
	//std::cout << "F equations are calculated" << std::endl;
#endif //DEBUG

}

//Calculate the the total cost function J = Phi(xN) + sumB(x)
double MPC::calcObjective(const Eigen::VectorXd& x_seq, const Eigen::VectorXd& u_seq, const int N, pos3d_t LShould, pos3d_t RShould, pos3d_t Body, double length_arm, double* BodyValues)
{
	double Phi, B, V, S, A;
	Eigen::Matrix<double, STATE_DIM / 2, 1> end_pos;
	pos3d_t EndEff;

	//Calculate the cost function PhiDash on x(N)
	Phi = calcObjectivePhiDash(x_seq.block<STATE_DIM, 1>((STATE_DIM * N), 0),N);

	//Calculate the cost function B for all the N intervals
	B = 0;
	for (int i = 0; i < N; ++i){
		B += calcObjectiveB(x_seq.block<STATE_DIM, 1>((STATE_DIM * i), 0), u_seq.block<INPUT_DOF, 1>((INPUT_DOF * i), 0));
	}
	V = 0; 
	for (int i = 0; i < N; ++i){
		V += calcVisibCost(x_seq.block<2, 1>((STATE_DIM * i), 0), BodyValues);
	}
	S = 0; 
	for (int i = 0; i < N; ++i){
		S += calcSafeCost(x_seq.block<2, 1>((STATE_DIM * i), 0), BodyValues);
	}
	A = 0;
	for (int i = 0; i < N; ++i) {
		JointAngle2EndPosition(end_pos, x_seq.block<2, 1>((STATE_DIM * i), 0));
		EndEff.x = end_pos(0);
		EndEff.y = end_pos(1);
		EndEff.z = 0;
		A += calcArmComfortR(RShould, length_arm, BodyValues[0], EndEff);
	}

	//Dealing with the infinite value
	if(!_finite(B))
		B = DBL_MAX;

#ifdef DEBUG
	std::cout << "phi : " << phi << std::endl;
	std::cout << "L : " << L << std::endl;
#endif //DEBUG 

	//std::cout<<"Phi "<<Phi<<" vs Safe "<<S<<std::endl; 

	return (Phi + B + V + S);
}


//Calculate the cost function Phi on final state constraint
double MPC::calcObjectivePhi(const Eigen::Matrix<double, STATE_DIM, 1>& goal_x, const int N)
{
	Eigen::Vector2d end_pos, end_vel;
	Eigen::Matrix<double, STATE_DIM, 1> x_end, diff;
	double temp;
	/*
	if(N>=20) {
		S(0, 0) = 150.0;
		S(1, 1) = 150.0;
		S(2, 2) = 50.0;
		S(3, 3) = 50.0; }
	else {
		S(0, 0) = 300.0;
		S(1, 1) = 300.0;
		S(2, 2) = 20.0;
		S(3, 3) = 20.0;}
		*/
	JointAngle2EndPosition(end_pos, goal_x.block<2, 1>(0, 0));
	JointVel2EndVel(end_vel, goal_x.block<2, 1>(0, 0), goal_x.block<2, 1>(2, 0));
	x_end.block<2, 1>(0, 0) = end_pos;
	x_end.block<2, 1>(2, 0) = end_vel;

	diff = x_end - xf_end;                    // Difference between actual and goal configuration   
	temp = diff.transpose() * S * diff;
	return  temp;
}

//Calculate the cost function PhiDash on final state constraint. Velocity only. Not position
// Final velocity should be zero so we do not need to take the difference.
double MPC::calcObjectivePhiDash(const Eigen::Matrix<double, STATE_DIM, 1>& goal_x, const int N)
{
	Eigen::Vector2d end_pos, end_vel;
	Eigen::Matrix<double, STATE_DIM, 1> x_end;
	double temp;
	/*					// Where do we get the value of S?
	if(N>=20) {
	S(0, 0) = 150.0;
	S(1, 1) = 150.0;
	S(2, 2) = 50.0;
	S(3, 3) = 50.0; }
	else {
	S(0, 0) = 300.0;
	S(1, 1) = 300.0;
	S(2, 2) = 20.0;
	S(3, 3) = 20.0;}
	*/
	
	JointVel2EndVel(end_vel, goal_x.block<2, 1>(0, 0), goal_x.block<2, 1>(2, 0));
	x_end(0, 0) = 0; // We are not considering position
	x_end(1, 0) = 0;		
	x_end.block<2, 1>(2, 0) = end_vel;
  
	temp = x_end.transpose() * S * x_end;
	return  temp;
}


//Calculate cost function B to penalize configurations exceeding limits
double MPC::calcObjectiveB(const Eigen::Matrix<double, STATE_DIM, 1>& cur_x, const Eigen::Matrix<double, INPUT_DOF, 1>& cur_u)
{
	double c1, c2, c3, c4;
	Eigen::Matrix<double, STATE_DIM, 1> dx;

	/* Calculate the Barrier function */
	calcStateFunc(dx, cur_x, cur_u);//Calculate the state equation
	//Consider limits of joints velocities and accelerations: if joint are in their limits no cost is applied
	c1 = ANVEL1_MAX - fabs(dx(0));
	if (c1 >= 0)
		c1 = 0;
	c2 = ANVEL2_MAX - fabs(dx(1));
	if (c2 >= 0)
		c2 = 0;
	//Considering the limits of angler accelaration
	c3 = ANACC1_MAX - fabs(dx(2));
	if (c3 >= 0)
		c3 = 0;
	c4 = ANACC2_MAX - fabs(dx(3));
	if (c4 >= 0)
		c4 = 0;

	
	return  ( R_C1 * (c1 * c1) / 2.0 + R_C2 * (c2 * c2) / 2.0
		+ R_C3 * (c3 * c3) / 2.0 + R_C4 * (c4 * c4) / 2.0);
	
}


//Calculate the Hamiltonian function
double MPC::calcHamiltonian(const Eigen::Matrix<double, STATE_DIM, 1>& cur_x, const Eigen::Matrix<double, INPUT_DOF, 1>& cur_u, const Eigen::Matrix<double, STATE_DIM, 1>& cur_lmd, pos3d_t LShould, pos3d_t RShould, pos3d_t Body, double length_arm, double* BodyValues, const int N)
{
	double objB, objV, objS, objA;
	Eigen::Matrix<double, 1, 1> objF;
	Eigen::Matrix<double, STATE_DIM, 1> f;
	Eigen::Matrix<double, STATE_DIM/2, 1> end_pos;
	pos3d_t EndEff;

	//Calculate the cost function B
	objB = calcObjectiveB(cur_x, cur_u);
#ifdef DEBUG
	std::cout << objB << std::endl;
#endif //DEBUG

	//Calculate Visibility Cost
	objV = calcVisibCost(cur_x.block<2, 1>(0, 0), BodyValues); 
    //Calculate Safety Cost
	objS = calcSafeCost(cur_x.block<2, 1>(0, 0), BodyValues);
	//Calculate Arm Comfort Cost
	JointAngle2EndPosition(end_pos, cur_x.block<2, 1>(0, 0));
	EndEff.x = end_pos(0);
	EndEff.y = end_pos(1);
	EndEff.z = 0;
	objA = calcArmComfortR(RShould, length_arm, BodyValues[0], EndEff);

	//Calculate the State equation
	calcStateFunc(f, cur_x, cur_u);
	objF = cur_lmd.transpose() * f;

#ifdef DEBUG
	std::cout << "ObjF is" << objF << std::endl;
	std::cout << "objF(0,0) is: "<<objF(0, 0) << std::endl;
#endif //DEBUG
	
	return (objS + objV + objA + objB + objF(0, 0));                  // H(x,u,lmd) = B(x) + lmd^T * f(x,u)
}

//Partial differentiate the Hamiltonian function by u
void MPC::calcdHdu(Eigen::Matrix<double, 1, INPUT_DOF>& dHdu, const Eigen::Matrix<double, STATE_DIM, 1>& cur_x,  const Eigen::Matrix<double, STATE_DIM, 1>& cur_lmd)
{
	Eigen::Matrix<double, STATE_DIM, INPUT_DOF> tempF;
	Eigen::Matrix<double, STATE_DIM, 1> f;

	//Partial differential of f w.r.t u is [0 M]'
	calcInertiaMat(cur_x.block<2, 1>(0, 0));  
	tempF.block<2, 2>(0, 0).setZero();
	tempF.block<2, 2>(2, 0) = M.inverse();

	dHdu = cur_lmd.transpose() * tempF;         // dH/du = -lmd^T * df/du

}


double MPC::calcVisibCost(const Eigen::Vector2d& joint_ang, double* BodyValues){

	Eigen::Vector2d end_pos; 
	JointAngle2EndPosition(end_pos, joint_ang);
	short k = 3; 
	
	pos2d_t Diff; 	
	Diff.x = (end_pos(0) - BodyValues[1])*cos(BodyValues[0]) + (end_pos(1) - BodyValues[2])*sin(BodyValues[0]); 
	Diff.y = -(end_pos(0) - BodyValues[1])*sin(BodyValues[0]) + (end_pos(1) - BodyValues[2])*cos(BodyValues[0]); 

	double RBangle = (atan2((Diff.y),(Diff.x)));                 // Orientation of vector from body to EndEffector wrt +x axis

	double VisCost = k*((RBangle)*(RBangle)); 
	
	return VisCost; 
}


double MPC::calcSafeCost(const Eigen::Vector2d& joint_ang, double* BodyValues){

	Eigen::Vector2d EE_pos, C_point1, C_point2; 
	JointAngle2EndPosition(EE_pos, joint_ang);

	// Safety cost obtained by inverse distance between End Effector and Body center
	int k1 = 3; 
	double k2 = 0.5; 
	double d_max = 0.7; 
	double SafeCost, SafeCost2; 

	double dist = sqrt(pow((EE_pos(0) - BodyValues[1]),2) + pow((EE_pos(1) - BodyValues[2]),2)); 
	C_point1(0) = L1*cos(joint_ang(0)) + 0.5*L2*cos(joint_ang(0)+joint_ang(1));
	C_point1(1) = L1*sin(joint_ang(0)) + 0.5*L2*sin(joint_ang(0)+joint_ang(1));
	C_point2(0) = L1*cos(joint_ang(0)) + 0.75*L2*cos(joint_ang(0)+joint_ang(1));
	C_point2(1) = L1*sin(joint_ang(0)) + 0.75*L2*sin(joint_ang(0)+joint_ang(1));
	double dist2 = sqrt(pow((C_point1(0) - BodyValues[1]),2) + pow((C_point1(1) - BodyValues[2]),2)); 
	double dist3 = sqrt(pow((C_point2(0) - BodyValues[1]),2) + pow((C_point2(1) - BodyValues[2]),2)); 

	if(dist < 0.20) {
		   SafeCost = 1000; }
	if(dist > d_max) 
		SafeCost = 0; 
	if(dist2 > d_max || dist3 > d_max)
		SafeCost2 = 0; 
	if(dist2 < 0.20 || dist3 < 0.20 ) { 
		//std::cout<<"Ouch!"<<std::endl; 
		SafeCost2 = 1000; }
	else {
		   SafeCost = k2*pow((1/dist - 1/d_max),k1)/k1;   
	       SafeCost2 = 0.5*k2*pow((1/dist2 - 1/d_max),k1)/k1 + 0.5*k2*pow((1/dist3 - 1/d_max),k1)/k1;}

	return SafeCost + SafeCost2; 
	//return 0; 
}

double MPC::calcArmComfortR(pos3d_t ShouldR, double Arm, double Bangle, pos3d_t EndEff){

	pos4d_t Q_r, dist_r;
	double r_q = pow(Arm, 2) - pow((ShouldR.z - EndEff.z), 2);
	double r_m = 0.1*0.1;
	double r_M = 0.7*0.7;
	double cost_1, cost_2, cost_3;

	// To compute the semicircle of reachability EndEffector Position must be rototraslated in ShoulderRF  (EE_S = R*(EE_W) - R*(Should_W)) on EE plane

	if (Bangle < 0) { Bangle = -(Bangle + 2 * M_PI - M_PI / 2); }
	else { Bangle = -(Bangle - M_PI / 2); }// Angle of rotation around z axis negative (CCW but inverse)

	// double Rot[4][4] = {// From RealWorld to ShoulderRF with negative angle
	//	{cos(Bangle), -sin(Bangle), 0, ShouldPos.x},
	//	{sin(Bangle), cos(Bangle), 0, ShouldPos.y},
	//	{0, 0, 0, EndEff.z},
	//	{0, 0, 0, 1} };

	pos3d_t EE_RS_t, EE_RS;
	EE_RS_t.x = EndEff.x - ShouldR.x;
	EE_RS_t.y = EndEff.y - ShouldR.y;

	EE_RS.x = EE_RS_t.x*cos(Bangle) - EE_RS_t.y*sin(Bangle);
	EE_RS.y = EE_RS_t.x*sin(Bangle) + EE_RS_t.y*cos(Bangle);
	EE_RS.z = EndEff.z - ShouldR.z;

	// Define joints range: range_1 = [-90 100]*deg  range_2 = [0,90]*deg range_3 = [180 -90]*deg range_4 = [0 135]*deg;
	Min_R.x = -M_PI / 2;
	Min_R.y = 0;
	Min_R.z = -M_PI / 2;
	Min_R.w = 0;
	Max_R.x = 100 * DEG2RAD;
	Max_R.y = M_PI / 2;
	Max_R.z = M_PI;
	Max_R.w = 3 * M_PI / 4;

	if ((EE_RS.y > 0.1) && (EE_RS.y*EE_RS.y + EE_RS.x*EE_RS.x <= r_q) && (EE_RS.y*EE_RS.y + EE_RS.x*EE_RS.x >= r_m) && (EE_RS.x >= -Arm / 2)) {             // Reachability set

		double phi = -M_PI / 3;                                                   // Choose 3 fixed swivel angles

		Q_r = InvKineArmR(EE_RS, ShouldR, Arm, phi);
		dist_r = JointDispR(Q_r, Max_R, Min_R);
		cost_1 = dist_r.x + dist_r.y + dist_r.z + dist_r.w;

		phi = -M_PI / 6;
		Q_r = InvKineArmR(EE_RS, ShouldR, Arm, phi);
		dist_r = JointDispR(Q_r, Max_R, Min_R);
		cost_2 = dist_r.x + dist_r.y + dist_r.z + dist_r.w;

		return std::min(cost_1, cost_2);

	}
	else if ((EE_RS.y > 0) && (EE_RS.y*EE_RS.y + EE_RS.x*EE_RS.x <= r_M)) {
		return 5.0;
	}
	else {
		return 10.0;
	}

}


pos4d_t MPC::InvKineArmR(pos3d_t EndEff_SR, pos3d_t ShouldR, double Arm, double phi) {

	pos3d_t n_r, u_r, v_r;
	pos3d_t Pc_r, Pe_r, Pw_r;
	double l1 = Arm / 2;
	double l2 = Arm / 2 + 0.05;
	double l3 = sqrt(EndEff_SR.x*EndEff_SR.x + EndEff_SR.y*EndEff_SR.y + EndEff_SR.z*EndEff_SR.z);
	Pw_r = EndEff_SR;
	//std::cout<<"Right RF "<<Pw_r.x<<","<<Pw_r.y<<","<<Pw_r.z<<std::endl;

	// unit vector to the EE position in the Shoulder RF
	n_r.x = EndEff_SR.x / l3;
	n_r.y = EndEff_SR.y / l3;
	n_r.z = EndEff_SR.z / l3;

	double cos_alfa = (l3*l3 + l1*l1 - l2*l2) / (2 * l3*l1);

	Pc_r.x = l1*cos_alfa*n_r.x;
	Pc_r.y = l1*cos_alfa*n_r.y;
	Pc_r.z = l1*cos_alfa*n_r.z;

	double R = l1*sin(acos(cos_alfa));

	// unit vector projection of -z axis (Shoulder RF) in the plane orthogonal to n
	u_r.x = n_r.x*n_r.z;
	u_r.y = n_r.y*n_r.z;
	u_r.z = n_r.z*n_r.z - 1;
	double mag = sqrt(u_r.x*u_r.x + u_r.y*u_r.y + u_r.z*u_r.z);
	u_r.x = u_r.x / mag;
	u_r.y = u_r.y / mag;
	u_r.z = u_r.z / mag;

	// unit vector given by the cross product: (n x v)
	v_r.x = n_r.y*u_r.z - n_r.z*u_r.y;
	v_r.y = n_r.z*u_r.x - n_r.x*u_r.z;
	v_r.z = n_r.x*u_r.y - n_r.y*u_r.x;

	// Impose position of the Elbow through Phi angle
	Pe_r.x = Pc_r.x + R*cos(phi)*u_r.x + R*sin(phi)*v_r.x;
	Pe_r.y = Pc_r.y + R*cos(phi)*u_r.y + R*sin(phi)*v_r.y;
	Pe_r.z = Pc_r.z + R*cos(phi)*u_r.z + R*sin(phi)*v_r.z;
	//std::cout<<"Pe:"<<Pe_r.x<<","<<Pe_r.y<<","<<Pe_r.z<<std::endl;

	// Find analitically Joint angles
	double teta_1, teta_2, teta_3, teta_4;
	double C, S;

	teta_2 = asin(-Pe_r.z / l1);
	teta_1 = atan2(Pe_r.y / cos(teta_2), Pe_r.x / cos(teta_2));

	C = (cos(teta_2)*(Pw_r.x*cos(teta_1) + Pw_r.y*sin(teta_1)) - Pw_r.z*sin(teta_2) - l1) / l2;
	S = sqrt(1 - C*C);
	teta_4 = atan2(S, C);

	S = (sin(teta_2)*(Pw_r.x*cos(teta_1) + Pw_r.y*sin(teta_1)) + Pw_r.z*cos(teta_2)) / (l2*sin(teta_4));
	C = (Pw_r.x*sin(teta_1) - Pw_r.y*cos(teta_1)) / (-L2*sin(teta_4));
	teta_3 = atan2(S, C);

	// Define Joint ranges to check feasibility

	pos4d_t Q_r;
	Q_r.x = teta_1;
	Q_r.y = teta_2;
	Q_r.z = teta_3;
	Q_r.w = teta_4;

	if ((Q_r.x > Max_R.x) || (Q_r.x < Min_R.x) || (Q_r.y < Min_R.y) || (Q_r.y > Max_R.y) || (Q_r.z < Min_R.z)
		|| (Q_r.z > Max_R.z) || (Q_r.w < Min_R.w) || (Q_r.w > Max_R.w)) {

		//std::cout<<"Exceeding Joint Limits!"<<std::endl;
		pos4d_t Q_err;
		Q_err.x = Q_err.y = Q_err.z = Q_err.w = 2 * M_PI;
		return Q_err;
	}


	return Q_r;

}

pos4d_t MPC::JointDispR(pos4d_t Joint_value, pos4d_t Max, pos4d_t Min) {

	// Define comfortable position as the mean value between joint limits and resting posture Q = [0, PI/2, PI/2, 0]
	pos4d_t Diff_r, Q_rest_r, Cost_r;

	double mean = (Min.x + Max.x) / 2;
	double range = (Max.x - Min.x) / 2;
	Diff_r.x = pow((mean - Joint_value.x) / range, 2);
	mean = (Min.y + Max.y) / 2;
	range = (Max.y - Min.y) / 2;
	Diff_r.y = pow((mean - Joint_value.y) / range, 2);
	mean = (Min.z + Max.z) / 2;
	range = (Max.z - Min.z) / 2;
	Diff_r.z = pow((mean - Joint_value.z) / range, 2);
	mean = (Min.w + Max.w) / 2;
	range = (Max.w - Min.w) / 2;
	Diff_r.w = pow((mean - Joint_value.w) / range, 2);

	//std::cout<<"Distance Mean "<<Diff_r.x<<" "<<Diff_r.y<<" "<<Diff_r.z<<" "<<Diff_r.w<<std::endl;

	Q_rest_r.x = 0;
	Q_rest_r.y = M_PI / 2;
	Q_rest_r.z = 0;
	Q_rest_r.w = M_PI / 6;

	range = (Max.x - Min.x);
	Cost_r.x = pow((Q_rest_r.x - Joint_value.x) / range, 2);
	range = (Max.y - Min.y);
	Cost_r.y = pow((Q_rest_r.y - Joint_value.y) / range, 2);
	range = (Max.z - Min.z);
	Cost_r.z = pow((Q_rest_r.z - Joint_value.z) / range, 2);
	range = (Max.w - Min.w);
	Cost_r.w = pow((Q_rest_r.w - Joint_value.w) / range, 2);

	//std::cout<<"Distance Rest "<<Cost_r.x<<" "<<Cost_r.y<<" "<<Cost_r.z<<" "<<Cost_r.w<<std::endl;

	int k1, k2, k3, k4;
	k1 = 2;
	k2 = 1;
	k3 = 0.8;
	k4 = 0.5;

	Cost_r.x = k2*Cost_r.x + k2*Diff_r.x;
	Cost_r.y = k1*Cost_r.y + k4*Diff_r.y;
	Cost_r.z = k3*Cost_r.z + k4*Diff_r.z;
	Cost_r.w = k4*Cost_r.w + k3*Diff_r.w;

	return Cost_r;

}


double MPC::calcSmoothCost(const Eigen::Matrix<double, STATE_DIM, 1>& cur_x,const Eigen::Matrix<double, INPUT_DOF, 1>& cur_u, const int N) {

	
	Eigen::Matrix<double, STATE_DIM, 1> x_end;
	Eigen::Vector2d diff;
	double SmthCost; 
	calcStateFunc(x_end, cur_x, cur_u);//Calculate the state equation

	diff = (x_end.block<2, 1>(2, 0) - x_end.block<2, 1>(0, 0))/(SAMPLING_TIME*0.001) ;

	if (N <=20) {
		SmthCost = 0.001*(diff(0)*diff(0) + diff(1)*diff(1)); }
	else {	SmthCost = 0.005*(diff(0)*diff(0) + diff(1)*diff(1)); }
	
	return 0; 
}



////Partial differentiate the f by u
//void MPC::calcdfdu(Eigen::Matrix<double, STATE_DIM, INPUT_DOF>& dfdu, const Eigen::Matrix<double, STATE_DIM, 1>& cur_x)
//{
//	calcInertiaMat(cur_x.block<2, 1>(0, 0));  
//	
//	dfdu.block<2, 2>(0, 0).setZero();
//	dfdu.block<2, 2>(2, 0) = M.inverse();
//}