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
	calcInertiaMat(cur_x.block<2, 1>(0, 0)); //Inertia(current position)
	calcCoriolisTerm(cur_x.block<2, 1>(0, 0), cur_x.block<2, 1>(2, 0));  //Coriolis(current position, current velocity)

	//Calculate X_dot vector using system equations and lagrangian dynamics 
	f.block<2, 1>(0, 0) = cur_x.block<2, 1>(2, 0); 
	f.block<2, 1>(2, 0) = -M.inverse()*C + M.inverse()*cur_u;
}

//Calculate the input by gradient method
int MPC::calcInputbyGradientMethod(Eigen::Matrix<double, INPUT_DOF, 1>& input, const Eigen::Matrix<double, STATE_DIM, 1>& cur_x, double T, double dt, double* BodyValues)
{
	int STEP_N;//Step size
	const double DEV_R = 0.618;                             //Golden Ratio
	const double EXP_STEP = 0.05;                           //Step size for linear search
	double J, pre_J, temp_J, J1, J2;                        //Cost values
	double scale, beta, alpha_l, alpha_u, alpha1, alpha2;
	Eigen::VectorXd pre_dir, dir;                           //Descending direction of gradient
	Eigen::VectorXd x;                                      //Group of state vector
	Eigen::VectorXd u;                                      //Group of input vector
	Eigen::VectorXd temp_u;
	Eigen::VectorXd pre_F, F;                               //Gradient of Hamiltonian function
	Eigen::Matrix<double, STATE_DIM, 1> f;                  //State equation
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

	//Resize vectors
	pre_dir.resize(INPUT_DOF * STEP_N);
	dir.resize(INPUT_DOF * STEP_N);
	x.resize(STATE_DIM * (STEP_N + 1));
	u.resize(INPUT_DOF * STEP_N);
	temp_u.resize(INPUT_DOF * STEP_N);
	pre_F.resize(INPUT_DOF * STEP_N);
	F.resize(INPUT_DOF * STEP_N);

	//Initialization
	u.setConstant(0.0);
	pre_J = 10000000000;

	/* Start gradient method for calculating the input */
	for (int cycle = 0; cycle < CYCLE_NEWTON; ++cycle)
	{
		#ifdef DEBUG
				std::cout << "cycle : " << cycle << std::endl;
		#endif //DEBUG

		/* Calculate the gradient vector */
		calcGradient(F, x, cur_x, u, STEP_N, dt, BodyValues); 
		
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
		temp_J = calcObjective(x, u, STEP_N,BodyValues);
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
			J = calcObjective(x, temp_u, STEP_N,BodyValues);

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
					J1 = calcObjective(x, temp_u, STEP_N,BodyValues);// Calculate the objective function
					/* J2 */
					temp_u = u - alpha2 * dir;//Update the input u
					//Calculte the group of state vectors by state equation
					x.block<STATE_DIM, 1>(0, 0) = cur_x;//Set the initial state
					for (int n = 0; n < STEP_N; ++n){
						calcStateFunc(f, x.block<STATE_DIM, 1>((n * STATE_DIM), 0), temp_u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0));
						x.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0) = x.block<STATE_DIM, 1>((n * STATE_DIM), 0) + f * dt;
					}
					J2 = calcObjective(x, temp_u, STEP_N,BodyValues);// Calculate the objective function

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
		J = calcObjective(x, temp_u, STEP_N,BodyValues);   
		
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
void MPC::calcGradient(Eigen::VectorXd& F, Eigen::VectorXd& x, const Eigen::VectorXd& cur_x, const Eigen::VectorXd& cur_u, const int N, const double dt, double* BodyValues)
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
	temp = calcObjectivePhi(x.block<STATE_DIM, 1>((STATE_DIM * N), 0),N);
	for (int i = 0; i < STATE_DIM; ++i){
		temp_x = x.block<STATE_DIM, 1>((STATE_DIM * N), 0);             //Take estimated last state vector
		temp_x(i) += EPS;                                               //Incremental ratio 
		dPhi(i) = (calcObjectivePhi(temp_x,N) - temp) / EPS;
	}
	//Calculate costate vectors imposing first the optimal constraint: lambda(N)=dPhi/dx
	lmd.block<STATE_DIM, 1>((STATE_DIM * N), 0) = dPhi;
	for (int n = (N - 1); n > 0; --n){
		//Partial differentiate the Hamiltonian
		temp = calcHamiltonian(x.block<STATE_DIM, 1>((n * STATE_DIM), 0), cur_u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0), lmd.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0), BodyValues,N);
		for (int i = 0; i < STATE_DIM; ++i){
			temp_x = x.block<STATE_DIM, 1>((n * STATE_DIM), 0);
			temp_x(i) += EPS;
			dHdx(i) = (calcHamiltonian(temp_x, cur_u.block<INPUT_DOF, 1>((n * INPUT_DOF), 0), lmd.block<STATE_DIM, 1>(((n + 1) * STATE_DIM), 0), BodyValues,N) - temp) / EPS;
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
double MPC::calcObjective(const Eigen::VectorXd& x_seq, const Eigen::VectorXd& u_seq, const int N, double* BodyValues)
{
	double Phi, B, V, S;

	//Calculate the cost function Phi on x(N)
	Phi = calcObjectivePhi(x_seq.block<STATE_DIM, 1>((STATE_DIM * N), 0),N);

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

	JointAngle2EndPosition(end_pos, goal_x.block<2, 1>(0, 0));
	JointVel2EndVel(end_vel, goal_x.block<2, 1>(0, 0), goal_x.block<2, 1>(2, 0));
	x_end.block<2, 1>(0, 0) = end_pos;
	x_end.block<2, 1>(2, 0) = end_vel;

	diff = x_end - xf_end;                    // Difference between actual and goal configuration   
	temp = diff.transpose() * S * diff;
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
double MPC::calcHamiltonian(const Eigen::Matrix<double, STATE_DIM, 1>& cur_x, const Eigen::Matrix<double, INPUT_DOF, 1>& cur_u, const Eigen::Matrix<double, STATE_DIM, 1>& cur_lmd, double* BodyValues, const int N)
{
	double objB, objV, objS;
	Eigen::Matrix<double, 1, 1> objF;
	Eigen::Matrix<double, STATE_DIM, 1> f;

	//Calculate the cost function B
	objB = calcObjectiveB(cur_x, cur_u);
#ifdef DEBUG
	std::cout << objB << std::endl;
#endif //DEBUG

	//Calculate Visibility Cost
	objV = calcVisibCost(cur_x.block<2, 1>(0, 0), BodyValues); 
    //Calculate Safety Cost
	objS = calcSafeCost(cur_x.block<2, 1>(0, 0), BodyValues);

	//Calculate the State equation
	calcStateFunc(f, cur_x, cur_u);
	objF = cur_lmd.transpose() * f;
#ifdef DEBUG
	std::cout << objF(0, 0) << std::endl;
#endif //DEBUG
	
	return (objS + objV + objB + objF(0, 0));                  // H(x,u,lmd) = B(x) + lmd^T * f(x,u)
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
	double k2 = 8; 
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