#include "GetState.h"
#include <direct.h>
#include <math.h>
#include <fstream>
#include "opencv2\core\core.hpp"

const double PI = 3.141592653589;
const double DEG = PI/180;
pos4d_t Max_R, Min_R, Min_L, Max_L;

double GetState::VisibilityCost(pos3d_t EndEff, double Bangle, pos3d_t BodyPos) {

	short k = 3;
	// Find the angle between the EndEffector position and the gaze line

	pos2d_t Diff_R;
	Diff_R.x = (EndEff.x - BodyPos.x)*cos(Bangle) + (EndEff.y - BodyPos.y)*sin(Bangle);
	Diff_R.y = -(EndEff.x - BodyPos.x)*sin(Bangle) + (EndEff.y - BodyPos.y)*cos(Bangle);

	double RBangle = (atan2((Diff_R.y),(Diff_R.x)));                 // Orientation of vector from body to EndEffector wrt +x axis

	double VisCost = k*((RBangle)*(RBangle));

	return VisCost ;
}

double GetState::SafetyCost(pos3d_t BodyPos, pos3d_t EndEff) {

	// Safety cost obtained by inverse distance between End Effector and Body center
	int k1 = 3;
	double k2 = 0.5;
    int d_max = 3;
	double SafeCost;
	double dist = sqrt(pow((EndEff.x - BodyPos.x),2) + pow((EndEff.y - BodyPos.y),2));
	if(dist < 0.07) {
		   SafeCost = 50; }
	else {
		   SafeCost = k2*pow((1/dist - 1/d_max),k1)/k1;   }

	return SafeCost;
}

double GetState::TotalCost(double VisibCost, double cost_L, double cost_R, double SafeCost){

	short P_hand = 1;                                                         //Left or Right handiness penality

	cost_L = P_hand + cost_L;
	double ArmCost = std::min(cost_L, cost_R);

	//std::cout<<"Visibility: "<<VisibCost<<std::endl;
	//std::cout<<"Safety: "<<SafeCost<<std::endl;
	//std::cout<<"Arm: "<<ArmCost<<std::endl;

	return SafeCost + VisibCost + ArmCost;
}

double GetState::ArmComfortR(pos3d_t ShouldR, double Arm, double Bangle,pos3d_t EndEff) {

	pos4d_t Q_r, dist_r;
	double r_q = pow(Arm,2) - pow((ShouldR.z - EndEff.z),2);
	double r_m = 0.1*0.1;
	double r_M = 0.7*0.7;
	double cost_1, cost_2, cost_3;

    // To compute the semicircle of reachability EndEffector Position must be rototraslated in ShoulderRF  (EE_S = R*(EE_W) - R*(Should_W)) on EE plane

	if(Bangle < 0) { Bangle = -(Bangle + 2*PI - PI/2); }
	else { Bangle = -(Bangle - PI/2); }                                                         // Angle of rotation around z axis negative (CCW but inverse)

	// double Rot[4][4] = {                                                                     // From RealWorld to ShoulderRF with negative angle
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
	Min_R.x = -PI/2;
	Min_R.y = 0;
	Min_R.z = -PI/2;
	Min_R.w = 0;
	Max_R.x = 100*DEG;
	Max_R.y = PI/2;
	Max_R.z = PI;
	Max_R.w = 3*PI/4;

	if((EE_RS.y > 0.1) && ( EE_RS.y*EE_RS.y + EE_RS.x*EE_RS.x <= r_q ) && ( EE_RS.y*EE_RS.y + EE_RS.x*EE_RS.x >= r_m ) && (EE_RS.x >= -Arm/2)) {             // Reachability set

		 double phi = -PI/3;                                                   // Choose 3 fixed swivel angles

		 Q_r = InvKineArmR(EE_RS,ShouldR,Arm, phi);
		 dist_r = JointDispR(Q_r, Max_R, Min_R);
		 cost_1 = dist_r.x + dist_r.y + dist_r.z + dist_r.w;

		 phi = -PI/6;
		 Q_r = InvKineArmR(EE_RS,ShouldR,Arm, phi);
		 dist_r = JointDispR(Q_r, Max_R, Min_R);
		 cost_2 = dist_r.x + dist_r.y + dist_r.z + dist_r.w;

		 /*phi = 0;
		 Q_r = InvKineArmR(EE_RS,ShouldR,Arm, phi);
		 dist_r = JointDispR(Q_r, Max_R, Min_R);
		 cost_3 = dist_r.x + dist_r.y + dist_r.z + dist_r.w;
		  */

		 //cost_1 = std::min(cost_1, cost_3);

		 return std::min(cost_1, cost_2);

	}
	else if ((EE_RS.y > 0) && ( EE_RS.y*EE_RS.y + EE_RS.x*EE_RS.x <= r_M )) {
		return 5.0; }
	else {
		 return 10.0;
	}

	 }


pos4d_t GetState::InvKineArmR(pos3d_t EndEff_SR, pos3d_t ShouldR, double Arm, double phi) {

	pos3d_t n_r, u_r, v_r;
	pos3d_t Pc_r, Pe_r, Pw_r;
	double L1 = Arm/2;
	double L2 = Arm/2 + 0.05;
	double L3 = sqrt(EndEff_SR.x*EndEff_SR.x + EndEff_SR.y*EndEff_SR.y + EndEff_SR.z*EndEff_SR.z);
	Pw_r = EndEff_SR;
	//std::cout<<"Right RF "<<Pw_r.x<<","<<Pw_r.y<<","<<Pw_r.z<<std::endl;

	// unit vector to the EE position in the Shoulder RF
	n_r.x = EndEff_SR.x/L3;
	n_r.y = EndEff_SR.y/L3;
	n_r.z = EndEff_SR.z/L3;

	double cos_alfa = (L3*L3 + L1*L1 - L2*L2)/(2*L3*L1);

	Pc_r.x = L1*cos_alfa*n_r.x;
	Pc_r.y = L1*cos_alfa*n_r.y;
	Pc_r.z = L1*cos_alfa*n_r.z;

	double R = L1*sin(acos(cos_alfa));

	// unit vector projection of -z axis (Shoulder RF) in the plane orthogonal to n
	u_r.x = n_r.x*n_r.z;
	u_r.y = n_r.y*n_r.z;
	u_r.z = n_r.z*n_r.z -1;
	double mag = sqrt(u_r.x*u_r.x + u_r.y*u_r.y + u_r.z*u_r.z);
	u_r.x = u_r.x/mag;
	u_r.y = u_r.y/mag;
	u_r.z = u_r.z/mag;

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

	teta_2 = asin(-Pe_r.z/L1);
	teta_1 = atan2(Pe_r.y/cos(teta_2), Pe_r.x/cos(teta_2));

	C = (cos(teta_2)*(Pw_r.x*cos(teta_1) + Pw_r.y*sin(teta_1)) -Pw_r.z*sin(teta_2) - L1)/L2;
	S = sqrt( 1 - C*C);
	teta_4 = atan2(S,C);

	S = (sin(teta_2)*(Pw_r.x*cos(teta_1) + Pw_r.y*sin(teta_1)) +Pw_r.z*cos(teta_2))/(L2*sin(teta_4));
	C = (Pw_r.x*sin(teta_1) - Pw_r.y*cos(teta_1))/(-L2*sin(teta_4));
	teta_3 = atan2(S,C);

	// Define Joint ranges to check feasibility

	pos4d_t Q_r;
	Q_r.x = teta_1;
	Q_r.y = teta_2;
	Q_r.z = teta_3;
	Q_r.w = teta_4;

	if((Q_r.x > Max_R.x) || (Q_r.x < Min_R.x) || (Q_r.y < Min_R.y) || (Q_r.y > Max_R.y) || (Q_r.z < Min_R.z)
		||(Q_r.z > Max_R.z) || (Q_r.w < Min_R.w) || (Q_r.w > Max_R.w)) {

			//std::cout<<"Exceeding Joint Limits!"<<std::endl;
			pos4d_t Q_err;
	        Q_err.x = Q_err.y = Q_err.z = Q_err.w = 2*PI;
	        return Q_err; }


	return Q_r;

}


pos4d_t GetState::JointDispR (pos4d_t Joint_value, pos4d_t Max, pos4d_t Min) {

	// Define comfortable position as the mean value between joint limits and resting posture Q = [0, PI/2, PI/2, 0]
	pos4d_t Diff_r, Q_rest_r, Cost_r;

	double mean = (Min.x + Max.x)/2;
	double range = (Max.x - Min.x)/2;
	Diff_r.x = pow((mean - Joint_value.x)/range,2);
	mean = (Min.y + Max.y)/2;
	range = (Max.y - Min.y)/2;
	Diff_r.y = pow((mean - Joint_value.y)/range,2);
	mean = (Min.z + Max.z)/2;
	range = (Max.z - Min.z)/2;
	Diff_r.z = pow((mean - Joint_value.z)/range,2);
	mean = (Min.w + Max.w)/2;
	range = (Max.w - Min.w)/2;
	Diff_r.w = pow((mean - Joint_value.w)/range,2);

	//std::cout<<"Distance Mean "<<Diff_r.x<<" "<<Diff_r.y<<" "<<Diff_r.z<<" "<<Diff_r.w<<std::endl;

	Q_rest_r.x = 0;
	Q_rest_r.y = PI/2;
	Q_rest_r.z = 0;
	Q_rest_r.w = PI/6;

	range = (Max.x - Min.x);
	Cost_r.x = pow((Q_rest_r.x - Joint_value.x)/range,2);
	range = (Max.y - Min.y);
	Cost_r.y = pow((Q_rest_r.y - Joint_value.y)/range,2);
	range = (Max.z - Min.z);
	Cost_r.z = pow((Q_rest_r.z - Joint_value.z)/range,2);
	range = (Max.w - Min.w);
	Cost_r.w = pow((Q_rest_r.w - Joint_value.w)/range,2);

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


double GetState::ArmComfortL(pos3d_t ShouldPos, double Arm, double Bangle ,pos3d_t EndEff) {

	pos4d_t Q, dist_l;
	double r_q = pow(Arm,2) - pow((ShouldPos.z - EndEff.z),2);
	double r_m = 0.1*0.1;
	double r_M = 0.7*0.7;
	double cost_1, cost_2, cost_3;

    // To compute the semicircle of reachability EndEffector Position must be rototraslated in ShoulderRF  (EE_S = R*(EE_W) - R*(Should_W)) on EE plane

	if(Bangle < 0) { Bangle = -(Bangle + 2*PI + PI/2); }
	else { Bangle = -(Bangle + PI/2); }                                                         // Angle of rotation around z axis negative (CCW but inverse)

	/* double Rot[4][4] = {                                                                     // From RealWorld to ShoulderRF with negative angle
		{cos(Bangle), -sin(Bangle), 0, ShouldPos.x},
		{sin(Bangle), cos(Bangle), 0, ShouldPos.y},
		{0, 0, 0, EndEff.z},
		{0, 0, 0, 1} };                                */

	pos3d_t EE_S_t, EE_SL;
	EE_S_t.x = EndEff.x - ShouldPos.x;
	EE_S_t.y = EndEff.y - ShouldPos.y;

	EE_SL.x = EE_S_t.x*cos(Bangle) - EE_S_t.y*sin(Bangle);
	EE_SL.y = EE_S_t.x*sin(Bangle) + EE_S_t.y*cos(Bangle);
	EE_SL.z = EndEff.z - ShouldPos.z;

	// Define Joint ranges:  range_1 = [-100 90]*deg  range_2 = [0,90]*deg range_3 = [-180 90]*deg range_4 = [0 135]*deg;
	Min_L.x = -100*DEG;
	Min_L.y = 0;
	Min_L.z = -PI;
	Min_L.w = 0;
	Max_L.x = PI/2;
	Max_L.y = PI/2;
	Max_L.z = PI/2;
	Max_L.w = 3*PI/4;

	if((EE_SL.y < -0.1) && ( EE_SL.y*EE_SL.y + EE_SL.x*EE_SL.x <= r_q ) && ( EE_SL.y*EE_SL.y + EE_SL.x*EE_SL.x >= r_m ) && (EE_SL.x >= -Arm/2)) {             // Reachability set
		 double phi = PI/6;                                                                      // Choose 3 fixed swivel angles
		 Q = InvKineArmL(EE_SL,ShouldPos,Arm, phi);
		 dist_l = JointDispL(Q, Max_L, Min_L);
		 cost_1 = dist_l.x + dist_l.y + dist_l.z + dist_l.w;

		 phi = PI/3;
		 Q = InvKineArmL(EE_SL,ShouldPos,Arm, phi);
		 dist_l = JointDispL(Q, Max_L, Min_L);
		 cost_2 = dist_l.x + dist_l.y + dist_l.z + dist_l.w;

		 /*phi = 0;
		 Q = InvKineArmL(EE_SL,ShouldPos,Arm, phi);
		 dist_l = JointDispL(Q, Max_L, Min_L);
		 cost_3 = dist_l.x + dist_l.y + dist_l.z + dist_l.w;
		*/

		 //cost_1 = std::min(cost_1, cost_3);

		 return std::min(cost_1, cost_2);
	}
	else if((EE_SL.y < 0) && ( EE_SL.y*EE_SL.y + EE_SL.x*EE_SL.x <= r_M )){
		return 5.0;  }
	else{
		return 10.0;
	     }
}


pos4d_t GetState::InvKineArmL(pos3d_t EndEff_S, pos3d_t Should, double Arm, double phi) {

	pos3d_t n, u, v;
	pos3d_t Pc, Pe, Pw;
	double L1 = Arm/2;
	double L2 = Arm/2 + 0.05;
	double L3 = sqrt(EndEff_S.x*EndEff_S.x + EndEff_S.y*EndEff_S.y + EndEff_S.z*EndEff_S.z);
	Pw = EndEff_S;
	//std::cout<<"Left RF "<<Pw.x<<","<<Pw.y<<","<<Pw.z<<std::endl;

	// unit vector to the EE position in the Shoulder RF
	n.x = EndEff_S.x/L3;
	n.y = EndEff_S.y/L3;
	n.z = EndEff_S.z/L3;

	double cos_alfa = (L3*L3 + L1*L1 - L2*L2)/(2*L3*L1);

	Pc.x = L1*cos_alfa*n.x;
	Pc.y = L1*cos_alfa*n.y;
	Pc.z = L1*cos_alfa*n.z;

	double R = L1*sin(acos(cos_alfa));

	// unit vector projection of -z axis (Shoulder RF) in the plane orthogonal to n
	u.x = n.x*n.z;
	u.y = n.y*n.z;
	u.z = n.z*n.z -1;
	double mag = sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
	u.x = u.x/mag;
	u.y = u.y/mag;
	u.z = u.z/mag;

	// unit vector given by the cross product: (n x v)
	v.x = n.y*u.z - n.z*u.y;
	v.y = n.z*u.x - n.x*u.z;
	v.z = n.x*u.y - n.y*u.x;

	// Impose position of the Elbow through Phi angle
	Pe.x = Pc.x + R*cos(phi)*u.x + R*sin(phi)*v.x;
	Pe.y = Pc.y + R*cos(phi)*u.y + R*sin(phi)*v.y;
	Pe.z = Pc.z + R*cos(phi)*u.z + R*sin(phi)*v.z;
	//std::cout<<"Pe:"<<Pe.x<<","<<Pe.y<<","<<Pe.z<<std::endl;

    // Find analitically Joint angles
	double teta_1, teta_2, teta_3, teta_4;
	double C, S;

	teta_2 = asin(-Pe.z/L1);
	teta_1 = atan2(Pe.y/cos(teta_2), Pe.x/cos(teta_2));

	C = (cos(teta_2)*(Pw.x*cos(teta_1) + Pw.y*sin(teta_1)) -Pw.z*sin(teta_2) - L1)/L2;
	S = sqrt( 1 - C*C);
	teta_4 = atan2(S,C);

	S = (sin(teta_2)*(Pw.x*cos(teta_1) + Pw.y*sin(teta_1)) +Pw.z*cos(teta_2))/(L2*sin(teta_4));
	C = (Pw.x*sin(teta_1) - Pw.y*cos(teta_1))/(-L2*sin(teta_4));
	if(S<=0) {teta_3 = atan2(S,C) + PI;}
	else { teta_3 = atan2(S,C) - PI; }

	// Check feasibility

	pos4d_t Q;
	Q.x = teta_1;
	Q.y = teta_2;
	Q.z = teta_3;
	Q.w = teta_4;

	if((Q.x > Max_L.x) || (Q.x < Min_L.x) || (Q.y < Min_L.y) || (Q.y > Max_L.y) || (Q.z < Min_L.z)
		||(Q.z > Max_L.z) || (Q.w < Min_L.w) || (Q.w > Max_L.w)) {

			//std::cout<<"Exceeding Joint Limits!"<<std::endl;
			pos4d_t Q_err_l;
	        Q_err_l.x = Q_err_l.y = Q_err_l.z = Q_err_l.w = 2*PI;
	        return Q_err_l; }


	return Q;

}


pos4d_t GetState::JointDispL (pos4d_t Joint_value, pos4d_t Max, pos4d_t Min) {

	// Define comfortable position as the mean value between joint limits and resting posture Q = [0, PI/2, PI/2, 0]
	pos4d_t Diff_l, Q_rest_l, Cost_l;

	double mean = (Min.x + Max.x)/2;
	double range = (Max.x - Min.x)/2;
	Diff_l.x = pow((mean - Joint_value.x)/range,2);
	mean = (Min.y + Max.y)/2;
	range = (Max.y - Min.y)/2;
	Diff_l.y = pow((mean - Joint_value.y)/range,2);
	mean = (Min.z + Max.z)/2;
	range = (Max.z - Min.z)/2;
	Diff_l.z = pow((mean - Joint_value.z)/range,2);
	mean = (Min.w + Max.w)/2;
	range = (Max.w - Min.w)/2;
	Diff_l.w = pow((mean - Joint_value.w)/range,2);

	//std::cout<<"Distance Mean "<<Diff_l.x<<" "<<Diff_l.y<<" "<<Diff_l.z<<" "<<Diff_l.w<<std::endl;

	Q_rest_l.x = 0;
	Q_rest_l.y = PI/2;
	Q_rest_l.z = 0;
	Q_rest_l.w = PI/6;

	range = (Max.x - Min.x);
	Cost_l.x = pow((Q_rest_l.x - Joint_value.x)/range,2);
	range = (Max.y - Min.y);
	Cost_l.y = pow((Q_rest_l.y - Joint_value.y)/range,2);
	range = (Max.z - Min.z);
	Cost_l.z = pow((Q_rest_l.z - Joint_value.z)/range,2);
	range = (Max.w - Min.w);
	Cost_l.w = pow((Q_rest_l.w - Joint_value.w)/range,2);

	//std::cout<<"Distance Rest "<<Cost_l.x<<" "<<Cost_l.y<<" "<<Cost_l.z<<" "<<Cost_l.w<<std::endl;

	int k1, k2, k3, k4;
	k1 = 2;
	k2 = 1;
	k3 = 0.8;
	k4 = 0.5;

	Cost_l.x = k2*Cost_l.x + k2*Diff_l.x;
	Cost_l.y = k1*Cost_l.y + k4*Diff_l.y;
	Cost_l.z = k3*Cost_l.z + k4*Diff_l.z;
	Cost_l.w = k4*Cost_l.w + k3*Diff_l.w;

	return Cost_l;

}


pos4d_t GetState::TRRT(pos3d_t LShould, pos3d_t RShould, pos3d_t Body, double length_R, double* BodyValues) {

	//std::ofstream cost_file("C:\\Users\\Fujitsu\\Documents\\Cost.txt",std::ios::app);      // write on file appending

	// double length_R = 0.5;
	double delta;
	pos2d_t m;
	pos4d_t Result;

	// Suppose EndEffector plane is always at z = 0;
	int height = 0;
	int width = 200;
	int depth = 400;

	// Start searching from body position
	pos3d_t EndEff = Body;
	/*EndEff.x = 0;
	EndEff.y = 0; */
	EndEff.z = height;
	pos3d_t New, Rand;

	/*EndEff.x = (std::rand() % depth);
	EndEff.x = EndEff.x/100;
	EndEff.y = (std::rand() % depth +(-width));
	EndEff.y = EndEff.y/100; */

	// std::cout<<"Selected Start: "<<EndEff.x<<","<<EndEff.y<<std::endl;

	// Orientation of the body in real world
	pos2d_t orient;
	orient.y = (LShould.y - RShould.y);
	orient.x = (LShould.x - RShould.x);
	double Bangle = atan2(-(orient.x),(orient.y));
	BodyValues[0] = Bangle;
	BodyValues[1] = Body.x;
	BodyValues[2] = Body.y;

	double VisCost_Base, cost_L_Base, cost_R_Base, cost_Base, SafeCost_Base;

	short i, alfa = 2, nFail = 0;
	int T = 0.1;

	for(i=0; i<140; i++)  // Why 140?
	{

		if(i<=10) {delta = 0.07;}

		VisCost_Base = VisibilityCost(EndEff, Bangle, Body);
		cost_L_Base = ArmComfortL(LShould, length_R, Bangle, EndEff);
		cost_R_Base = ArmComfortR(RShould, length_R, Bangle, EndEff);
		SafeCost_Base = SafetyCost(Body, EndEff);
		cost_Base = TotalCost(VisCost_Base, cost_L_Base, cost_R_Base,SafeCost_Base);

		// Take a random position on the plane
		Rand.x = (std::rand() % depth);
		Rand.x = Rand.x/100;
		Rand.y = (std::rand() % depth +(-width));
		Rand.y = Rand.y/100;
		Rand.z = height;

		// Extract the new EndEffector position at distance delta in the direction of Rand
		m.x = (EndEff.y - Rand.y);
		m.y = (EndEff.x - Rand.x);
		New.x = EndEff.x + delta*m.x/sqrt(m.x*m.x + m.y*m.y);
		New.y = EndEff.y - delta*m.y/sqrt(m.x*m.x + m.y*m.y);
		New.z = height;

		double VisCost = VisibilityCost(New, Bangle, Body);
		double cost_L = ArmComfortL(LShould, length_R, Bangle, New);
		double cost_R = ArmComfortR(RShould, length_R, Bangle, New);
		double SafeCost = SafetyCost(Body, New);
		double cost_New = TotalCost(VisCost, cost_L, cost_R, SafeCost);
	
		int v = TransitionTest(cost_Base, cost_New, delta, T);

		if (v == 1) 
		{
			EndEff = New;
			T = 0.1;
		    nFail = 0;
		    delta = 0.03;
		}
		else 
		{
			nFail = nFail +1;
			if(nFail > 6 && cost_Base > 5) 
			{
				delta = delta*alfa;
				T = alfa*T;
				delta = 0.5;
			    /*nFail = 0;*/
			}
		}

		//cost_file<<cost_Base<<" "<<delta<<std::endl;
		//std::cout<<EndEff.x<<","<<EndEff.y<<","<<EndEff.z<<std::endl;

	}

	//cost_file<<std::endl;
	//std::cout<<"Cost "<<cost_Base<<std::endl;

	Result.x = EndEff.x;
	Result.y = EndEff.y;
	Result.z = EndEff.z;
	Result.w = cost_Base;

	return Result;

}


int GetState::TransitionTest(double cost_i, double cost_j, double delta, int T) {

	if (cost_j < cost_i) {
	    return 1;   }
	else {
		double p = exp( -(cost_j-cost_i)/(delta*100*T));

		if((std::rand() %100)/100 < p) {
			return 1; }
		return 0;}

}


