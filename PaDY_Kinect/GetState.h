#pragma once

#include <Kinect.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>


//Euclide coordinate 2d
#ifndef _POS2D
#define _POS2D
typedef struct pos2d{
	double x;
	double y;
}pos2d_t;
#endif

//Euclide coordinate 3d
#ifndef _POS3D
#define _POS3D
typedef struct pos3d{
	double x;
	double y;
	double z;
}pos3d_t;
#endif

//Homogeneous coordinate (for quaternion)
#ifndef _POS4D
#define _POS4D
typedef struct pos4d{
	double x;
	double y;
	double z;
	double w;
}pos4d_t;
#endif


//Smart pointer
template<typename T>
class ComPtr{
private:

	T* ptr;

public:

	~ComPtr(){
		if (ptr != nullptr){
			ptr->Release();
			ptr = nullptr;
		}
	}
	T** operator & (){ return &ptr; }
	T* operator -> (){ return ptr; }
	operator T* () { return ptr; }
};

class GetState
{
public:

	GetState();// Constructor
	~GetState();// Destructor

	static const int SEQ_LENGTH = 919; //Length of sample data

	Joint joints[JointType::JointType_Count];//Joint positions
	JointOrientation jointsOri[JointType::JointType_Count];//Joint orientations

	int CTrackingState;//Tracking state of body
	int FTrackingState;//Tracking state of face

	int Initialize();// Initialize Kinect
	void UpdateColorFrame();  //Update ColorFrame    
	void UpdateBodyFrame();   //Update BodyFrame
	
	pos3d_t TransQuaterToEuler(Vector4 Quater);//Transform Quaternion to Euler angle
	pos3d_t TransKinectToRealWorld3D(CameraSpacePoint Pos);//Transform Kinect cordinate to Real world coordinate

	//pos2d_t BodyToScreen(CameraSpacePoint bodyPoint, int width, int height);  // Trabsforms kinect coordinates to a Color image point

	double Distance (pos3d_t joint_1, pos3d_t joint_2); 
	double ArmLength (pos3d_t Shoulder, pos3d_t Elbow, pos3d_t Hand); 
	cv::Mat Draw(pos2d_t Hand, pos2d_t Elbow, pos2d_t Shoulder, cv::Mat image); 

	void UpSideVision (pos3d_t SpineBase, pos3d_t ShouldLeft, pos3d_t ShouldRight, pos3d_t EE,  double EE_noF_x, double EE_noF_y, pos3d_t EndPoint, pos2d_t FirstLink);

	double VisibilityCost(pos3d_t EndEff, double BodyOrient, pos3d_t BodyPosition); 
	double ArmComfortR(pos3d_t ShouldPos, double Arm, double BodyOrient, pos3d_t EndEff); 
	double ArmComfortL(pos3d_t ShouldPos, double Arm, double BodyOrient, pos3d_t EndEff);;
	double SafetyCost(pos3d_t BodyPos, pos3d_t EndEff); 

	double TotalCost(double VisibCost, double cost_left, double cost_right, double SafeCost);
    
	pos4d_t InvKineArmR(pos3d_t EndEff_S, pos3d_t Should, double Arm, double phi); 
	pos4d_t InvKineArmL(pos3d_t EndEff_S, pos3d_t Should, double Arm, double phi);

	pos4d_t JointDispR (pos4d_t Joint_value, pos4d_t Max, pos4d_t Min); 
	pos4d_t JointDispL (pos4d_t Joint_value, pos4d_t Max, pos4d_t Min);

	pos4d_t TRRT(pos3d_t LShould, pos3d_t RShould, pos3d_t Body, double length_arm, double* Bangle);
	int TransitionTest(double Cost_i, double Cost_j, double delta, int T);
	void SwapShoulder(pos3d_t LShould, pos3d_t RShould, double* Left, double* Right);

	void ReadSampleData();//add 20200705

	pos3d_t dataBody[SEQ_LENGTH];
	pos3d_t dataLShoulder[SEQ_LENGTH];
	pos3d_t dataRShoulder[SEQ_LENGTH];
	pos3d_t dataRElbow[SEQ_LENGTH];
	pos3d_t dataRHand[SEQ_LENGTH];

private:

	IKinectSensor* Kinect;

	ICoordinateMapper* mapper;
	
	IColorFrameReader* ColorFrameReader;
	
	IBodyFrameReader* BodyFrameReader;
	IBody* Bodies[BODY_COUNT];

	cv::Scalar colors[6];
	
};

double Min(double a[][4], int n); 