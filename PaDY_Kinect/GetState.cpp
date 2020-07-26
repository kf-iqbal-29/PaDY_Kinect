#include "GetState.h"
#include <direct.h> 
#include "stdafx.h"
#include<iostream>
#include<sstream>
#include<fstream>
#include<string>
#include <Kinect.Face.h>
#include<math.h>

BOOLEAN BodyisTracked; 


// Constructor
GetState::GetState()
{
	// Set colors
	colors[0] = cv::Scalar(255, 0, 0);
	colors[1] = cv::Scalar(0, 255, 0);
	colors[2] = cv::Scalar(0, 0, 255);
	colors[3] = cv::Scalar(255, 255, 0);
	colors[4] = cv::Scalar(255, 0, 255);
	colors[5] = cv::Scalar(0, 255, 255);

}

// Destructor
GetState::~GetState()
{
	/*
	//Release all objects
	for (int i = 0; i < BODY_COUNT; i++){
		Bodies[i]->Release();
	}

	BodyFrameReader->Release();
	ColorFrameReader->Release();

	mapper->Release();

	//Close Kinect
	if (Kinect){
		Kinect->Close();
	}
	Kinect->Release();
	*/
}


//Initialize Kinect
int GetState::Initialize()
{
	HRESULT hr;
	ReadSampleData();
	
	//Get a default Kinect data
	hr = GetDefaultKinectSensor(&Kinect);
	if (FAILED(hr)){
		std::cout << "Kinect can not be connected" << std::endl;
		getchar();
		return -1;  }

	//Open Kinect
	Kinect->Open();

	BOOLEAN isOpen = false;
	Kinect->get_IsOpen(&isOpen);
	if (!isOpen){
		std::cout << "Kinect can not be opend" << std::endl;
		return -1;
	}

	//For cordinate transformation
	Kinect->get_CoordinateMapper(&mapper);

	//Get ColorReader
	ComPtr<IColorFrameSource> ColorFrameSource;
	Kinect->get_ColorFrameSource(&ColorFrameSource);
	ColorFrameSource->OpenReader(&ColorFrameReader);

	//Get BodyReader
	ComPtr<IBodyFrameSource> BodyFrameSource;
	Kinect->get_BodyFrameSource(&BodyFrameSource);
	hr = BodyFrameSource->OpenReader(&BodyFrameReader);

	//Check Kinect state
	if (!Kinect || FAILED(hr)){
		std::cout << "No ready Kinect found!" << std::endl;
		getchar();
		return -1;
	}
	
	return 1;
}

//Update BodyFrame
void GetState::UpdateBodyFrame()
{
	int first_trig = 1;

	CTrackingState = 0;

	if (!BodyFrameReader){
		std::cout << "Error : BodyFrameReader" << std::endl;
		return;
	}

	//Get BodyFrame
	ComPtr<IBodyFrame> BodyFrame;
	auto ret = BodyFrameReader->AcquireLatestFrame(&BodyFrame);

	if (ret == S_OK){
		//Update body data
		BodyFrame->GetAndRefreshBodyData(BODY_COUNT, Bodies);

		first_trig = 1;
		//Get body data
		for (auto body : Bodies)
		{
			if (!first_trig) break;
			if (body == nullptr) continue;

			BodyisTracked = false;
			body->get_IsTracked(&BodyisTracked);//Check the tracking state

			if (!BodyisTracked) continue;
			CTrackingState = 1;

			//Get the joint positions
			body->GetJoints(JointType::JointType_Count, joints);
			//Get the joint orientations
			body->GetJointOrientations(JointType::JointType_Count, jointsOri);

			first_trig = 0;		
		}
	
	}
}

pos3d_t GetState::TransKinectToRealWorld3D(CameraSpacePoint Pos)
{
	pos3d_t TransedPos;
	
	TransedPos.x = Pos.Z;
	TransedPos.y = Pos.X;
	TransedPos.z = Pos.Y;

	//TransedPos.x = -Pos.Z;
	//TransedPos.y = -Pos.X;
	//TransedPos.z = Pos.Y;

	return TransedPos;
}

pos3d_t GetState::TransLogToRealWorld3D(pos3d Pos)
{
	pos3d_t TransedPos;
	TransedPos.x = Pos.x;
	TransedPos.y = Pos.y;
	TransedPos.z = Pos.z;

	return TransedPos;
}


//Transform quaternion to Euler angle
pos3d_t GetState::TransQuaterToEuler(Vector4 Quater)
{
	const double PI = 3.141592653589;
	pos3d_t Euler;

	const double w2 = Quater.w*Quater.w;
	const double x2 = Quater.x*Quater.x;
	const double y2 = Quater.y*Quater.y;
	const double z2 = Quater.z*Quater.z;
	const double unitLength = w2 + x2 + y2 + z2; // Normalised == 1, otherwise correction divisor.
	const double abcd = Quater.w*Quater.x + Quater.y*Quater.z;
	const double eps = 1e-7; // TODO: pick from your math lib instead of hardcoding.

	if (abcd > (0.5 - eps)*unitLength){
		Euler.z = 2 * atan2(Quater.y, Quater.w);
		Euler.x = PI;
		Euler.y = 0;
	}
	else if (abcd < (-0.5 + eps)*unitLength){
		Euler.z = -2 * atan2(Quater.y, Quater.w);
		Euler.x = -PI;
		Euler.y = 0;
	}
	else{
		const double adbc = Quater.w*Quater.z - Quater.x*Quater.y;
		const double acbd = Quater.w*Quater.y - Quater.x*Quater.z;
		Euler.z = atan2(2 * adbc, 1 - 2 * (z2 + x2));
		Euler.x = asin(2 * abcd / unitLength);              
		Euler.y = atan2(2 * acbd, 1 - 2 * (y2 + x2));
	}

	return Euler;
}


double GetState::Distance(pos3d_t joint_1, pos3d_t joint_2) {

	pos3d_t diff; 

	diff.x = joint_1.x - joint_2.x;
	diff.y = joint_1.y - joint_2.y;
	diff.z = joint_1.z - joint_2.z;

	double temp = std::hypot( diff.x, diff.y);   // returns square root of (x^2 + y^2) 
	double distance = std::hypot(diff.z, temp); 

	return distance; 

}

double GetState::ArmLength (pos3d_t Shoulder, pos3d_t Elbow, pos3d_t Hand) {

	double temp = Distance(Shoulder, Elbow); 
	double temp2 = Distance(Elbow, Hand); 

	double length = temp + temp2; 

	return length; 

}

cv::Mat GetState::Draw(pos2d_t Hand, pos2d_t Elbow, pos2d_t Shoulder, cv::Mat image) {

	int rad = 7; 

	cv::line(image, cv::Point(Hand.x/2, Hand.y/2), cv::Point(Elbow.x/2, Elbow.y/2), cv::Scalar(0,0,255), 1,8,0);
	cv::line(image, cv::Point(Elbow.x/2, Elbow.y/2), cv::Point(Shoulder.x/2, Shoulder.y/2), cv::Scalar(0,0,255), 1,8,0);

	cv::circle(image, cv::Point(Hand.x/2, Hand.y/2), rad, cv::Scalar(0,0,255), 1,8,0);                     // In drawing position must be scaled...
	cv::circle(image, cv::Point(Elbow.x/2, Elbow.y/2), rad, cv::Scalar(0,0,255), 1,8,0); 
	cv::circle(image, cv::Point(Shoulder.x/2, Shoulder.y/2), rad, cv::Scalar(0,0,255), 1,8,0);

	return image; 

} 

void GetState::UpSideVision(pos3d_t SpineBase, pos3d_t ShouldLeft, pos3d_t ShouldRight, pos3d_t EE, double EE_noF_x, double EE_noF_y, pos3d_t EndPoint, pos2d_t FirstLink) {

	int width = 400,  cell_size = 10; 
	int body = 5, hand = 5;	

	// Plot the XY plane
	cv::Mat plane(300, width, CV_8UC3, cv::Scalar(0,0,0)); 

	int cell_n = width/cell_size; 
	pos2d_t orient; 

	// General Reference frame shifted wrt CameraRF on y and z axis

	cv::circle(plane, cv::Point(SpineBase.y*100+width/2, SpineBase.x*100), body, cv::Scalar(255,255,255), 1, 8, 0); 
	cv::circle(plane, cv::Point(ShouldLeft.y*100+width/2, ShouldLeft.x*100), hand, cv::Scalar(0,0,255), 1, 8, 0); 
	cv::circle(plane, cv::Point(ShouldRight.y*100+width/2, ShouldRight.x*100), hand, cv::Scalar(0,0,255), 1, 8, 0);
	//std::cout << "SpinBase : " << cv::Point(SpineBase.y * 100 + width / 2, SpineBase.x * 100) << std::endl;
	// Orientation of the body wrt Camera RF evaluated from the position of the shoulders (components of the vector normal to shoulders line) 
	orient.x = (ShouldLeft.y - ShouldRight.y);                         
	orient.y = -(ShouldLeft.x - ShouldRight.x); 

	cv::circle(plane, cv::Point(EE.y*100+width/2, EE.x*100), hand, cv::Scalar(255, 0, 255), 1, 8, 0); 
	//cv::circle(plane, cv::Point(EE_noF_y*100+width/2, EE_noF_x*100), hand, cv::Scalar(255, 0, 0), 1, 8, 0); 
	//std::cout << "EE : " << cv::Point(EE.y * 100 + width / 2, EE.x * 100) << std::endl;

	cv::circle(plane, cv::Point(EndPoint.y*100 + width / 2, EndPoint.x * 100), hand, cv::Scalar(255, 255, 0), 1, 8, 0);
	cv::circle(plane, cv::Point(FirstLink.y*100 + width / 2, FirstLink.x * 100), hand, cv::Scalar(255, 255, 0), 1, 8, 0);
	cv::line(plane, cv::Point(FirstLink.y*100 + width / 2, FirstLink.x * 100), cv::Point(EndPoint.y*100 + width / 2, EndPoint.x * 100), cv::Scalar(255, 255, 0), 1, 8, 0); 
	cv::line(plane, cv::Point(FirstLink.y*100 + width / 2, FirstLink.x * 100), cv::Point(0 + width / 2, 0), cv::Scalar(255, 255, 0), 1, 8, 0); 

	cv::arrowedLine(plane, cv::Point(SpineBase.y*100+width/2, SpineBase.x*100), cv::Point(SpineBase.y*100+width/2 + orient.y*100, SpineBase.x*100+orient.x*100), cv::Scalar(255,255,255),1,8,0,0.1);    

	cv::imshow("XY plane", plane);

	return; 
}

/*
pos2d_t GetState::BodyToScreen(CameraSpacePoint bodyPoint, int width, int height)
{
    // Calculate the body's position on the screen
	ColorSpacePoint colorPoint = {0}; 
	mapper->MapCameraPointToColorSpace(bodyPoint, &colorPoint);
	
	pos2d_t coord;
	float coordX, coordY; 

	coord.x = static_cast<float>(colorPoint.X);
	coord.y = static_cast<float>(colorPoint.Y); 
	
	return coord; 
}
*/

void GetState::SwapShoulder(pos3d_t LShould, pos3d_t RShould, double* Left, double* Right) {

	    Left[0] = LShould.x; 
		Left[1] = LShould.y; 
		Left[2] = LShould.z; 
		Right[0] = RShould.x; 
		Right[1] = RShould.y; 
		Right[2] = RShould.z;

	pos3d_t Temp; 
	Temp.x = Left[0];
	Temp.y = Left[1];
	Temp.z = Left[2]; 
	
	Left[0] = Right[0];
	Left[1] = Right[1];
	Left[2] = Right[2];

	Right[0] = Temp.x;
	Right[1] = Temp.y;
	Right[2] = Temp.z;


	return; 
}

double Min(double a[][4], int n) {
  int j, min, k=0;
 
  min = a[0][3];
 
  for (j = 0; j < n; j++) {
    if (a[j][3] < min) {
       k = j; 
       min = a[j][3];
    }
  }
 
  return k;
}

void GetState::ReadSampleData() {

	int NowCount = 0;

	//観測時系列データの読み込み
	//一時保存用ストリームを用意
	std::string str;
	std::stringstream ss;
	//csvファイルを読み込み
	std::ifstream ifs("Input/LogKinect.csv");

	if (!ifs) {
		std::cout << "Error:Input data file not found" << std::endl;
		return;
	}

	for (int i = 0; i < SEQ_LENGTH; ++i) {

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);
		//std::cout << "label : " << SampleData[i].label << std::endl;

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataBody[i].x;
		//std::cout << "Bodyx : " << dataBody[i].x << std::endl;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);
		
		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataBody[i].y;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataBody[i].z;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataLShoulder[i].x;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataLShoulder[i].y;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataLShoulder[i].z;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataRShoulder[i].x;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataRShoulder[i].y;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataRShoulder[i].z;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataRElbow[i].x;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataRElbow[i].y;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataRElbow[i].z;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataRHand[i].x;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//1列読み込み
		getline(ifs.seekg(0, std::ios_base::cur), str, ',');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataRHand[i].y;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//改行コードまで読み込む
		getline(ifs.seekg(0, std::ios_base::cur), str, '\n');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//stringstreamから配列に流す
		ss >> dataRHand[i].z;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);

		//std::cout << "current sample : " << i << std::endl;

	}
	//CSVファイルを閉じてファイルへのアクセス権を開放
	ifs.close();

	std::cout << "Finish reading the sample data" << std::endl;


}