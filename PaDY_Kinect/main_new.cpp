#include <Windows.h>
#include <process.h>
#include <conio.h>
#include<sstream>
#include<fstream>
#include<string>
#include "GetState.h"
#include "filter.h"
#include "MPC.h"
#include "ArmParams.h"
//#include "GMM.h" //20200628

#define SAMPLING_TIME 30 //(ms)
#define RSTART_TIMING 100

#define INITIAL_J1 (90 * DEG2RAD)		// joint angles of the two arms of paddy 
#define INITIAL_J2 (-180 * DEG2RAD)
#define DELIVERY_TIME 2.0  //(s)
#define OFFSET_X 0.7       //(m)
#define OFFSET_Y 0.5       //(m)
#define LOG_SIZE 100000

#define distance_for_termination 0.2  // The distance between end-effector and worker body that would fulfill the termination condition
#define velocity_for_termination 0.05 // The velocity difference between end-effector velocity and desired velocity that would fulfill the termination condition

GetState* kinect;//Object for Kinect and linear filter
Linear* MAF;	// MAF = Moving Average Filter
Linear* MAF2; 
MPC mpc;

unsigned long Start_t, End_t; 
double Arm_Length; 
pos2d_t FirstLink; 
pos3d_t Delivery_Pos;
pos3d_t RShould, LShould, LHand, RHand, LElbow, RElbow, Body; ;
pos3d_t End_Pos;
pos4d_t Output;


// For writing the log of robot
int log_counter = 0;
Eigen::Vector2d RobotAngleLog[LOG_SIZE];
Eigen::Vector2d RobotAnVelLog[LOG_SIZE];
Eigen::Vector2d RobotEndPosLog[LOG_SIZE];
Eigen::Vector2d RobotEndVelLog[LOG_SIZE];
Eigen::Vector2d RobotDelPosLog[LOG_SIZE];
Eigen::Vector2d RobotAnAccLog[LOG_SIZE]; 

void WriteLog()
{

	std::cout << "----------  Log Start  ----------" << std::endl;

	//Open the csv file
	std::ofstream ofs("LogData.csv");
	ofs << "Time" << "," << "Angle1 [deg]" << "," << "Angle2 [deg]"
		<< "," << "Angle1 vel [deg/s]" << "," << "Angle2 vel [deg/s]"
		<< "," << "Angle1 acc [deg/s^2]" << "," << "Angle2 acc [deg/s^2]"
		<< "," << "Endpoint pos x [m]" << "," << "Endpoint pos y [m]"
		<< "," << "Endpoint vel x [m/s]" << "," << "Endpoint vel y [m/s]"
		<< "," << "Delivery pos x [m]" << "," << "Delivery pos y [m]";
	ofs << std::endl;
	for (int i = 0; i < log_counter; ++i){
		ofs << i * SAMPLING_TIME * 0.001 << ","
			<< RobotAngleLog[i](0) << "," << RobotAngleLog[i](1) << ","
			<< RobotAnVelLog[i](0) << "," << RobotAnVelLog[i](1) << ","
			<< RobotAnAccLog[i](0) << "," << RobotAnAccLog[i](1) << ","
			<< RobotEndPosLog[i](0) << "," << RobotEndPosLog[i](1) << ","
			<< RobotEndVelLog[i](0) << "," << RobotEndVelLog[i](1) << ","
			<< RobotDelPosLog[i](0) << "," << RobotDelPosLog[i](1) << ","
			<< std::endl;
	}
	ofs.close();


	std::cout << "----------  Log Finish  ----------" << std::endl;
}

// For writing the Log of Kinect
unsigned int kinect_counter = 0;
pos3d_t KinectData_LShould[LOG_SIZE];
pos3d_t KinectData_RShould[LOG_SIZE];
pos3d_t KinectData_RElbow[LOG_SIZE];
pos3d_t KinectData_RHand[LOG_SIZE];
pos3d_t KinectData_Body[LOG_SIZE];

void WriteKinectLog()
{
	std::cout << "----------  Kinect Log Start  ----------" << std::endl;

	//Open the csv file
	std::ofstream ofs("KinectLogData.csv");
	ofs << "Time" << "," << "Body[x]" << "," << "Body[y]" << "," << "Body[z]"
		<< "," << "Left Shoulder [x]" << "," << "Left Shoulder [y]" << "," << "Left Shoulder [z]"
		<< "," << "Right Shoulder [x]" << "," << "Right Shoulder [y]" << "," << "Right Shoulder [z]"
		<< "," << "Right Elbow [x]" << "," << "Right Elbow [y]" << "," << "Right Elbow [z]"
		<< "," << "Right Hand [x]" << "," << "Right Hand [y]" << "," << "Right Hand [z]";

	ofs << std::endl;
	for (int i = 0; i < kinect_counter; ++i) {
		ofs << i * SAMPLING_TIME * 0.001 << ","
			<< KinectData_Body[i].x << "," << KinectData_Body[i].y << "," << KinectData_Body[i].z << ","
			<< KinectData_LShould[i].x << "," << KinectData_LShould[i].y << "," << KinectData_LShould[i].z << ","
			<< KinectData_RShould[i].x << "," << KinectData_RShould[i].y << "," << KinectData_RShould[i].z << ","
			<< KinectData_RElbow[i].x << "," << KinectData_RElbow[i].y << "," << KinectData_RElbow[i].z << ","
			<< KinectData_RHand[i].x << "," << KinectData_RHand[i].y << "," << KinectData_RHand[i].z << ","
			<< std::endl;
	}
	ofs.close();


	std::cout << "----------  Kinect Log Finish  ----------" << std::endl;

}

//Show the color image
void Display(LPVOID	pParam)
{
	while (1){

		cv::Mat frame; 

		//kinect->UpdateBodyFrame(); //Update BodyFrame

		kinect->UpdateColorFrame();//Update ColorFrame

		int key = cv::waitKey(30);
		if (key == 27)	break;
	}
}


void GetState::UpdateColorFrame()   
{
	/*
	if (!ColorFrameReader){
		std::cout << "Error : ColorFrameReader" << std::endl;
		return;
	}
	
	//Get ColorFrame
	ComPtr<IColorFrame> ColorFrame;
	auto ret = ColorFrameReader->AcquireLatestFrame(&ColorFrame);
	*/
	bool ret = S_OK;
	if (ret == S_OK){
		/*
		ComPtr<IFrameDescription> ColorFrameDescription;

		//Get the size of color image
		ColorFrame->get_FrameDescription(&ColorFrameDescription);
		ColorFrameDescription->get_Width(&Width);
		ColorFrameDescription->get_Height(&Height);
		*/
		//Make the buffer (for RGBA)
		//std::vector<BYTE> Buffer;
		//int Width = 1920;
		//int Height = 1080;
		//Buffer.resize(Width * Height * 4);

		//Get the color image data
		//ColorFrame->CopyConvertedFrameDataToArray(Buffer.size(), &Buffer[0], ColorImageFormat::ColorImageFormat_Bgra);

		//color image
		//cv::Mat colorImage(Height, Width, CV_8UC4, &Buffer[0]);
		
		//Make the window
		//cv::resize(colorImage, colorImage, cv::Size(), 0.5, 0.5);//Resize

		// Extract position of the hands and shoulders in the camera image
		//pos2d_t Left_H, Right_H, Left_S, Right_S, Left_E, Right_E;
		
		//CameraSpacePoint HandLeft = kinect->joints[JointType_HandLeft].Position; 
		//CameraSpacePoint HandRight = kinect->joints[JointType_HandRight].Position; 
		//Left_H = kinect->BodyToScreen(HandLeft, Width, Height); 
		//Right_H = kinect->BodyToScreen(HandRight, Width, Height);
		//Right_E = kinect->BodyToScreen(kinect->joints[JointType_ElbowRight].Position, Width, Height);
		//Left_E = kinect->BodyToScreen(kinect->joints[JointType_ElbowLeft].Position, Width, Height);
		//Left_S = kinect->BodyToScreen(kinect->joints[JointType_ShoulderLeft].Position, Width, Height);
		//Right_S = kinect->BodyToScreen(kinect->joints[JointType_ShoulderRight].Position, Width, Height);
		
		// Draw arms 
		//colorImage = kinect->Draw(Left_H, Left_E, Left_S, colorImage); 
		//colorImage = kinect->Draw(Right_H, Right_E, Right_S, colorImage); 

		//Show the color image
		//cv::imshow("GetState", colorImage);
    

		// Show the position of the worker on the grid
		//pos3d_t Position = kinect->TransKinectToRealWorld3D(kinect->joints[JointType_SpineBase].Position); 
		//pos3d_t HandPos = kinect->TransKinectToRealWorld3D(HandLeft); 
		//pos3d_t HandPosR = kinect->TransKinectToRealWorld3D(HandRight);

		//Offset the values
		//Position.x = Position.x + OFFSET_X;
		//Position.y = Position.y + OFFSET_Y;


		kinect->UpSideVision(Body, LShould, RShould, Delivery_Pos, Output.x, Output.y, End_Pos, FirstLink); 

		return; 
	}
}


int main()
{
	bool useKinect = false; // do you want to use Kinect in real time? or do you want to use Kinect Log Data?
	int j, i = 0;
	int cnt = 0;
	double Filt[2], Array[20][4], BodyValues[3], R_temp[3], L_temp[3];
	double desTime, outTime = 0;
	Eigen::Vector4d des_x;  //Target state
	Eigen::Vector4d f;      // f is system state [x,y,x_dot,y_dot] 
	Eigen::Vector4d cur_x;  //Current state
	Eigen::Vector2d input;  //Input vector
	Eigen::Vector2d end_pos, end_vel;// EndEffector state

	/*bool use_GMM = true; //20200628 Do you want to use GMM?

	GMM Worker_GMM; //20200628 
	int sample_no = 0; //20200628 The samples under consideration for GMMM
	sample_no = Worker_GMM.ReadModelInfo();  // 20200628 Read Model info to initialize GMM
	Eigen::Vector2d pos; //20200628 Current position of worker (Taken from Kinect sensor reading)
	static std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> temp;	//20200628 For temporary storage of observation points
	*/
	//int faces[BODY_COUNT]; 
	int BackP = 0; 
	double FTrackState = 0;//Tracking state of face
	double FaceFilt[2]; 

	HANDLE hThread;//Handle of thread for showing the color image

	//Generate the object for Kinect and filters
	kinect = new GetState();  
	MAF = new Linear(); 
	MAF2 = new Linear();	


	//Make the thread for showing the color image
	hThread = (HANDLE)_beginthread(Display, 0, NULL);

	//Initialize the robot position
	cur_x(0) = INITIAL_J1;
	cur_x(1) = INITIAL_J2;
	cur_x(2) = 0.0;
	cur_x(3) = 0.0;
	desTime = DELIVERY_TIME;
	End_Pos.x = 0.0;
	End_Pos.y = 0.128;
	FirstLink.x = L1*cos(INITIAL_J1); 
	FirstLink.y = L1*sin(INITIAL_J1); 

	cnt = 0; 

	// If you are using kinect then initialize it. Otherwise load data from a log file.
	if (useKinect)
	{
		if (!kinect->Initialize()) return -1;
	}
	else
		kinect->ReadSampleData();

	// Get the data as long as a key is not pressed 
	while (!_kbhit()){

		Start_t = timeGetTime();

		if (useKinect)
		{
			kinect->UpdateBodyFrame(); //Update BodyFrame

		   /* Left arm position */
		   //LHand = kinect->TransKinectToRealWorld3D(kinect->joints[JointType_HandLeft].Position);
		   LShould = kinect->TransKinectToRealWorld3D(kinect->joints[JointType_ShoulderLeft].Position);
		   //LElbow = kinect->TransKinectToRealWorld3D(kinect->joints[JointType_ElbowLeft].Position); 

		   /* Right arm position */
		   RHand = kinect->TransKinectToRealWorld3D(kinect->joints[JointType_HandRight].Position);
		   RShould = kinect->TransKinectToRealWorld3D(kinect->joints[JointType_ShoulderRight].Position); 
		   RElbow = kinect->TransKinectToRealWorld3D(kinect->joints[JointType_ElbowRight].Position); 

		   Body = kinect->TransKinectToRealWorld3D(kinect->joints[JointType_SpineBase].Position);

		   // Write Log for Kinect Data
		   KinectData_Body[kinect_counter] = Body;
		   KinectData_LShould[kinect_counter] = LShould;
		   KinectData_RShould[kinect_counter] = RShould;
		   KinectData_RElbow[kinect_counter] = RElbow;
		   KinectData_RHand[kinect_counter] = RHand;
		   kinect_counter++;
		   if (kinect_counter >= LOG_SIZE)
			   kinect_counter = LOG_SIZE - 1;
		   
		}
		else
		{
			LShould = kinect->dataLShoulder[cnt];
			RHand = kinect->dataRHand[cnt];
			RShould = kinect->dataRShoulder[cnt];
			RElbow = kinect->dataRElbow[cnt];
			Body = kinect->dataBody[cnt];
		}


		/* Evaluate Arm length just at the beginning */ 
		double length_R = kinect->ArmLength(RShould, RElbow, RHand); 
		Arm_Length = length_R; 
	 
		if(BackP>60)
		FTrackState = 0; 
		else
		FTrackState = 1; 

		for(j=0; j<10; j++) {
			MAF2->Start(FTrackState, 0); 
		    MAF2->Out(FaceFilt);        }

		FTrackState = FaceFilt[0]; 

		if(FTrackState < 0.2) {					// Why 0.2?
		//std::cout<<"No Face"<<std::endl; 
		kinect->SwapShoulder(LShould, RShould,L_temp,R_temp);	    
		LShould.x = L_temp[0];
		LShould.y = L_temp[1];
		LShould.z = L_temp[2];
		RShould.x = R_temp[0];
		RShould.y = R_temp[1];
		RShould.z = R_temp[2]; }

        /* Analize different branches of the RRT for the same detected body position */
		/*for(j=0; j<20; j++) {
			Output = kinect->TRRT(LShould, RShould, Body, Arm_Length, BodyValues); 
		Array[j][0] = Output.x; 
		Array[j][1] = Output.y; 
		Array[j][2] = Output.z; 
		Array[j][3] = Output.w; 
		}

		int k = Min(Array, 20); 
		Output.x = Array[k][0] ;//+ OFFSET_X;
		Output.y = Array[k][1] ;//+ OFFSET_Y;
		Output.z = Array[k][2];
		Output.w = Array[k][3];
		*/
		//std::cout<<"Minimum in array "<<Output.x<<","<<Output.y<<","<<Output.z<<" Cost "<<Output.w<<std::endl;
		double Elaps_ms = timeGetTime() - Start_t; 

		//Filt[0] = Output.x;
		//Filt[1] = Output.y;
		
		/* Apply Moving Average Filter */
		/*for(j=0; j<20; j++) {
			MAF->Start(Filt[0], Filt[1]); 
			MAF->Out(Filt);                 
		}
		*/
		//Elaps_ms = timeGetTime() - Start_t; 
		//std::cout<<"Milliseconds = "<<Elaps_ms<<std::endl; 

		//Delivery_Pos.x = Filt[0];
		//Delivery_Pos.y = Filt[1];

		/*
		// 20200628 For GMM
		if (use_GMM == true)
		{
			pos(0) = Output.x;
			pos(1) = Output.y;
			temp.push_back(pos);
			Worker_GMM.SetSampleData(pos);
			Worker_GMM.EMAlgorithm();

			Delivery_Pos.x = Worker_GMM.Model[1].mu(0);
			Delivery_Pos.y = Worker_GMM.Model[1].mu(1);
		}
		*/

		/* Show data */
		if (i % 10 == 0){
			//std::cout << "Left hand position x:" << LHand.x << ", y:" << LHand.y << ", z:" << LHand.z << std::endl;
			//std::cout << "Right hand position x:" << RHand.x << ", y:" << RHand.y << ", z:" << RHand.z << std::endl;
			//std::cout << "Left shoulder position x:" << LShould.x<<", y:" <<LShould.y <<", z:"<< LShould.z<< std::endl;
			//std::cout << "Right shoulder position x:" << RShould.x<<", y:" <<RShould.y <<", z:"<< RShould.z<< std::endl;
	        std::cout<<"Arm Lenght "<<Arm_Length<<std::endl;
			//std::cout<<"Milliseconds = "<<Elaps_ms<<std::endl;
			//std::cout<<"Delivery Position x:"<<Delivery_Pos.x<<" y:"<<Delivery_Pos.y<<std::endl; 
			//std::cout<< "Total cost: "<<Output.w<<std::endl;
			std::cout << "cnt: " << cnt << std::endl;
			std::cout<<std::endl; 
		}
		
		if(i == RSTART_TIMING) { std::cout<<"Planning Started"<<std::endl<<std::endl; } 

		/****** Trajectory planning part  ******/
		if ((i > RSTART_TIMING) && (desTime > 0)){

			/* Robot starts moving */
			//std::cout << "Robot Moving" << std::endl;
			//Set the target position
			//des_x(0) = Delivery_Pos.x;
			//des_x(1) = Delivery_Pos.y;
			des_x(0) = 0.0;
			des_x(1) = 0.0;
			des_x(2) = 0.0;
			des_x(3) = 0.0;
			mpc.setDestination(des_x);

			//Check the convergence: don't stop the algorithm if destination is not reached
			mpc.JointAngle2EndPosition(end_pos, cur_x.block<2, 1>(0, 0));  
			// Block command: cur_x.block<2, 1>(0, 0): Start from (0,0) of the matrix cur_x and take <2,1> elements i.e. 2 rows 1 column
			mpc.JointVel2EndVel(end_vel, cur_x.block<2, 1>(0, 0), cur_x.block<2, 1>(2, 0));

			// Check the Termination conditions
			if ((fabs(Body.x - end_pos(0)) < distance_for_termination) && (fabs(Body.y - end_pos(1)) < distance_for_termination)
				&& (fabs(des_x(2) - end_vel(0)) < velocity_for_termination) && (fabs(des_x(3) - end_vel(1)) < velocity_for_termination))
			{ 
				desTime = -1;
				std::cout<<"Task completed"<<std::endl; 
				std::cout<<"Extra time (s): "<<(outTime*SAMPLING_TIME*0.001 - 0.15)<<std::endl<<std::endl; 
			}
			else
			{
				desTime = DELIVERY_TIME - (i - RSTART_TIMING) * SAMPLING_TIME * 0.001;
				if (desTime < 0.15)
				{
					 desTime = 0.15;	
				     outTime++; 
				}
			}

			//Generate the trajectory by model predictive control
			if (mpc.calcInputbyGradientMethod(input, cur_x, desTime, SAMPLING_TIME * 0.001, LShould, RShould, Body, Arm_Length, BodyValues))
			{
				mpc.calcStateFunc(f, cur_x, input);			// f is system state X = [x,y,x_dot,y_dot] 
				cur_x = cur_x + f * SAMPLING_TIME * 0.001;      
			}
			else
			{
				input.setZero();
			} 

			//Check the current position of endpoint
			mpc.JointAngle2EndPosition(end_pos, cur_x.block<2, 1>(0, 0));
			mpc.JointVel2EndVel(end_vel, cur_x.block<2, 1>(0, 0), cur_x.block<2, 1>(2, 0));
			

			End_Pos.x = end_pos(0);
			End_Pos.y = end_pos(1);
			FirstLink.x = L1*cos(cur_x(0)); 
			FirstLink.y = L1*sin(cur_x(0));

			//Get the log data
			RobotAngleLog[log_counter] = cur_x.block<2, 1>(0, 0);  //Robot's joint angles
			RobotAnVelLog[log_counter] = cur_x.block<2, 1>(2, 0);  //Robot's joint angular velocities
			RobotEndPosLog[log_counter] = end_pos;                 //Robot's endeffector position
			RobotEndVelLog[log_counter] = end_vel;                 //Robot's endeffector velocity
			RobotAnAccLog[log_counter] = f.block<2, 1>(2, 0);      //Robot's joint angular acceleration
			RobotDelPosLog[log_counter] = des_x.block<2, 1>(0, 0); //Target delivery position 
			log_counter++;
			if (log_counter >= LOG_SIZE)
				log_counter = LOG_SIZE - 1;
		}



		/* For keeping the constant calculation time */
		double CurrentTime;
		while (1){
			End_t = timeGetTime();
			CurrentTime = End_t - Start_t;
			if (SAMPLING_TIME <= CurrentTime)
				break;
		}
		
		i++;

		cnt++;
		if (cnt >= kinect->SEQ_LENGTH) cnt = kinect->SEQ_LENGTH - 1;
		Sleep(SAMPLING_TIME);
	}

	//Delete the object for Kinect
	delete kinect;
	delete MAF;
	delete MAF2; 
	mpc.~MPC();
	
	std::cout<<"Destroyed"<<std::endl; 

	//Write the log data in csv file
	WriteLog();

	if (useKinect)
		WriteKinectLog();

	return 0;	
}



