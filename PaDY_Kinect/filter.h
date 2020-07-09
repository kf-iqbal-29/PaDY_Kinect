#pragma once

#include <time.h>
#include <opencv2/opencv.hpp>


class Filter {
 private:
	 cv::Matx22d A;
	 cv::Matx<double, 1,2> H;
	 cv::Matx22d P;
	 cv::Matx22d V;
	 cv::Matx<double,1,1> W;
	 cv::Matx<double, 2, 1> x_cur;
	 cv::Matx<double, 2, 1> x_old;
	 cv::Matx<double,1,1> z; 
     cv::Matx<double, 2, 1> K;


 public:
	Filter::Filter() { 
   
    // State Vector: X Y Xdot Ydot
	
	H = cv::Matx<double, 1,2>(1, 0);   //Measurement matrix
	P = cv::Matx22d(100000, 0,
                     0, 100000);          //Initial covariance matrix
	V = cv::Matx22d(1, 0,
                    0, 1);           //Process noise matrix
	W = 10;                 //measurement noise matrix
	x_cur = cv::Matx<double, 2, 1>(0, 0);           //state vector

	 }

	Filter::~Filter() {}

	void Filter::Prediction(double DeltaT) {

    x_old = x_cur;                                 //Save old state

	A = cv::Matx22d (1, DeltaT, 
                     0, 1);                  //Current step dynamic matrix 

	x_cur = A*x_old;                               //Predict next step state
         
    P = (A*P)*A.t() + V;                           //Predict next step covariance
  }


	void Filter::Update(double x) {     //Input new measurements
    
    K = (P*H.t())*((H*P)*H.t() + W).inv();         //Compute kalman gain matrix
	
    z = x;       
    
    x_cur = x_cur + K*(z - H*x_cur);              //Update state prediction

	P = P - (K*H)*P;                              //Update covariance
	   /* std::cout << "p. " << P(0, 0) << P(1, 0) << P(2, 0) << P(3, 0) << "\n";
		std::cout << "p. " << P(0, 1) << P(1, 1) << P(2, 1) << P(3, 1) << "\n";
		std::cout << "p. " << P(0, 2) << P(1, 2) << P(2, 2) << P(3, 2) << "\n";
		std::cout << "p. " << P(0, 3) << P(1, 3) << P(2, 3) << P(3, 3) << "\n";*/

  }

	void Filter::State(double * data) {           //pass array to write final estimate

    data[0] = x_cur(0);
   // data[1] = x_cur(1);

  }

};



class Linear {

private: 
	cv::Vec6d X;
	cv::Vec4d X2; 
	cv::Vec2d Y;
	cv::Vec6d Z;
	cv::Vec4d Z2; 

public:

	Linear::Linear() {                           // Linear Filter with 10 elements for both X and Y coordinate
	    Y[0]=Y[1] = 0; 
	    X[0]=X[1]=X[2]=X[3]=X[4]=X[5]=0;  
	    X2[0]=X2[1]=X2[2]=X2[3] = 0; 
	    Z2[0]=Z2[1]=Z2[2]=Z2[3] = 0; }
	
	Linear::~Linear() {}

	void Linear::Start(double x_new, double z_new) {
 
		int j; 

		for(j = 0; j<6; j++) {         // Shift values
				X[j] = X[j+1];  }
		X[5] = X2[0]; 
		for(j = 0; j<4; j++) {         // Shift values
				X2[j] = X2[j+1];  }

		for(j = 0; j<6; j++) {         // Shift values
				Z[j] = Z[j+1];  }
		Z[5] = Z2[0]; 
		for(j = 0; j<4; j++) {         // Shift values
				Z2[j] = Z2[j+1];  }
		

		X2[3] = x_new;
		Z2[3] = z_new; 

		Y[1] = (X[0]+X[1]+X[2]+X[3]+X[4]+X[5]+X2[0]+X2[1]+X2[2]+X2[3])/10 ; 
		Y[0] = (Z[0]+Z[1]+Z[2]+Z[3]+Z[4]+Z[5]+Z2[0]+Z2[1]+Z2[2]+Z2[3])/10 ; 

			
	}

	void Linear::Out(double * data) {           //pass array to write final estimate

		data[0] = Y[1]; 
		data[1] = Y[0]; 
		//std::cout<<Y[1]<<std::endl;
	}
 
};



//class Filter {
// private:
//	 cv::Matx44d A;
//	 cv::Matx<double, 2, 4> H;
//	 cv::Matx44d P;
//	 cv::Matx44d V;
//	 cv::Matx22d W;
//	 cv::Vec4d x_cur;
//	 cv::Vec4d x_old;
//	 cv::Vec2d z; 
//     cv::Matx<double, 4, 2> K;
//
//
// public:
//	Filter::Filter( double x_1, double x_2, double x_3, double x_4) { 
//   
//    // State Vector: X Y Xdot Ydot
//	
//	H = cv::Matx<double, 2, 4>(1, 0, 0, 0,
//                               0, 1, 0, 0);   //Measurement matrix
//	P = cv::Matx44d(10, 0, 0, 0,
//                     0, 10, 0, 0,
//                     0, 0, 10, 0,
//                     0, 0, 0, 10);          //Initial covariance matrix
//	V = cv::Matx44d(5, 0, 0, 0,
//                    0, 5, 0, 0,
//                    0, 0, 10, 0,
//                    0, 0, 0, 10);           //Process noise matrix
//	W = cv::Matx22d(1, 0,
//                    0, 1);                 //measurement noise matrix
//	x_cur = cv::Vec4d(x_1, x_2, x_3, x_4);           //state vector
//	
//	 }
//
//	Filter::~Filter() {}
//
//	void Filter::Prediction(double DeltaT) {
//
//    x_old = x_cur;                                 //Save old state
//
//	A = cv::Matx44d (1, DeltaT, 0, 0,
//                     0, 1, 0, DeltaT,
//                     0, 0, 1, 0,
//				     0, 0, 0, 1);                  //Current step dynamic matrix 
//
//	x_cur = A*x_old;                               //Predict next step state
//         
//    P = (A*P)*A.t() + V;                           //Predict next step covariance
//  }
//
//
//	void Filter::Update(double x, double y) {     //Input new measurements
//    
//    K = (P*H.t())*((H*P)*H.t() + W).inv();         //Compute kalman gain matrix
//	
//    z[0] = x; 
//    z[1] = y;      
//    
//    x_cur = x_cur + K*(z - H*x_cur);              //Update state prediction
//
//	P = P - (K*H)*P;                              //Update covariance
//	   // std::cout << "p. " << P(0, 0) << P(1, 0) << P(2, 0) << P(3, 0) << "\n";
//		//std::cout << "p. " << P(0, 1) << P(1, 1) << P(2, 1) << P(3, 1) << "\n";
//		//std::cout << "p. " << P(0, 2) << P(1, 2) << P(2, 2) << P(3, 2) << "\n";
//		//std::cout << "p. " << P(0, 3) << P(1, 3) << P(2, 3) << P(3, 3) << "\n";
//
//  }
//
//	void Filter::State(double * data) {           //pass array to write final estimate
//
//    data[0] = x_cur[0];
//    data[1] = x_cur[1];
//
//  }
//
//};