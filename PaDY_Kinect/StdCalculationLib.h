/* 基本計算用ライブラリ */
#pragma once
#define _USE_MATH_DEFINES
#define EIGEN_NO_DEBUG
#include <Eigen/Core>
#include <Eigen/LU>


#define REGRESSION_ORDER 4//自己回帰モデルの次数
#define MODEL_DOF ((REGRESSION_ORDER)*2)//状態変数の次元(x,yのため必ず偶数にすること)

#ifndef DEG2RAD
#define DEG2RAD    (M_PI/180.0)
#endif

#ifndef Deg2Rad
#define Deg2Rad    (M_PI/180.0)
#endif

#ifndef RAD2DEG
#define RAD2DEG    (180.0/M_PI)
#endif

#ifndef Rad2Deg
#define Rad2Deg    (180.0/M_PI)
#endif

namespace StdCalc{

	//乱数の初期化
	void InitializeRandomNumber();

	//(0,1)の一様乱数生成
	double Uniform();
	
	//正規分布に従う乱数取得
	double NormalRand(double mu, double sigma);

	//2次元のユークリッド距離の計算
	double CalcEuclideDistance2d(const Eigen::Vector2d& x1, const Eigen::Vector2d& x2);

	//4次元のユークリッド距離の計算
	double CalcEuclideDistance4d(const Eigen::Vector4d& x1, const Eigen::Vector4d& x2);

	//X-2次元のユークリッド距離の計算
	double CalcEuclideDistanceX_2d(const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x1, const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x2);

	//X次元のユークリッド距離の計算
	double CalcEuclideDistanceXd(const Eigen::Matrix<double, MODEL_DOF, 1>& x1, const Eigen::Matrix<double, MODEL_DOF, 1>& x2);

	//2次元のマハラノビス距離の計算
	double CalcMahalanobisDistance2d(const Eigen::Vector2d& x1, const Eigen::Vector2d& x2, const Eigen::Matrix2d& sigma);

	//4次元のマハラノビス距離の計算
	double CalcMahalanobisDistance4d(const Eigen::Vector4d& x1, const Eigen::Vector4d& x2, const Eigen::Matrix4d& sigma);

	//X-2次元のマハラノビス距離の計算
	double CalcMahalanobisDistanceX_2d(const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x1, const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x2, const Eigen::Matrix<double, (MODEL_DOF - 2), (MODEL_DOF - 2)>& sigma);

	//X次元のマハラノビス距離の計算
	double CalcMahalanobisDistanceXd(const Eigen::Matrix<double, MODEL_DOF, 1>& x1, const Eigen::Matrix<double, MODEL_DOF, 1>& x2, const Eigen::Matrix<double, MODEL_DOF, MODEL_DOF>& sigma);

	//1次元正規分布の尤度計算
	double Calc1dGaussian(double x, double mu, double sigma);

	//2次元正規分布の尤度計算(共分散行列が対角行列と仮定した場合)←今後拡張必須
	double Calc2dGaussian(const Eigen::Vector2d& x, const Eigen::Vector2d& mu, const Eigen::Matrix2d& sigma);
	
	//4次元正規分布の尤度計算(共分散行列が対角行列と仮定した場合)←今後拡張必須
	double Calc4dGaussian(const Eigen::Vector4d& x, const Eigen::Vector4d& mu, const Eigen::Matrix4d& sigma);

	//X-2次元正規分布の尤度計算
	double CalcX_2dGaussian(const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x, const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& mu, const Eigen::Matrix<double, (MODEL_DOF - 2), (MODEL_DOF - 2)>& sigma);

	//X次元正規分布の尤度計算
	double CalcXdGaussian(const Eigen::Matrix<double, MODEL_DOF, 1>& x, const Eigen::Matrix<double, MODEL_DOF, 1>& mu, const Eigen::Matrix<double, MODEL_DOF, MODEL_DOF>& sigma);
	
	//計算可能かどうか判定
	int INFCheck(double x, int mode);

	double max(double a, double b);

};//namespace StdCalc

