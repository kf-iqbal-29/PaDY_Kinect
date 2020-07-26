#include "StdCalculationLib.h"
#include<math.h>
#include<iostream>
#include <time.h>
#include <random>


//乱数の初期化
void StdCalc::InitializeRandomNumber()
{
	srand((unsigned int)time(NULL));
}


//(0,1)の一様乱数生成
double StdCalc::Uniform()
{
	return ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
}

//正規分布に従う乱数取得
double StdCalc::NormalRand(double mu, double sigma)
{
	double z = sqrt(-2.0*log(Uniform())) * sin(2.0*M_PI*Uniform());
	return mu + sigma*z;
}

//2次元ユークリッド距離の計算
double StdCalc::CalcEuclideDistance2d(const Eigen::Vector2d& x1, const Eigen::Vector2d& x2)
{
	Eigen::Vector2d temp;
	temp = x1 - x2;
	return temp.norm();
}

//4次元ユークリッド距離の計算
double StdCalc::CalcEuclideDistance4d(const Eigen::Vector4d& x1, const Eigen::Vector4d& x2)
{
	Eigen::Vector4d temp;
	temp = x1 - x2;
	return temp.norm();
}

//X-2次元ユークリッド距離の計算
double StdCalc::CalcEuclideDistanceX_2d(const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x1, const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x2)
{
	Eigen::Matrix<double, (MODEL_DOF - 2), 1> temp;
	temp = x1 - x2;
	return temp.norm();
}

//X次元ユークリッド距離の計算
double StdCalc::CalcEuclideDistanceXd(const Eigen::Matrix<double, MODEL_DOF, 1>& x1, const Eigen::Matrix<double, MODEL_DOF, 1>& x2)
{
	Eigen::Matrix<double, MODEL_DOF, 1> temp;
	temp = x1 - x2;
	return temp.norm();
}

//2次元マハラノビス距離の計算
double StdCalc::CalcMahalanobisDistance2d(const Eigen::Vector2d& x1, const Eigen::Vector2d& x2, const Eigen::Matrix2d& sigma)
{
	Eigen::Vector2d temp;
	double val;
	temp = x1 - x2;
	val = temp.transpose() * sigma.inverse() * temp;
	return val;

}

//4次元マハラノビス距離の計算
double StdCalc::CalcMahalanobisDistance4d(const Eigen::Vector4d& x1, const Eigen::Vector4d& x2, const Eigen::Matrix4d& sigma)
{
	Eigen::Vector4d temp;
	double val;
	temp = x1 - x2;
	val = temp.transpose() * sigma.inverse() * temp;
	return val;

}

//X-2次元のマハラノビス距離の計算
double StdCalc::CalcMahalanobisDistanceX_2d(const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x1, const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x2, const Eigen::Matrix<double, (MODEL_DOF - 2), (MODEL_DOF - 2)>& sigma)
{
	Eigen::Matrix<double, (MODEL_DOF - 2), 1> temp;
	double val;
	temp = x1 - x2;
	val = temp.transpose() * sigma.inverse() * temp;
	return val;
}

//X次元のマハラノビス距離の計算
double StdCalc::CalcMahalanobisDistanceXd(const Eigen::Matrix<double, MODEL_DOF, 1>& x1, const Eigen::Matrix<double, MODEL_DOF, 1>& x2, const Eigen::Matrix<double, MODEL_DOF, MODEL_DOF>& sigma)
{
	Eigen::Matrix<double, MODEL_DOF, 1> temp;
	double val;
	temp = x1 - x2;
	val = temp.transpose() * sigma.inverse() * temp;
	return val;
}

//1次元正規分布の尤度計算
double StdCalc::Calc1dGaussian(double x, double mu, double sigma)
{
	double val;
	val = (x - mu) * (x - mu) / (2.0 * sigma * sigma);
	return exp(-val) / (sqrt(2.0 * M_PI) * sigma);
}

//2次元正規分布の尤度計算
double StdCalc::Calc2dGaussian(const Eigen::Vector2d& x, const Eigen::Vector2d& mu, const Eigen::Matrix2d& sigma)
{
	double val;
	val = CalcMahalanobisDistance2d(x, mu, sigma);//マハラノビス距離の計算
	return exp(-val / 2) / (2.0 * M_PI *  sqrt(sigma.determinant()));
}

//4次元正規分布の尤度計算
double StdCalc::Calc4dGaussian(const Eigen::Vector4d& x, const Eigen::Vector4d& mu, const Eigen::Matrix4d& sigma)
{
	double val;
	val = CalcMahalanobisDistance4d(x, mu, sigma);//マハラノビス距離の計算
	return exp(-val / 2) / (4.0 * M_PI * M_PI * sqrt(sigma.determinant()));
}

//X-2次元正規分布の尤度計算
double StdCalc::CalcX_2dGaussian(const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x, const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& mu, const Eigen::Matrix<double, (MODEL_DOF - 2), (MODEL_DOF - 2)>& sigma)
{
	double val;
	val = CalcMahalanobisDistanceX_2d(x, mu, sigma);//マハラノビス距離の計算
	return exp(-val / 2) / (pow((2 * M_PI), ((MODEL_DOF - 2) / 2.0)) * sqrt(sigma.determinant()));
}

//X次元正規分布の尤度計算
double StdCalc::CalcXdGaussian(const Eigen::Matrix<double, MODEL_DOF, 1>& x, const Eigen::Matrix<double, MODEL_DOF, 1>& mu, const Eigen::Matrix<double, MODEL_DOF, MODEL_DOF>& sigma)
{
	double val;
	val = CalcMahalanobisDistanceXd(x, mu, sigma);//マハラノビス距離の計算
	return exp(-val / 2) / (pow((2 * M_PI), (MODEL_DOF / 2.0)) * sqrt(sigma.determinant()));
}

//計算可能かどうか判定
int StdCalc::INFCheck(double x, int mode)
{
	if (!isfinite(x)){
		if(mode) std::cout << "error,  val = " << x << std::endl;
		return 1;
	}

	return 0;
}

double StdCalc::max(double a, double b)
{
	if (a>b)
		return a;
	else
		return b;
}