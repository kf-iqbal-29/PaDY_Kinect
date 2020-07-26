/* ��{�v�Z�p���C�u���� */
#pragma once
#define _USE_MATH_DEFINES
#define EIGEN_NO_DEBUG
#include <Eigen/Core>
#include <Eigen/LU>


#define REGRESSION_ORDER 4//���ȉ�A���f���̎���
#define MODEL_DOF ((REGRESSION_ORDER)*2)//��ԕϐ��̎���(x,y�̂��ߕK�������ɂ��邱��)

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

	//�����̏�����
	void InitializeRandomNumber();

	//(0,1)�̈�l��������
	double Uniform();
	
	//���K���z�ɏ]�������擾
	double NormalRand(double mu, double sigma);

	//2�����̃��[�N���b�h�����̌v�Z
	double CalcEuclideDistance2d(const Eigen::Vector2d& x1, const Eigen::Vector2d& x2);

	//4�����̃��[�N���b�h�����̌v�Z
	double CalcEuclideDistance4d(const Eigen::Vector4d& x1, const Eigen::Vector4d& x2);

	//X-2�����̃��[�N���b�h�����̌v�Z
	double CalcEuclideDistanceX_2d(const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x1, const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x2);

	//X�����̃��[�N���b�h�����̌v�Z
	double CalcEuclideDistanceXd(const Eigen::Matrix<double, MODEL_DOF, 1>& x1, const Eigen::Matrix<double, MODEL_DOF, 1>& x2);

	//2�����̃}�n���m�r�X�����̌v�Z
	double CalcMahalanobisDistance2d(const Eigen::Vector2d& x1, const Eigen::Vector2d& x2, const Eigen::Matrix2d& sigma);

	//4�����̃}�n���m�r�X�����̌v�Z
	double CalcMahalanobisDistance4d(const Eigen::Vector4d& x1, const Eigen::Vector4d& x2, const Eigen::Matrix4d& sigma);

	//X-2�����̃}�n���m�r�X�����̌v�Z
	double CalcMahalanobisDistanceX_2d(const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x1, const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x2, const Eigen::Matrix<double, (MODEL_DOF - 2), (MODEL_DOF - 2)>& sigma);

	//X�����̃}�n���m�r�X�����̌v�Z
	double CalcMahalanobisDistanceXd(const Eigen::Matrix<double, MODEL_DOF, 1>& x1, const Eigen::Matrix<double, MODEL_DOF, 1>& x2, const Eigen::Matrix<double, MODEL_DOF, MODEL_DOF>& sigma);

	//1�������K���z�̖ޓx�v�Z
	double Calc1dGaussian(double x, double mu, double sigma);

	//2�������K���z�̖ޓx�v�Z(�����U�s�񂪑Ίp�s��Ɖ��肵���ꍇ)������g���K�{
	double Calc2dGaussian(const Eigen::Vector2d& x, const Eigen::Vector2d& mu, const Eigen::Matrix2d& sigma);
	
	//4�������K���z�̖ޓx�v�Z(�����U�s�񂪑Ίp�s��Ɖ��肵���ꍇ)������g���K�{
	double Calc4dGaussian(const Eigen::Vector4d& x, const Eigen::Vector4d& mu, const Eigen::Matrix4d& sigma);

	//X-2�������K���z�̖ޓx�v�Z
	double CalcX_2dGaussian(const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& x, const Eigen::Matrix<double, (MODEL_DOF - 2), 1>& mu, const Eigen::Matrix<double, (MODEL_DOF - 2), (MODEL_DOF - 2)>& sigma);

	//X�������K���z�̖ޓx�v�Z
	double CalcXdGaussian(const Eigen::Matrix<double, MODEL_DOF, 1>& x, const Eigen::Matrix<double, MODEL_DOF, 1>& mu, const Eigen::Matrix<double, MODEL_DOF, MODEL_DOF>& sigma);
	
	//�v�Z�\���ǂ�������
	int INFCheck(double x, int mode);

	double max(double a, double b);

};//namespace StdCalc

