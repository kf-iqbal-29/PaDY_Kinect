/* 作業者位置モデル化と作業進度推定のためのGMM用クラス (Class for GMM for worker position modeling and work progress estimation)*/

#pragma once
#include <vector>
#include "StdCalculationLib.h"


#ifndef _GMM
#define _GMM
typedef struct gmm{
	Eigen::Vector2d mu;
	Eigen::Matrix2d sigma;
	double c;
}gmm_t;
#endif

class GMM
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	GMM();//コンストラクタ (Constructor)
	~GMM();//デストラクタ (Destructor)

	static const int TASK_NUM = 3;//GMMのクラスタの数 (Number of GMM clusters) // It was 3 in Kanazawa's program
	static const double INITIAL_POSITION[TASK_NUM + 1][2];//初期作業位置 (Initial working position)

	std::vector<gmm_t, Eigen::aligned_allocator<gmm_t>> Model;//モデル用基本オブジェクト (Basic object for model)

	void InitializeSampleData(){ SampleData.resize(0); }//サンプルデータの初期化 (Initialization of sample data)
	void SetSampleData(const Eigen::Vector2d& x){ SampleData.push_back(x);}//サンプルデータのセット (Set of sample data)

	void EMAlgorithm();//EMアルゴリズムによるGMMの更新 (GMM update by EM algorithm)
	double CalcMixtureGaussian2d(const Eigen::Vector2d& x, const std::vector<gmm_t, Eigen::aligned_allocator<gmm_t>>& gmm);//2次元の混合ガウス分布の尤度を返す (Give the likelihood of the two-dimensional Gaussian mixture distribution)
	double CalcGaussianProb(const Eigen::Vector2d& x, int No);
	int ReadModelInfo();//モデル情報の読み込み (Load model information)
	void WriteModelInfo(int cur_sample);//モデル情報の書き込み (Writing model information)

private:

	static const double INI_SIGMA;//初期分散 (Initial variance)
	static const double ALPHA;//EMアルゴリズムのパラメータ更新用 (For updating the parameters of the EM algorithm)

	int counter = 0; //ループ回数のカウント用 (For counting the number of loops)
	std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> SampleData;//作業位置モデル化用のサンプルデータ (Sample data for modeling work position)

};

