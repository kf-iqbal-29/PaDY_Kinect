/* ��Ǝ҈ʒu���f�����ƍ�Ɛi�x����̂��߂�GMM�p�N���X (Class for GMM for worker position modeling and work progress estimation)*/

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

	GMM();//�R���X�g���N�^ (Constructor)
	~GMM();//�f�X�g���N�^ (Destructor)

	static const int TASK_NUM = 3;//GMM�̃N���X�^�̐� (Number of GMM clusters) // It was 3 in Kanazawa's program
	static const double INITIAL_POSITION[TASK_NUM + 1][2];//������ƈʒu (Initial working position)

	std::vector<gmm_t, Eigen::aligned_allocator<gmm_t>> Model;//���f���p��{�I�u�W�F�N�g (Basic object for model)

	void InitializeSampleData(){ SampleData.resize(0); }//�T���v���f�[�^�̏����� (Initialization of sample data)
	void SetSampleData(const Eigen::Vector2d& x){ SampleData.push_back(x);}//�T���v���f�[�^�̃Z�b�g (Set of sample data)

	void EMAlgorithm();//EM�A���S���Y���ɂ��GMM�̍X�V (GMM update by EM algorithm)
	double CalcMixtureGaussian2d(const Eigen::Vector2d& x, const std::vector<gmm_t, Eigen::aligned_allocator<gmm_t>>& gmm);//2�����̍����K�E�X���z�̖ޓx��Ԃ� (Give the likelihood of the two-dimensional Gaussian mixture distribution)
	double CalcGaussianProb(const Eigen::Vector2d& x, int No);
	int ReadModelInfo();//���f�����̓ǂݍ��� (Load model information)
	void WriteModelInfo(int cur_sample);//���f�����̏������� (Writing model information)

private:

	static const double INI_SIGMA;//�������U (Initial variance)
	static const double ALPHA;//EM�A���S���Y���̃p�����[�^�X�V�p (For updating the parameters of the EM algorithm)

	int counter = 0; //���[�v�񐔂̃J�E���g�p (For counting the number of loops)
	std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> SampleData;//��ƈʒu���f�����p�̃T���v���f�[�^ (Sample data for modeling work position)

};

