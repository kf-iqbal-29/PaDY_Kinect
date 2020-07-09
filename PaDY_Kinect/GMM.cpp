#include "GMM.h"
#include <iostream>
#include <direct.h>
#include <sstream>
#include <fstream>
#include <string>

//#define DEBUG //�f�o�b�N�p (for debugging)

//////* �ÓI�����o�ϐ��̒�` *////// (Defining static member variables)
const double GMM::INI_SIGMA = 0.02;//�������U (Initial variance)
const double GMM::ALPHA = 0.8;//EM�A���S���Y���̃p�����[�^�X�V�p (For updating the parameters of the EM algorithm)
const double GMM::INITIAL_POSITION[TASK_NUM + 1][2] = //������ƈʒu (Initial working position)
//{
//	0.6, -0.6,
//	1.9, 0.8,
//	0.9, 1.4,
//	0.6, 0.8,
//};
{
	0.6, -0.6,
		1.9, 0.7,
		0.9, 1.2,
		0.6, 0.6,
};

//�R���X�g���N�^ (constructor)
GMM::GMM()
{
	//��ƈʒu���f���̏����ݒ� (Initial setting of work position model)
	Model.resize(TASK_NUM + 1);
	for (int i = 0; i < (TASK_NUM + 1); ++i){
		//�����W�� (Mixing coefficient)
		Model[i].c = 1.0 / static_cast<double>(TASK_NUM + 1);
		//���σx�N�g�� (Mean vector)
		Model[i].mu(0) = INITIAL_POSITION[i][0];
		Model[i].mu(1) = INITIAL_POSITION[i][1];
		//�����U�s�� (Covariance matrix)
		Model[i].sigma(0, 0) = INI_SIGMA;
		Model[i].sigma(0, 1) = 0.0;
		Model[i].sigma(1, 0) = 0.0;
		Model[i].sigma(1, 1) = INI_SIGMA;
	}

}

//�f�X�g���N�^ (Destructor)
GMM::~GMM()
{
	SampleData.clear();//�T���v���f�[�^�p�o�b�t�@�̉��
	Model.clear();//��ƈʒu���f���̉��
}


//2�����̍����K�E�X���z�̖ޓx��Ԃ� (Give the likelihood of the two-dimensional Gaussian mixture distribution)
double GMM::CalcMixtureGaussian2d(const Eigen::Vector2d& x, const std::vector<gmm_t, Eigen::aligned_allocator<gmm_t>>& gmm)
{
	double val;

	val = 0;
	for (int i = 0; i < static_cast<int>(gmm.size()); ++i){
		val += gmm[i].c
			* StdCalc::Calc2dGaussian(x, gmm[i].mu, gmm[i].sigma);
	}
#ifdef DEBUG
	std::cout << "GMM Likelihood : " << val << std::endl;
#endif// DEBUG
	return val;
}

//GMM���̂������̃K�E�X���z�ɑ�����m����Ԃ� (Give the probability of belonging to a particular Gaussian distribution in GMM)
double GMM::CalcGaussianProb(const Eigen::Vector2d& x, int No)
{
	return (Model[No].c * StdCalc::Calc2dGaussian(x, Model[No].mu, Model[No].sigma) / CalcMixtureGaussian2d(x, Model));
}


//EM�A���S���Y���ɂ��GMM�̍X�V (GMM update by EM algorithm)
void GMM::EMAlgorithm()
{
	double sum, mix_yudo;
	double eta;
	std::vector<std::vector<double>> gamma;//���S�� (Burden rate)
	std::vector<gmm_t, Eigen::aligned_allocator<gmm_t>> temp, sum_temp;
	gmm_t zero;//vector��0�������p (for 0 initialization of vector)
	gmm_t new_clust;
	Eigen::Vector2d diff;
	Eigen::Matrix2d mat_diff;

	
#ifdef DEBUG
		std::cout << "Start EM, Sample Num : " << static_cast<int>(SampleData[task].size()) << std::endl;
#endif// DEBUG

	if (static_cast<int>(SampleData.size()) < 10) return;//�f�[�^�������Ȃ�������X�L�b�v (Skip if there is too little data)


	/* gamma�̏����� (gamma initialization) */ 
	gamma.resize(static_cast<int>(Model.size()));
	for (int i = 0; i < static_cast<int>(Model.size()); ++i){
		gamma[i].clear();
		gamma[i].resize(static_cast<int>(SampleData.size()));
	}
	temp = Model;//�ΏۂƂ��鍬���K�E�X���z���Z�b�g (Set the mixture Gaussian distribution of interest)

	/* E Step (���S���̌v�Z) (Burden rate calculation)*/
	for (int i = 0; i < static_cast<int>(SampleData.size()); ++i){
		mix_yudo = CalcMixtureGaussian2d(SampleData[i], temp);
		for (int j = 0; j < static_cast<int>(Model.size()); ++j){
			gamma[j][i] = temp[j].c * StdCalc::Calc2dGaussian(SampleData[i], temp[j].mu, temp[j].sigma) / mix_yudo;
		}
	}
	/* M Step (Stepwise���Ń��f���p�����[�^�̍X�V) (Update model parameters with Stepwise formula)*/
	//������ (Initialization)
	memset(&zero, 0, sizeof(gmm_t));
	sum_temp.clear();
	sum_temp.resize(static_cast<int>(Model.size()), zero);
	//�e�ϑ��f�[�^���琄��l���v�Z (Calculate estimated value from each observation data)
	for (int j = 0; j < static_cast<int>(Model.size()); ++j){
		sum = 0;
		for (int i = 0; i < static_cast<int>(SampleData.size()); ++i){
			sum_temp[j].mu.array() += gamma[j][i] * SampleData[i].array();
			diff = SampleData[i] - temp[j].mu;
			mat_diff = diff * diff.transpose();
			sum_temp[j].sigma.array() += gamma[j][i] * mat_diff.array();
			sum += gamma[j][i];
		}

		//�e�p�����[�^���X�V (Update each parameter)
		temp[j].mu.array() = sum_temp[j].mu.array() / sum;
		temp[j].sigma.array() = sum_temp[j].sigma.array() / sum;
		temp[j].c = sum / static_cast<double>(SampleData.size());
	} 
	
	/* GMM�̃p�����[�^�X�V (GMM parameter update) */
	eta = pow(static_cast<double>(counter + 2), (-ALPHA));
	for (int j = 0; j < static_cast<int>(Model.size()); ++j){
		Model[j].c = (1.0 - eta) * Model[j].c + eta * temp[j].c;
		Model[j].mu.array() = (1.0 - eta) * Model[j].mu.array() + eta * temp[j].mu.array();
		Model[j].sigma.array() = (1.0 - eta) * Model[j].sigma.array() + eta * temp[j].sigma.array();
	}

#ifdef DEBUG
	std::cout << "c : " << std::endl;
	for (int i = 0; i < static_cast<int>(Model.size()); ++i){
		std::cout << Model[i].c << std::endl;
	}
	std::cout << "mu : " << std::endl;
	for (int i = 0; i < static_cast<int>(Model.size()); ++i){
		std::cout << Model[i].mu << std::endl;
	}
	std::cout << "sigma : " << std::endl;
	for (int i = 0; i < static_cast<int>(Model.size()); ++i){
		std::cout << Model[i].sigma << std::endl;
	}
	std::cout << "Mixture Num : " << static_cast<int>(Model.size()) << std::endl;
	getchar();
#endif// DEBUG

	counter++;//�J�E���^���C���N�������g (Increment the counter)
	InitializeSampleData();//�T���v���f�[�^�������� (Initialize sample data)

}

//���f�����̓ǂݍ���
int GMM::ReadModelInfo()
{
	int no;

	//�ꎞ�ۑ��p�X�g���[����p��
	std::string str, file_name;
	std::stringstream ss;
	std::ifstream ifs;

	//��ƈʒu���f���̏���ǂݍ���
	//csv�t�@�C����W�J
	std::ostringstream oss;
	oss << "Input/GMM/model.csv";
	file_name = oss.str();
	ifs.clear();
	ifs.open(file_name);
	if (!ifs){
		std::cout << "Error: " << file_name << " cannot open" << std::endl;
		return -1;
	}
	//���s�R�[�h�܂œǂݍ���
	getline(ifs.seekg(0, std::ios_base::cur), str, '\n');
	//stringstream�ɓǂ݂�����string�𗬂�
	ss.str(str);
	//�T���v���ԍ���ǂݍ���
	ss >> no;
	//stringstream���ȉ��̂Q�s�̃R�[�h�ŃN���A����
	ss.str("");
	ss.clear(std::stringstream::goodbit);
	//���������������[�v
	for (int i = 0; i < (GMM::TASK_NUM + 1); ++i){
		/* �����W�� */
		//���s�R�[�h�܂œǂݍ���
		getline(ifs.seekg(0, std::ios_base::cur), str, '\n');
		//stringstream�ɓǂ݂�����string�𗬂�
		ss.str(str);
		//�����W����ǂݍ���
		ss >> Model[i].c;
		//stringstream���ȉ��̂Q�s�̃R�[�h�ŃN���A����
		ss.str("");
		ss.clear(std::stringstream::goodbit);
		/* ���σx�N�g�� */
		for (int col = 0; col < 2; ++col){
			//1��ǂݍ���
			getline(ifs.seekg(0, std::ios_base::cur), str, ',');
			//stringstream�ɓǂ݂�����string�𗬂�
			ss.str(str);
			//���σx�N�g���̓ǂݍ���
			ss >> Model[i].mu(col);
			//stringstream���ȉ��̂Q�s�̃R�[�h�ŃN���A����
			ss.str("");
			ss.clear(std::stringstream::goodbit);
		}
		//���s�R�[�h�܂œǂݍ���
		getline(ifs.seekg(0, std::ios_base::cur), str, '\n');
		/* �����U�s�� */
		for (int col = 0; col < 2; ++col){
			for (int row = 0; row < 2; ++row){
				//1��ǂݍ���
				getline(ifs.seekg(0, std::ios_base::cur), str, ',');
				//stringstream�ɓǂ݂�����string�𗬂�
				ss.str(str);
				//�����U�s��̓ǂݍ���
				ss >> Model[i].sigma(col, row);
				//stringstream���ȉ��̂Q�s�̃R�[�h�ŃN���A����
				ss.str("");
				ss.clear(std::stringstream::goodbit);
			}
		}
		//���s�R�[�h�܂œǂݍ���
		getline(ifs.seekg(0, std::ios_base::cur), str, '\n');//�ǂݍ��݂����܂������Ȃ��ꍇ�̓R�����g�A�E�g
	}
	ifs.close();

	return no;
}

//���f�����̏�������
void GMM::WriteModelInfo(int cur_sample)
{

	//��ƈʒu���f���̏��ۑ��p�t�H���_�쐬
	_mkdir("LogFiles");
	_mkdir("LogFiles/GMM");
	//�ꎞ�ۑ��p�X�g���[����p��
	std::ofstream ofs;
	std::string file_name;
	//csv�t�@�C����W�J
	std::ostringstream oss;
	oss << "LogFiles/GMM/model.csv";
	file_name = oss.str();
	ofs.clear();
	ofs.open(file_name);
	//���݂̃T���v���ԍ�����������
	ofs << cur_sample << std::endl;
	//�����K�E�X���z�̊e�p�����[�^����������
	for (int i = 0; i < (GMM::TASK_NUM + 1); ++i){
		//�����W��
		ofs << Model[i].c << std::endl;
		//���σx�N�g��
		for (int col = 0; col < 2; ++col){
			ofs << Model[i].mu(col) << ",";
		}
		ofs << std::endl;
		//�����U�s��
		for (int col = 0; col < 2; ++col){
			for (int row = 0; row < 2; ++row){
				ofs << Model[i].sigma(col, row) << ",";
			}
		}
		ofs << std::endl;
	}
	ofs.close();

}