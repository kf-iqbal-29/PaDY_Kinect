#include "GMM.h"
#include <iostream>
#include <direct.h>
#include <sstream>
#include <fstream>
#include <string>

//#define DEBUG //デバック用 (for debugging)

//////* 静的メンバ変数の定義 *////// (Defining static member variables)
const double GMM::INI_SIGMA = 0.02;//初期分散 (Initial variance)
const double GMM::ALPHA = 0.8;//EMアルゴリズムのパラメータ更新用 (For updating the parameters of the EM algorithm)
const double GMM::INITIAL_POSITION[TASK_NUM + 1][2] = //初期作業位置 (Initial working position)
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

//コンストラクタ (constructor)
GMM::GMM()
{
	//作業位置モデルの初期設定 (Initial setting of work position model)
	Model.resize(TASK_NUM + 1);
	for (int i = 0; i < (TASK_NUM + 1); ++i){
		//混合係数 (Mixing coefficient)
		Model[i].c = 1.0 / static_cast<double>(TASK_NUM + 1);
		//平均ベクトル (Mean vector)
		Model[i].mu(0) = INITIAL_POSITION[i][0];
		Model[i].mu(1) = INITIAL_POSITION[i][1];
		//共分散行列 (Covariance matrix)
		Model[i].sigma(0, 0) = INI_SIGMA;
		Model[i].sigma(0, 1) = 0.0;
		Model[i].sigma(1, 0) = 0.0;
		Model[i].sigma(1, 1) = INI_SIGMA;
	}

}

//デストラクタ (Destructor)
GMM::~GMM()
{
	SampleData.clear();//サンプルデータ用バッファの解放
	Model.clear();//作業位置モデルの解放
}


//2次元の混合ガウス分布の尤度を返す (Give the likelihood of the two-dimensional Gaussian mixture distribution)
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

//GMM内のある特定のガウス分布に属する確率を返す (Give the probability of belonging to a particular Gaussian distribution in GMM)
double GMM::CalcGaussianProb(const Eigen::Vector2d& x, int No)
{
	return (Model[No].c * StdCalc::Calc2dGaussian(x, Model[No].mu, Model[No].sigma) / CalcMixtureGaussian2d(x, Model));
}


//EMアルゴリズムによるGMMの更新 (GMM update by EM algorithm)
void GMM::EMAlgorithm()
{
	double sum, mix_yudo;
	double eta;
	std::vector<std::vector<double>> gamma;//負担率 (Burden rate)
	std::vector<gmm_t, Eigen::aligned_allocator<gmm_t>> temp, sum_temp;
	gmm_t zero;//vectorの0初期化用 (for 0 initialization of vector)
	gmm_t new_clust;
	Eigen::Vector2d diff;
	Eigen::Matrix2d mat_diff;

	
#ifdef DEBUG
		std::cout << "Start EM, Sample Num : " << static_cast<int>(SampleData[task].size()) << std::endl;
#endif// DEBUG

	if (static_cast<int>(SampleData.size()) < 10) return;//データ数が少なすぎたらスキップ (Skip if there is too little data)


	/* gammaの初期化 (gamma initialization) */ 
	gamma.resize(static_cast<int>(Model.size()));
	for (int i = 0; i < static_cast<int>(Model.size()); ++i){
		gamma[i].clear();
		gamma[i].resize(static_cast<int>(SampleData.size()));
	}
	temp = Model;//対象とする混合ガウス分布をセット (Set the mixture Gaussian distribution of interest)

	/* E Step (負担率の計算) (Burden rate calculation)*/
	for (int i = 0; i < static_cast<int>(SampleData.size()); ++i){
		mix_yudo = CalcMixtureGaussian2d(SampleData[i], temp);
		for (int j = 0; j < static_cast<int>(Model.size()); ++j){
			gamma[j][i] = temp[j].c * StdCalc::Calc2dGaussian(SampleData[i], temp[j].mu, temp[j].sigma) / mix_yudo;
		}
	}
	/* M Step (Stepwise式でモデルパラメータの更新) (Update model parameters with Stepwise formula)*/
	//初期化 (Initialization)
	memset(&zero, 0, sizeof(gmm_t));
	sum_temp.clear();
	sum_temp.resize(static_cast<int>(Model.size()), zero);
	//各観測データから推定値を計算 (Calculate estimated value from each observation data)
	for (int j = 0; j < static_cast<int>(Model.size()); ++j){
		sum = 0;
		for (int i = 0; i < static_cast<int>(SampleData.size()); ++i){
			sum_temp[j].mu.array() += gamma[j][i] * SampleData[i].array();
			diff = SampleData[i] - temp[j].mu;
			mat_diff = diff * diff.transpose();
			sum_temp[j].sigma.array() += gamma[j][i] * mat_diff.array();
			sum += gamma[j][i];
		}

		//各パラメータを更新 (Update each parameter)
		temp[j].mu.array() = sum_temp[j].mu.array() / sum;
		temp[j].sigma.array() = sum_temp[j].sigma.array() / sum;
		temp[j].c = sum / static_cast<double>(SampleData.size());
	} 
	
	/* GMMのパラメータ更新 (GMM parameter update) */
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

	counter++;//カウンタをインクリメント (Increment the counter)
	InitializeSampleData();//サンプルデータを初期化 (Initialize sample data)

}

//モデル情報の読み込み
int GMM::ReadModelInfo()
{
	int no;

	//一時保存用ストリームを用意
	std::string str, file_name;
	std::stringstream ss;
	std::ifstream ifs;

	//作業位置モデルの情報を読み込み
	//csvファイルを展開
	std::ostringstream oss;
	oss << "Input/GMM/model.csv";
	file_name = oss.str();
	ifs.clear();
	ifs.open(file_name);
	if (!ifs){
		std::cout << "Error: " << file_name << " cannot open" << std::endl;
		return -1;
	}
	//改行コードまで読み込む
	getline(ifs.seekg(0, std::ios_base::cur), str, '\n');
	//stringstreamに読みだしたstringを流す
	ss.str(str);
	//サンプル番号を読み込み
	ss >> no;
	//stringstreamを以下の２行のコードでクリアする
	ss.str("");
	ss.clear(std::stringstream::goodbit);
	//混合数分だけループ
	for (int i = 0; i < (GMM::TASK_NUM + 1); ++i){
		/* 混合係数 */
		//改行コードまで読み込む
		getline(ifs.seekg(0, std::ios_base::cur), str, '\n');
		//stringstreamに読みだしたstringを流す
		ss.str(str);
		//混合係数を読み込み
		ss >> Model[i].c;
		//stringstreamを以下の２行のコードでクリアする
		ss.str("");
		ss.clear(std::stringstream::goodbit);
		/* 平均ベクトル */
		for (int col = 0; col < 2; ++col){
			//1列読み込み
			getline(ifs.seekg(0, std::ios_base::cur), str, ',');
			//stringstreamに読みだしたstringを流す
			ss.str(str);
			//平均ベクトルの読み込み
			ss >> Model[i].mu(col);
			//stringstreamを以下の２行のコードでクリアする
			ss.str("");
			ss.clear(std::stringstream::goodbit);
		}
		//改行コードまで読み込む
		getline(ifs.seekg(0, std::ios_base::cur), str, '\n');
		/* 共分散行列 */
		for (int col = 0; col < 2; ++col){
			for (int row = 0; row < 2; ++row){
				//1列読み込み
				getline(ifs.seekg(0, std::ios_base::cur), str, ',');
				//stringstreamに読みだしたstringを流す
				ss.str(str);
				//共分散行列の読み込み
				ss >> Model[i].sigma(col, row);
				//stringstreamを以下の２行のコードでクリアする
				ss.str("");
				ss.clear(std::stringstream::goodbit);
			}
		}
		//改行コードまで読み込む
		getline(ifs.seekg(0, std::ios_base::cur), str, '\n');//読み込みがうまくいかない場合はコメントアウト
	}
	ifs.close();

	return no;
}

//モデル情報の書き込み
void GMM::WriteModelInfo(int cur_sample)
{

	//作業位置モデルの情報保存用フォルダ作成
	_mkdir("LogFiles");
	_mkdir("LogFiles/GMM");
	//一時保存用ストリームを用意
	std::ofstream ofs;
	std::string file_name;
	//csvファイルを展開
	std::ostringstream oss;
	oss << "LogFiles/GMM/model.csv";
	file_name = oss.str();
	ofs.clear();
	ofs.open(file_name);
	//現在のサンプル番号を書き込み
	ofs << cur_sample << std::endl;
	//混合ガウス分布の各パラメータを書き込み
	for (int i = 0; i < (GMM::TASK_NUM + 1); ++i){
		//混合係数
		ofs << Model[i].c << std::endl;
		//平均ベクトル
		for (int col = 0; col < 2; ++col){
			ofs << Model[i].mu(col) << ",";
		}
		ofs << std::endl;
		//共分散行列
		for (int col = 0; col < 2; ++col){
			for (int row = 0; row < 2; ++row){
				ofs << Model[i].sigma(col, row) << ",";
			}
		}
		ofs << std::endl;
	}
	ofs.close();

}