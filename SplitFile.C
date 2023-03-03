#include "ana_const.h"
void SplitFile(){
	 ROOT::EnableImplicitMT();
	//source data
	auto df1 = ROOT::RDataFrame("result_mc", "/Users/forhadhossain/SeaQuest/Unfold/Unfolding/Unfolding1D/data/DY_DUMP_4pi_GMC_Jan08_LD2.root"); //Train data for response matrix
	auto df2 = ROOT::RDataFrame("result_mc", "/Users/forhadhossain/SeaQuest/Unfold/Unfolding/Unfolding1D/Y_4pi_GMC_Jan03_LD2_Dataset1and3.root"); // combined test data1
	auto df3 = ROOT::RDataFrame("result_mc", "/Users/forhadhossain/SeaQuest/Unfold/Unfolding/Unfolding1D/Y_4pi_GMC_Jan03_LD2_Dataset2and4.root"); //combined test data2
	auto df4 = ROOT::RDataFrame("result_mc", "/Users/forhadhossain/SeaQuest/Unfold/Unfolding/Unfolding1D/data/GMC_Nov19_2022_LD2_Acc_DY_67_5to7GeV.root"); // combined test data3
	auto df5 = ROOT::RDataFrame("result_mc", "/Users/forhadhossain/SeaQuest/Unfold/Unfolding/Unfolding1D/combined.root"); // Independent 4pi (from NONE target); for acc factor

	df1.Filter(cutTrue).Snapshot("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Train.root");  //for response matrix
	df5.Filter(cutTrue).Snapshot("result_mc", "data/DY_4pi_GMC_4piAccfactor.root");  // denominator 4pi accceptance factor
	df2.Filter(cutRecoMC).Snapshot("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Test1.root");
	df3.Filter(cutRecoMC).Snapshot("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Test2.root");
	df4.Filter(cutRecoMC).Snapshot("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Test3.root");
}
