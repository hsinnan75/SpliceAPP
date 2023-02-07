#include "structure.h"
#include <dirent.h>
#include <sys/types.h>
#include <unistd.h>

string DNAbase = "ACGTacgt";

typedef struct
{
	int n_snps;
	int iFinished;
	fstream file;
	FILE* out_file;
} JobInfo_t;

bool CheckSeq(string seq)
{
	if (DNAbase.find(seq[0]) != string::npos) return true;
	else return false;
}

int CountQuerySNPs(string filename)
{
	int n = 0;
	string str;
	fstream file;

	file.open(filename.c_str(), ios_base::in);
	
	while (!file.eof())
	{
		getline(file, str); 
		if (str == "" || str[0] == '#' || str.find("chr") == string::npos) continue;
		else n++;
	}
	file.close();

	return n;
}

bool CheckValidPrediction(vector<Prediction_Info_t>& PredictionVec)
{
	bool bRet = false;

	for (vector<Prediction_Info_t>::iterator iter = PredictionVec.begin(); iter != PredictionVec.end(); iter++)
	{
		if (iter->bValid)
		{
			bRet = true;
			break;
		}
	}
	return bRet;
}

void LocalMode(char* filename)
{
	FILE* outfp;
	stringstream ss;
	fstream in_file;
	bool bValidPrediction;
	Prediction_Info_t Prediction;
	int i, pos, n_snps, iFinish = 0;
	vector<Prediction_Info_t> PredictionVec;
	string str, tmp, in_filename, out_filename, chr, ref, alt;

	in_filename = filename;
	out_filename = in_filename.substr(0, in_filename.find_last_of('.')) + ".out";

	n_snps = CountQuerySNPs(in_filename); fprintf(stderr, "(%d / %d completed)\n", iFinish, n_snps);
	in_file.open(in_filename.c_str(), ios_base::in); outfp = fopen(out_filename.c_str(), "w");
	fprintf(outfp, "chr\tpos\tref\talt\tgene\tintron\tstrand\tmodel\tmt_dist\tmt_eff\teffect\n");
	while (!in_file.eof())
	{
		getline(in_file, str); if (str == "" || str[0] == '#' || str.find("chr") == string::npos) continue;
		if (str[(int)str.length() - 1] == '\r') str.resize((int)str.length() - 1);
		ss.clear(); ss.str(str); ss >> chr >> pos >> tmp >> ref >> alt;
		if (ref != alt && ref.length() == 1 && RefSeqMap[chr][pos] == ref[0] && CheckSeq(ref) && alt.length() == 1 && CheckSeq(alt))
		{
			ref = ToUpperStr(ref); alt = ToUpperStr(alt);
			//fprintf(stderr, "%s:%d %s --> %s\n", chr.c_str(), pos, ref.c_str(), alt.c_str());
			PredictionVec = SplicingEffectAnalysis(chr, pos, ref, alt);
			if (PredictionVec.size() > 0)
			{
				bValidPrediction = CheckValidPrediction(PredictionVec);
				for (i = 0; i < (int)PredictionVec.size(); i++)
				{
					Prediction = PredictionVec[i];
					if (Prediction.bValid) fprintf(outfp, "%s\t%d\t%s\t%s\t%s\t%d-%d\t%c\t%c'ss\t%d\t%.4f\t%s\n", chr.c_str(), pos, (Prediction.intron.strand ? ref.c_str() : GetReverseSeq(ref).c_str()), (Prediction.intron.strand ? alt.c_str() : GetReverseSeq(alt).c_str()), Prediction.gene.c_str(), Prediction.intron.beg, Prediction.intron.end, Prediction.intron.strand ? '+' : '-', Prediction.pred_type, Prediction.mt_dist, Prediction.pred_type == '5' ? Prediction.model0_exp : (Prediction.bNovelAG ? Prediction.model3_exp : Prediction.model4_exp), Prediction.effect ? "significant" : "neutral");
					else if (!bValidPrediction) fprintf(outfp, "%s\t%d\t%s\t%s\t%s\t%d-%d\t%c\t-\t-\t-\t-\n", chr.c_str(), pos, (Prediction.intron.strand ? ref.c_str() : GetReverseSeq(ref).c_str()), (Prediction.intron.strand ? alt.c_str() : GetReverseSeq(alt).c_str()), Prediction.gene.c_str(), Prediction.intron.beg, Prediction.intron.end, Prediction.intron.strand ? '+' : '-');
					fflush(outfp);

				}
			}
		}
		iFinish++;
	}		
}
