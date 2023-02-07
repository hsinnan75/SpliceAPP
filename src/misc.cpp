#include "structure.h"

bool HGMD_Query(string chr, int pos)
{
	map<pair<string, int>, bool>::iterator iter;

	if ((iter = HGMDmap.find(make_pair(chr, pos))) != HGMDmap.end()) return true;
	else return false;

}

bool ClinVar_Query(string chr, int pos)
{
	map<int, int>::iterator iter;
	if((iter = ClinVarMap[chr].lower_bound(pos))!= ClinVarMap[chr].end() && pos >= iter->second && pos <= iter->first)
	{
		//printf("match!!\n");
		return true;
	}
	else return false;
}

bool dbSNPmap_Query(string chr, int pos, bool strand)
{
	map<pair<string, int>, bool>::iterator iter;

	if ((iter = dbSNPmap.find(make_pair(chr, pos))) == dbSNPmap.end() || iter->second != strand) return false;
	else return true;
}

bool AltSpliceRegion_Query(string chr, int pos, bool strand)
{
	bool bRet = false;
	map<int, AltSpliceInfo_t>::iterator iter;

	for (iter = AltSpliceMap[chr].lower_bound(pos); iter != AltSpliceMap[chr].end(); iter++)
	{
		if (strand == iter->second.strand && pos >= iter->second.beg && pos <= iter->second.end)
		{
			//printf("match ucsc_alt_db %s:%d-%d %c\n", chr.c_str(), iter->second.beg, iter->second.end, iter->second.strand? '+':'-');
			bRet = true;
			break;
		}
	}
	return bRet;
}

int Count_Adenine(string seq)
{
	int i, len, n = 0;

	for (len = (int)seq.length(), i = 0; i < len; i++)
	{
		if (toupper(seq[i]) == 'A') n++;
	}
	return n;
}

map<string, int> MotifScan(string seq)
{
	int count;
	fstream file;
	char cmd[1024];
	stringstream ss;
	string str, motif_id;
	map<string, int> MotifMap;
	map<string, string>::iterator iter;

	sprintf(cmd, "Rscript script/countpwm.R %s > /dev/null 2>&1", seq.c_str()); fprintf(stderr, "run command:%s\n", cmd);
	if (system(cmd) == 0)
	{
		file.open("pwm_motif.tsv", ios_base::in);
		while (!file.eof())
		{
			getline(file, str); if (str == "") break;
			ss.clear(); ss.str(str); ss >> motif_id >> count;
			iter = AttractMap.find(motif_id);
			if (iter != AttractMap.end() && MotifMap[iter->second] < count) MotifMap[iter->second] = count;
		}
		file.close();
	}
	return MotifMap;
}

float AvgPhyloP100_5(string chr, intron_t intron)
{
	int pos;
	float avg = 0;

	for (pos = intron.beg - 3; pos < intron.beg + 6; pos++) avg += PhyloP100Map[chr][pos];

	return avg /= 9;
}

float AvgPhyloP100_5(string chr, int beg, int end)
{
	int pos;
	float avg = 0;

	for (pos = beg; pos <= end; pos++) avg += PhyloP100Map[chr][pos];
	//fprintf(stderr, "avg = %f / %d = %f\n", avg, (end - beg + 1), avg/(end - beg + 1));
	return (avg / (end - beg + 1));
}
