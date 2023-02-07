#include "structure.h"

char boundary[41];

string ToUpperStr(string seq)
{
	string ret_seq;

	for (string::iterator iter = seq.begin(); iter != seq.end(); iter++) ret_seq.push_back(toupper(*iter));

	return ret_seq;
}

string GetReverseSeq(string seq)
{
	string rseq;
	int i, j, len;

	len = (int)seq.length(); rseq.resize(len);
	for (i = 0, j = len - 1; i < len; i++, j--)
	{
		switch (seq[i])
		{
		case 'A':
		case 'a': rseq[j] = 'T'; break;
		case 'C':
		case 'c': rseq[j] = 'G'; break;
		case 'G':
		case 'g': rseq[j] = 'C'; break;
		case 'T':
		case 't': rseq[j] = 'A'; break;
		default: rseq[j] = 'N';
		}
	}
	return rseq;
}

// void InitWebKitFormBoundary(char* boundary)
// {
// 	int i, j;
// 	char random_str[17];

// 	srand((unsigned int)time(NULL));
// 	for (i = 0; i < 16; i++)
// 	{
// 		j = rand() % 62;
// 		if (j < 10) random_str[i] = j + 48; // 0-9
// 		else if (j < 36) random_str[i] = j + 55; // A-Z
// 		else random_str[i] = j + 61;
// 	}
// 	random_str[16] = '\0';
// 	sprintf(boundary, "------WebKitFormBoundary%s", random_str); boundary[40] = '\0';
// }

int Count_AG_Occurrences(string str)
{
	int i, j, len, ag = 0;

	len = (int)str.length();
	for (i = 0, j = 1; j < len; i++, j++)
	{
		if (toupper(str[i]) == 'A' && toupper(str[j]) == 'G') ag++;
	}
	return ag;
}

int FindNewAG(int pos, string wt3seq, string mt3seq)
{
	int i, j, new3ss = -1;
	
	for (i = 0, j = 1; j < (int)mt3seq.length(); i++, j++)
	{
		if (toupper(mt3seq[i]) == 'A' && toupper(mt3seq[j]) == 'G')
		{
			if (j < (int)wt3seq.length() && (toupper(wt3seq[i]) != 'A' || toupper(wt3seq[j]) != 'G'))
			{
				new3ss = obs_intron_len + left_primer_len - pos + j;
				break;
			}
		}
	}
	return new3ss;
}

float GC_count(string& seq)
{
	char c;
	int i, len, count = 0;

	for (len = seq.length(), i = 0; i < len; i++)
	{
		c = toupper(seq[i]);
		if (c == 'C' || c == 'G') count++;
	}
	//printf("GC%%: %s\n%d/%d = %f\n", seq.c_str(), count, len, 1.0 * count / len);
	return 1. * count / len;
}

GC_info_t Get_GC_Info(int x3ss, string seq)
{
	string sub_seq;
	GC_info_t gc_info;

	gc_info.gc = GC_count(seq);
	sub_seq = seq.substr(x3ss, obs_exon_len); gc_info.ex_gc = GC_count(sub_seq);
	sub_seq = seq.substr(left_primer_len, x3ss - left_primer_len); gc_info.in_gc = GC_count(sub_seq);

	return gc_info;
}

GC_info_t Get_GC_Info(int new3ss, int mt_x3ss, string seq)
{
	string sub_seq;
	GC_info_t gc_info;

	sub_seq = seq.substr(new3ss); gc_info.ex_gc = GC_count(sub_seq); //printf("novel_exon  : %s\n", sub_seq.c_str());
	sub_seq = seq.substr(new3ss, mt_x3ss - new3ss); gc_info.short_ex_gc = GC_count(sub_seq); //printf("novel_short : %s\n", sub_seq.c_str());
	sub_seq = seq.substr(0, new3ss); gc_info.in_gc = GC_count(sub_seq);//printf("novel_intron: %s\n", sub_seq.c_str());
	
	return gc_info;
}

FreeEnergy_Info_t Get_FreeEnergy_Info(int x3ss, string seq)
{
	string sub_seq;
	FreeEnergy_Info_t dG_info;

	dG_info.dG = RNAfold(seq);
	sub_seq = seq.substr(x3ss); dG_info.ex_dG = RNAfold(sub_seq); //printf("ex.dG = %f\n", dG_info.ex_dG);
	sub_seq = seq.substr(0, x3ss); dG_info.in_dG = RNAfold(sub_seq); //printf("in.dG = %f\n", wt_indG);

	return dG_info;
}

float Get_FreeEnergy_Info(string seq)
{
	float dG;

	dG = RNAfold(seq);

	return dG;
}

int CountOpenness(string str)
{
	int i, len, openness = 0;

	for (len = (int)str.length(), i = 0; i < len; i++)
	{
		if (str[i] == '.') openness++;
	}
	return openness;
}

bool shellCmd(const string& cmd, string& result)
{
	char buffer[512];
	
	result = "";
	// Open pipe to file
	FILE* pipe = popen(cmd.c_str(), "r");
	if (!pipe) {
		return false;
	}
	// read till end of process:
	while (!feof(pipe)) {
		// use buffer to read and add to result
		if (fgets(buffer, sizeof(buffer), pipe) != NULL)
			result += buffer;
	}
	pclose(pipe);
	return true;
}
