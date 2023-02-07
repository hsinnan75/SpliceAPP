#include "structure.h"

RNAcofold_info_t Run_RNAcofold(string bp_seq)
{
	FILE* fp = NULL;
	string result, tmp;
	RNAcofold_info_t pred;

	fp = fopen("rna_fold_seq.in", "w");
	fprintf(fp, ">query\n%s&GTGTAGTA\n(((((x(((&))))))))\n", bp_seq.c_str()); fclose(fp);
	shellCmd("/usr/local/bin/RNAcofold --noPS -C --output-format=D < rna_fold_seq.in", result);
	stringstream ss(result);
	getline(ss, tmp, '\n'); getline(ss, pred.str, '\n');
	tmp = pred.str.substr(pred.str.find_first_of(' ') + 1); pred.mfe = atof(tmp.c_str());

	return pred;
}

float RNAfold(string query_seq)
{
	int p1, p2;
	float ret = 0;
	FILE* fp = NULL;
	string result, str, tmp;

	fp = fopen("rna_fold_seq.in", "w"); fprintf(fp, ">query\n%s\n", query_seq.c_str()); fclose(fp);
	shellCmd("/usr/local/bin/RNAfold --noPS < rna_fold_seq.in", result);
	stringstream ss(result);
	getline(ss, str, '\n'); getline(ss, str, '\n'); getline(ss, str, '\n');
	p1 = str.find_last_of('(') + 1; p2 = str.find_last_of(')');
	tmp = str.substr(p1, p2 - p1); ret = atof(tmp.c_str());
	//printf("ret=%f\n", ret);

	return ret;
}

string RNAfold_str(string query_seq)
{
	int len;
	FILE* fp = NULL;
	string result, str, ret_str;

	fp = fopen("rna_fold_seq.in", "w"); fprintf(fp, ">query\n%s\n", query_seq.c_str()); fclose(fp);
	shellCmd("/usr/local/bin/RNAfold --noPS < rna_fold_seq.in", result);
	stringstream ss(result);
	getline(ss, str, '\n'); getline(ss, str, '\n'); getline(ss, str, '\n');
	
	len = (int)query_seq.length(); ret_str = str.substr(0, len);

	return ret_str;
}
