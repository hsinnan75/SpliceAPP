#include "structure.h"

Hollywood_info_t Run_Hollywood3(string seq)
{
	string result;
	stringstream ss;
	char cmd[MAXLINE];
	Hollywood_info_t pred;

	sprintf(cmd, "perl script/hollywood_score3.pl %s", seq.c_str()); //fprintf(stderr, "%s\n", cmd); fflush(stderr);
	shellCmd(cmd, result);

	if (result != "")
	{
		pred.bValid = true;
		ss.str(result); ss >> pred.maxent >> pred.mm >> pred.wmm;
		//fprintf(stderr, "%s\n%s:%.2f %.2f %.2f\n", seq.c_str(), result.c_str(), pred.maxent, pred.mm, pred.wmm);
	}
	return pred;
}

Hollywood_info_t Run_Hollywood5(string seq)
{
	string result;
	stringstream ss;
	char cmd[MAXLINE];
	Hollywood_info_t pred;

	sprintf(cmd, "perl script/hollywood_score5.pl %s", seq.c_str()); //fprintf(stderr, "%s\n", cmd); fflush(stderr);
	shellCmd(cmd, result);

	if (result != "")
	{
		pred.bValid = true;
		ss.str(result); ss >> pred.maxent;
		//fprintf(stderr, "%s\n%.2f\n", seq.c_str(),  pred.maxent);
	}
	return pred;
}
