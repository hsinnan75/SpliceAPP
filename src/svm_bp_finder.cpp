#include "structure.h"

svm_bp_info_t Run_SVM_BP_Finder_Docker(int mt_dist, string intron_seq)
{
	string result;
	FILE* in_file = fopen("bp_query.fa", "w");
	fprintf(in_file, ">query\n%s\n", intron_seq.c_str());
	fclose(in_file);

	fprintf(stderr, "run SVM_bp.sh\n");
	shellCmd("script/SVM_bp.sh bp_query.fa Hsap 100", result);
	//printf("%s\n", result.c_str());

	svm_bp_info_t pred;
	string tmp, line, bp_seq;
	int ppt_scr, ss_dist, bp;
	stringstream parser, ss(result);
	float bp_scr, svm_scr, max_svm_scr = -100;

	pred.mut_bp_diff = -4; pred.bp_seq = ""; pred.svm_bp_ct = 0; pred.bValid = pred.bp0 = pred.bp2 = pred.bp3 = false;
	while (getline(ss, line, '\n'))
	{
		if (line.find("query") != string::npos)
		{
			//printf("%s\n\n", line.c_str());
			parser.clear(); parser.str(line); parser >> tmp >> tmp >> ss_dist >> bp_seq >> bp_scr >> tmp >> tmp >> tmp >> ppt_scr >> svm_scr;
			bp = ss_dist - mt_dist; //printf("bp = %d - %d = %d\n", ss_dist, mt_dist, bp);
			if (bp == 0)  pred.bp0 = true;
			if (bp == -2) pred.bp2 = true;
			if (bp == -3) pred.bp3 = true;
			if ((bp == 0 || bp == -2 || bp == -3) && bp > pred.mut_bp_diff) pred.mut_bp_diff = bp;
			if (svm_scr > 0)  pred.svm_bp_ct++;
			if (svm_scr > max_svm_scr)
			{
				pred.bValid = true;
				max_svm_scr = svm_scr;
				pred.ss_dist = ss_dist;
				pred.bp_seq = bp_seq;
				pred.bp_scr = bp_scr;
				pred.ppt_scr = ppt_scr;
				pred.svm_scr = svm_scr;
			}
		}
	}
	return pred;
}

