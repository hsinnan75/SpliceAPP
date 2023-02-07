#include "structure.h"

#define MaxMutDistance 50
#define X5_Threshold 0.13
#define AG_Threshold 0.1837238
#define NonAG_Threshold 0.07784796

bool debug = false;

intron_t FindPrevIntron(int idx, vector<intron_t>& intron_vec)
{
	int i; 
	intron_t intron;

	if (intron_vec[idx].strand)
	{
		for (i = idx - 1; i >= 0; i--)
		{
			if (intron_vec[i].end < intron_vec[idx].beg) break;
		}
	}
	else
	{
		for (i = idx + 1; i < (int)intron_vec.size(); i++)
		{
			if (intron_vec[i].beg > intron_vec[idx].end) break;
		}
	}
	if (i >= 0) return intron_vec[i];
	else
	{
		intron.beg = intron.end = 0;
		return intron;
	}
}

string GetObservedSeq5(string chr, int pos, int intron_idx, vector<intron_t>& intron_vec, int& exon_beg, int& exon_end, int& intron1_len, int& intron2_len)
{
	string seq;
	int beg, end;
	intron_t PrevIntron, NextIntron;

	if (intron_vec[intron_idx].strand)
	{
		if (intron_idx > 0)
		{
			if ((PrevIntron = FindPrevIntron(intron_idx, intron_vec)).beg > 0)
			{
				intron1_len = PrevIntron.end - PrevIntron.beg + 1;
				intron2_len = intron_vec[intron_idx].end - intron_vec[intron_idx].beg + 1;
				beg = PrevIntron.end - 39;
				end = intron_vec[intron_idx].beg + 30;
				seq = RefSeqMap[chr].substr(beg, end - beg);
				exon_beg = PrevIntron.end + 1;
				exon_end = intron_vec[intron_idx].beg - 1;
				//fprintf(stderr, "intron[%d-%d]-exon[%d-%d]-intron[%d-%d]\n", PrevIntron.beg, PrevIntron.end, exon_beg, exon_end, intron_vec[intron_idx].beg, intron_vec[intron_idx].end);
				//fprintf(stderr, "wt_seq [%d - %d] = %d\n%s\n", beg + 1, end, (int)seq.length(), seq.c_str());
			}
		}
	}
	else
	{
		if(intron_idx < (int)intron_vec.size() - 1)
		{
			if ((NextIntron = FindPrevIntron(intron_idx, intron_vec)).beg > 0)
			{
				intron1_len = NextIntron.end - NextIntron.beg + 1;
				intron2_len = intron_vec[intron_idx].end - intron_vec[intron_idx].beg + 1;
				beg = intron_vec[intron_idx].end - 29;
				end = NextIntron.beg + 40;
				seq = RefSeqMap[chr].substr(beg, end - beg);
				seq = GetReverseSeq(seq);
				exon_beg = intron_vec[intron_idx].end + 1;
				exon_end = NextIntron.beg - 1;
				//fprintf(stderr, "intron[%d-%d]-exon[%d-%d]-intron[%d-%d]\n", intron_vec[intron_idx].beg, intron_vec[intron_idx].end, exon_beg, exon_end, NextIntron.beg, NextIntron.end);
				//fprintf(stderr, "wt_seq [%d - %d] = %d\n%s\n", beg + 1, end, seq.length(), seq.c_str());
			}
		}
	}
	return seq;
}

string GetObservedSeq3(string chr, int pos, intron_t intron)
{
	string seq;

	if (intron.strand)
	{
		if (intron.end - pos <= obs_intron_len)
		{
			seq = RefSeqMap[chr].substr(intron.end - obs_intron_len + 1, obs_intron_len);
			seq += RefSeqMap[chr].substr(intron.end + 1, obs_exon_len);
		}
	}
	else
	{
		if (pos - intron.beg <= obs_intron_len)
		{
			seq = RefSeqMap[chr].substr(intron.beg - obs_exon_len, obs_exon_len);
			seq += RefSeqMap[chr].substr(intron.beg , obs_intron_len);
			seq = GetReverseSeq(seq);
		}
	}
	return seq;
}

int Find_Mutant_5SS_Distance(int pos, intron_t intron)
{
	if (intron.strand) return (pos - intron.beg + 1);
	else return (intron.end - pos + 1);
}

int Find_Mutant_3SS_Distance(int pos, intron_t intron)
{
	if (intron.strand) return (intron.end - pos + 1);
	else return (pos - intron.beg + 1);
}

string GetMutantSeq5(int mt_x5_dist, string alt, string wt_seq)
{
	string mt_seq;
	int len, mt_pos;

	len = (int)wt_seq.length(); mt_seq = wt_seq;
	mt_pos = (len - 30) + (mt_x5_dist - 1);
	mt_seq[mt_pos] = alt[0];

	return mt_seq;
}

string GetMutantSeq3(int cut_pos, bool strand, string ref, string alt, string wt_seq)
{
	string mt_seq;

	if (strand) mt_seq = wt_seq.substr(0, cut_pos) + alt + wt_seq.substr(cut_pos + (int)ref.length());
	else mt_seq = wt_seq.substr(0, cut_pos) + GetReverseSeq(alt) + wt_seq.substr(cut_pos + (int)ref.length());
	//string str = string().assign(cut_pos, ' ') + string().assign(alt.length(), '*'); printf("   %s\n", str.c_str());
	return mt_seq;
}

/// VCF ---- begins here
bp_info_t Find_Mut_BP_Dist(string chr, int pos, int mt_dist, svm_bp_info_t svm_bp_info)
{
	bp_info_t bp_info;
	map<pair<string, int>, bp_info_t>::iterator iter;

	if ((iter = BPmap.find(make_pair(chr, pos))) == BPmap.end())
	{
		//printf("no bp_exp data!\n");
		bp_info.bp0 = svm_bp_info.bp0;
		bp_info.bp2 = svm_bp_info.bp2;
		bp_info.bp3 = svm_bp_info.bp3;

		if (svm_bp_info.mut_bp_diff == -4) bp_info.bp = 0;
		else bp_info.bp = svm_bp_info.mut_bp_diff;
	}
	else bp_info = iter->second;

	return bp_info;
}

void AddIntercept(Prediction_Info_t& Prediction)
{
	Prediction.model0_scr += Model0["Intercept"];
	Prediction.model1_scr += Model1["Intercept"];
	Prediction.model2_scr += Model2["Intercept"];
	Prediction.model3_scr += Model3["Intercept"];
	Prediction.model4_scr += Model4["Intercept"];
	Prediction.model5_scr += Model5["Intercept"];
}

vector<Prediction_Info_t> SplicingEffectAnalysis(string chr, int pos, string ref, string alt)
{
	intron_t intron;
	bp_info_t bp_info;
	Prediction_Info_t Prediction;
	map<pair<int, int>, bool> CacheMap;
	map<int, string>::iterator GeneIter;
	RNAcofold_info_t u2_mfe_c, u2_mfe_cmt;
	vector<Prediction_Info_t> PredictionVec;
	FreeEnergy_Info_t wt_dG_info, mt_dG_info;
	GC_info_t wt_gc_info, mt_gc_info, novel_gc_info;
	svm_bp_info_t wt_svm_bp, mt_svm_bp, novel_svm_bp;
	map<pair<string, int>, bp_info_t>::iterator exp_bp_iter;
	map<string, int> wt_MotifMap, mt_MotifMap, novel_MotifMap;
	Hollywood_info_t wt_hollywood_info, mt_hollywood_info, nss_hollywood_info;
	float f, wt_dG, mt_dG, wt_exdG, mt_exdG, wt23_dG, mt23_dG, wt_snp23_dG, mt_snp23_dG;
	string wt3seq, mt3seq, wt_subseq, mt_subseq, geneID, wt_seq, mt_seq, tmp, full_wt_seq, full_mt_seq, bp_wt_intron_seq, wt_2ndstr, mt_2ndstr;
	int i, exon_beg, exon_end, mt_x5_dist, mt_x3_dist, bp_right, wt_x3ss, mt_x3ss, new3ss, wt_Adenine, mt_Adenine, wt_openness, mt_openness, intron1_len, intron2_len;

	fprintf(stderr, "%s:%d=%s->%s\n", chr.c_str(), pos, ref.c_str(), alt.c_str());

	for (GeneIter = GeneRegionMap[chr].upper_bound(pos); GeneIter != GeneRegionMap[chr].end(); GeneIter--)
	{
		geneID = GeneIter->second; if (GeneMap[chr][geneID].beg > pos || GeneMap[chr][geneID].end < pos) continue;
		//fprintf(stderr, "check gene: %s [%d-%d]\n", geneID.c_str(), GeneMap[chr][geneID].beg, GeneMap[chr][geneID].end); 
		for (map<string, vector<intron_t> >::iterator TrantIter = GeneMap[chr][geneID].TranscriptIntronMap.begin(); TrantIter != GeneMap[chr][geneID].TranscriptIntronMap.end(); TrantIter++)
		{
			//fprintf(stderr, "\tcheck transcript: %s\n", TrantIter->first.c_str());
			for(i = 0; i<(int)TrantIter->second.size();i++)
			{
				intron = TrantIter->second[i];
				//fprintf(stderr, "\tcheck %s intron: [%d-%d]\n", TrantIter->first.c_str(), intron.beg, intron.end);
				//if (intron.beg <= (pos + 3) && intron.end >= (pos - 3)) fprintf(stderr, "****\n");
				if (CacheMap.find(make_pair(intron.beg, intron.end)) != CacheMap.end()) continue;
				CacheMap[make_pair(intron.beg, intron.end)] = true;
				
				//fprintf(stderr, "mt_x5_dist=%d, mt_x3_dist=%d\n", mt_x5_dist, mt_x3_dist);
				if (intron.beg <= (pos + 3) && intron.end >= (pos - 3) && (mt_x5_dist = Find_Mutant_5SS_Distance(pos, intron)) >= -2 && mt_x5_dist <= 30)
				{
					Prediction.pred_type = '5'; Prediction.bValid = true; Prediction.model0_scr = 0;
					Prediction.gene = geneID; Prediction.transcript = TrantIter->first; Prediction.intron = intron;
					Prediction.mt_dist = mt_x5_dist;
					fprintf(stderr, "check %s-%s[%d-%d] mt_dist5=%d\n", geneID.c_str(), TrantIter->first.c_str(), intron.beg, intron.end, mt_x5_dist);
					if ((wt_seq = GetObservedSeq5(chr, pos, i, TrantIter->second, exon_beg, exon_end, intron1_len, intron2_len)).length() == 0) continue;
	
					if (intron.strand) mt_seq = GetMutantSeq5(mt_x5_dist, alt, wt_seq);
					else mt_seq = GetMutantSeq5(mt_x5_dist, GetReverseSeq(alt), wt_seq);

					//printf("wt=%s\nmt=%s\n", wt_seq.c_str(), mt_seq.c_str());
					//printf("in1.len=%d, in2.len=%d\n", intron1_len, intron2_len);
					Prediction.model0_scr += Model0["in1.len"] * intron1_len;
					Prediction.model0_scr += Model0["in2.len"] * intron2_len;

					if ((intron.strand && ref[0] == 'C') || (!intron.strand && ref[0] == 'G')) Prediction.model0_scr += Model0["mt.snpC"];
					if ((intron.strand && alt[0] == 'G') || (!intron.strand && alt[0] == 'C')) Prediction.model0_scr += Model0["mt.snpG"];
					if ((intron.strand && ref[0] == 'T') || (!intron.strand && ref[0] == 'A')) Prediction.model0_scr += Model0["mt.snpT"];

					f = AvgPhyloP100_5(chr, intron);
					Prediction.model0_scr += Model0["X5ss_mean.phylop_100_per_base."] * f;
					f = AvgPhyloP100_5(chr, exon_beg, exon_end);
					Prediction.model0_scr += Model0["exon_mean.phylop_100_per_base."] * f;

					if (mt_x5_dist < 0) Prediction.model0_scr += mt_x5_dist * Model0["X5ss_position_left"];
					else Prediction.model0_scr += mt_x5_dist * Model0["X5ss_position_right"];

					wt_dG = Get_FreeEnergy_Info(wt_seq); mt_dG = Get_FreeEnergy_Info(mt_seq); f = mt_dG - wt_dG;
					Prediction.model0_scr += Model0["delta.dG"] * f;
					//Prediction.model0_scr += Model0["abs.delta.dG"] * abs(f);
					//printf("wt_dG=%f, mt_dG=%f\n", wt_dG, mt_dG);
					wt_subseq = wt_seq.substr(40, (int)wt_seq.length() - 70); wt_exdG = Get_FreeEnergy_Info(wt_subseq);
					mt_subseq = mt_seq.substr(40, (int)mt_seq.length() - 70); mt_exdG = Get_FreeEnergy_Info(mt_subseq);
					f = mt_exdG - wt_exdG;
					Prediction.model0_scr += Model0["mt.exdG"] * mt_exdG;
					Prediction.model0_scr += Model0["delta.exdG"] * f;
					Prediction.model0_scr += Model0["abs.delta.exdG"] * abs(f);
					//printf("wt_exdG=%f, mt_exdG=%f\n", wt_exdG, mt_exdG);

					wt_subseq = wt_seq.substr((int)wt_seq.length() - 33, 9); wt_hollywood_info = Run_Hollywood5(wt_subseq);
					mt_subseq = mt_seq.substr((int)mt_seq.length() - 33, 9); mt_hollywood_info = Run_Hollywood5(mt_subseq);
					//printf("wt_subseq=%s\nmaxent=%f\n", wt_subseq.c_str(), wt_hollywood_info.maxent);
					//printf("mt_subseq=%s\nmaxent=%f\n", mt_subseq.c_str(), mt_hollywood_info.maxent);
					//printf("maxent=%f, mdd=%f, mm=%f, wmm=%f\n", wt_hollywood_info.maxent, wt_hollywood_info.mdd, wt_hollywood_info.mm, wt_hollywood_info.wmm);
					if (wt_hollywood_info.bValid && mt_hollywood_info.bValid)
					{
						f = mt_hollywood_info.maxent - wt_hollywood_info.maxent;
						Prediction.model0_scr += Model0["delta_maxent"] * f;
					}
					wt_subseq = wt_seq.substr(20, 23); wt_hollywood_info = Run_Hollywood3(wt_subseq);
					Prediction.model0_scr += Model0["in1.3ss.wmm"] * wt_hollywood_info.wmm;

					Prediction.model0_scr += Model0["variant_snp.phylop_100."] * PhyloP100Map[chr][pos];
					Prediction.model0_scr += Model0["variant_snp.phastcons_100."] * PhastCons100Map[chr][pos];

					if (ClinVar_Query(chr, pos) || HGMD_Query(chr, pos)) Prediction.model0_scr += Model0["patho.sPathogenic/Likely pathogenic"];
					else Prediction.model0_scr += Model0["patho.sBenign/Likely benign"];

					AddIntercept(Prediction);
					f = exp(Prediction.model0_scr); Prediction.model0_exp = f / (1 + f);
					Prediction.effect = Prediction.model0_exp >= X5_Threshold ? true : false;
					PredictionVec.push_back(Prediction);
				}
				else if (intron.beg <= pos && intron.end >= pos && (mt_x3_dist = Find_Mutant_3SS_Distance(pos, intron)) > 3 && obs_intron_len > mt_x3_dist)
				{
					Prediction.pred_type = '3'; Prediction.bValid = true; 
					Prediction.gene = geneID; Prediction.transcript = TrantIter->first; Prediction.intron = intron;
					Prediction.model1_scr = Prediction.model2_scr = Prediction.model3_scr = Prediction.model4_scr = Prediction.model5_scr = 0;
					Prediction.mt_dist = mt_x3_dist;
					if ((wt_seq = GetObservedSeq3(chr, pos, intron)).length() == 0) continue;
					mt_seq = GetMutantSeq3((obs_intron_len - mt_x3_dist), intron.strand, ref, alt, wt_seq);
					//fprintf(stderr, "wt=%s\nmt=%s\n\n", wt_seq.c_str(), mt_seq.c_str());
					fprintf(stderr, "check %s-%s[%d-%d] mt_dist3=%d\n", geneID.c_str(), TrantIter->first.c_str(), intron.beg, intron.end, mt_x3_dist);

					full_wt_seq = LeftPrimer + wt_seq + RightPrimer; 
					full_mt_seq = LeftPrimer + mt_seq + RightPrimer;
						
					wt_x3ss = (int)full_wt_seq.length() - obs_exon_len - right_primer_len; // 1-base
					mt_x3ss = (int)full_mt_seq.length() - obs_exon_len - right_primer_len; // 1-base
					wt3seq = full_wt_seq.substr(wt_x3ss - mt_x3_dist - ref.length(), ref.length() + 2); //printf("wt3seq=%s (%s)\n", wt3seq.c_str(), (wt3seq == vec[7] ? "True" : "False"));
					mt3seq = full_mt_seq.substr(wt_x3ss - mt_x3_dist - alt.length(), alt.length() + 2); //printf("mt3seq=%s (%s)\n", mt3seq.c_str(), (mt3seq == vec[8] ? "True" : "False"));
					new3ss = FindNewAG(mt_x3_dist, wt3seq, mt3seq); Prediction.bNovelAG = (new3ss != -1 ? true : false);
					//wt_ag = Count_AG_Occurrences(wt3seq); mt_ag = Count_AG_Occurrences(mt3seq); //printf("wt.ag=%d, mt.ag=%d, new3ss=%d\n", wt_ag, mt_ag, new3ss);
						
					if (mt_seq[obs_intron_len - mt_x3_dist] == 'G') { Prediction.model2_scr += Model2["mt.snpG"]; if (debug) printf("mt.snpG\n"); }
					if (mt_seq[obs_intron_len - mt_x3_dist] == 'A') Prediction.model4_scr += Model4["mt.snpA"];
					if (mt_seq[obs_intron_len - mt_x3_dist] == 'G') Prediction.model4_scr += Model4["mt.snpG"];
					if (mt_seq[obs_intron_len - mt_x3_dist] == 'T') { Prediction.model2_scr += Model2["mt.snpT"]; if (debug) printf("mt.snpT\n"); }
					if (mt_seq[obs_intron_len - mt_x3_dist] == 'T') Prediction.model4_scr += Model4["mt.snpT"];

					bp_right = (0 - mt_x3_dist) < -25 ? 25 : mt_x3_dist;
					Prediction.model1_scr += Model1["bp.right"] * bp_right;
					Prediction.model2_scr += Model2["wt.mt_3ss_distance.right."] * mt_x3_dist; if (debug) printf("distance.right=%d\n", mt_x3_dist);
					Prediction.model3_scr += Model3["bp.right"] * bp_right;
					Prediction.model4_scr += Model4["bp.right"] * bp_right;
					Prediction.model4_scr += Model4["wt.mt_3ss_distance.right."] * mt_x3_dist;
					//A's count
					wt_subseq = wt_seq.substr(37, 23); wt_Adenine = Count_Adenine(wt_subseq);
					mt_subseq = mt_seq.substr(37, 23); mt_Adenine = Count_Adenine(mt_subseq);
					Prediction.model3_scr += Model3["A.wt18.40"] * wt_Adenine;
					Prediction.model4_scr += Model4["A.mt18.40"] * mt_Adenine;
					Prediction.model5_scr += Model5["A.wt18.40"] * wt_Adenine;
					//A's count around mt position
					wt_subseq = full_wt_seq.substr(left_primer_len + obs_intron_len - mt_x3_dist - 10, 21); wt_subseq[10] = 'N'; wt_Adenine = Count_Adenine(wt_subseq);
					mt_subseq = full_mt_seq.substr(left_primer_len + obs_intron_len - mt_x3_dist - 10, 21); mt_subseq[10] = 'N'; mt_Adenine = Count_Adenine(mt_subseq);
					if (mt_Adenine > 0) Prediction.model3_scr += Model3["Anumber.10yes"]; // to_be_checked
					//hollywood
					wt_hollywood_info = Run_Hollywood3(full_wt_seq.substr(wt_x3ss - 20, 23));
					mt_hollywood_info = Run_Hollywood3(full_mt_seq.substr(mt_x3ss - 20, 23));
					if (!wt_hollywood_info.bValid || !mt_hollywood_info.bValid)
					{
						fprintf(stderr, "Warning! hollywood failed!\n");
						Prediction.bValid = false;
					}
					//fprintf(stderr, "maxent=%f\nmm=%f\nwmm=%f\nmaxent.m=%f\nmm.m=%f\nwmm.m=%f\n", wt_hollywood_info.maxent, wt_hollywood_info.mm, wt_hollywood_info.wmm, mt_hollywood_info.maxent, mt_hollywood_info.mm, mt_hollywood_info.wmm);
					Prediction.model1_scr += Model1["mm.m"] * mt_hollywood_info.mm;
					Prediction.model1_scr += Model1["wmm.m"] * mt_hollywood_info.wmm;
					Prediction.model1_scr += Model1["delta_mm"] * (mt_hollywood_info.mm - wt_hollywood_info.mm);
					Prediction.model2_scr += Model2["mm.m"] * mt_hollywood_info.mm; if (debug) printf("mm.m=%f\n", mt_hollywood_info.mm);
					Prediction.model2_scr += Model2["maxent"] * wt_hollywood_info.maxent;; if (debug) printf("maxent=%f\n", wt_hollywood_info.maxent);
					Prediction.model2_scr += Model2["maxent.m"] * mt_hollywood_info.maxent; if (debug) printf("maxent.m=%f\n", mt_hollywood_info.maxent);
					Prediction.model2_scr += Model2["delta_mm"] * (mt_hollywood_info.mm - wt_hollywood_info.mm); if (debug) printf("delta_mm=%f\n", (mt_hollywood_info.mm - wt_hollywood_info.mm));
					Prediction.model3_scr += Model3["delta_mm"] * (mt_hollywood_info.mm - wt_hollywood_info.mm);
					Prediction.model3_scr += Model3["maxent.m"] * mt_hollywood_info.maxent;
					Prediction.model3_scr += Model3["mm.m"] * mt_hollywood_info.mm;
					Prediction.model4_scr += Model4["wmm.m"] * mt_hollywood_info.wmm;
					Prediction.model4_scr += Model4["delta_mm"] * (mt_hollywood_info.mm - wt_hollywood_info.mm);
					Prediction.model4_scr += Model4["maxent.m"] * mt_hollywood_info.maxent;
					Prediction.model4_scr += Model4["mm.m"] * mt_hollywood_info.mm;
					Prediction.model5_scr += Model5["wmm"] * wt_hollywood_info.wmm;
					Prediction.model5_scr += Model5["maxent"] * wt_hollywood_info.maxent;
					Prediction.model5_scr += Model5["mm"] * wt_hollywood_info.mm;
					//printf("delta_mm=%f\n", (mt_hollywood_info.mm - wt_hollywood_info.mm));
					if (Prediction.bNovelAG && (new3ss - 20) >= 0)
					{
						nss_hollywood_info = Run_Hollywood3(full_mt_seq.substr(new3ss - 20, 23)); 
						Prediction.model1_scr += Model1["novel.wmm"] * nss_hollywood_info.wmm;
						Prediction.model3_scr += Model3["novel_mt.wmm"] * (nss_hollywood_info.wmm - mt_hollywood_info.wmm);
						Prediction.model3_scr += Model3["novel_mt.maxent"] * (nss_hollywood_info.maxent - mt_hollywood_info.maxent);
					}
					//svm_bp_finder
					wt_svm_bp = Run_SVM_BP_Finder_Docker(mt_x3_dist, full_wt_seq.substr(wt_x3ss - 100, 100)); //printf("wt_bp_seq=%s\n", full_wt_seq.substr(wt_x3ss - 100, 100).c_str());
					if (!wt_svm_bp.bValid)
					{
						fprintf(stderr, "warning! svm_bp failed!\n");
						Prediction.bValid = false;
					}
					if (!Prediction.bNovelAG) mt_svm_bp = Run_SVM_BP_Finder_Docker(mt_x3_dist, full_mt_seq.substr(mt_x3ss - 100, 100)); //printf("mt_bp_seq=%s\n", full_mt_seq.substr(mt_x3ss - 100, 100).c_str());

					//printf("ppt_scr.m=%d, svm_scr.m=%f\n", mt_svm_bp.ppt_scr, mt_svm_bp.svm_scr);
					bp_info = Find_Mut_BP_Dist(chr, pos, mt_x3_dist, wt_svm_bp); //printf("bp = %d\n", bp_info.bp);
					if ((exp_bp_iter = BPmap.find(make_pair(chr, pos))) != BPmap.end())
					{
						if (exp_bp_iter->second.bp0) Prediction.model1_scr += Model1["bp_position0"];
						if (exp_bp_iter->second.bp0) { Prediction.model2_scr += Model2["bp_position0"]; if (debug) printf("bp0\n"); }
						if (exp_bp_iter->second.bp3) { Prediction.model2_scr += Model2["bp_position-3"]; if (debug) printf("bp3\n"); }
						if ((exp_bp_iter->second.bp0 || exp_bp_iter->second.bp2 || exp_bp_iter->second.bp3) && full_wt_seq[wt_x3ss - mt_x3_dist - exp_bp_iter->second.bp] == 'G')
						{
							Prediction.model2_scr += Model2["bp.snpG"]; if (debug) printf("bp.snpG\n");
							Prediction.model4_scr += Model4["bp.snpG"]; //printf("bp=%d, bp.snp=%c\n", bp_info.bp, full_wt_seq[wt_x3ss - mt_dist - bp_info.bp]);
						}
					}
					if (!Prediction.bNovelAG)
					{
						Prediction.model2_scr += Model2["bp_scr"] * wt_svm_bp.bp_scr; if (debug) printf("bp_scr=%f\n", wt_svm_bp.bp_scr);
						Prediction.model2_scr += Model2["ppt_scr"] * wt_svm_bp.ppt_scr; if (debug) printf("ppt_scr=%d\n", wt_svm_bp.ppt_scr);
						Prediction.model2_scr += Model2["ppt_scr.m"] * mt_svm_bp.ppt_scr; if (debug) printf("ppt_scr.m=%d\n", mt_svm_bp.ppt_scr);
						Prediction.model2_scr += Model2["dsvm"] * (mt_svm_bp.svm_scr - wt_svm_bp.svm_scr); if (debug) printf("dsvm=%f\n", (mt_svm_bp.svm_scr - wt_svm_bp.svm_scr));
						Prediction.model2_scr += Model2["ss_dist.m"] * (mt_svm_bp.ss_dist); if (debug) printf("ss_dist.m=%d\n", (mt_svm_bp.ss_dist));
						Prediction.model4_scr += Model4["bp_scr.m"] * mt_svm_bp.bp_scr;
						Prediction.model4_scr += Model4["dsvm"] * (mt_svm_bp.svm_scr - wt_svm_bp.svm_scr);
						Prediction.model4_scr += Model4["ss_dist.m"] * mt_svm_bp.ss_dist;
						Prediction.model4_scr += Model4["ppt_scr"] * wt_svm_bp.ppt_scr;
						Prediction.model4_scr += Model4["svm_bp_ct"] * wt_svm_bp.svm_bp_ct;
					}
					Prediction.model3_scr += Model3["ppt_scr"] * wt_svm_bp.ppt_scr;
					Prediction.model5_scr += Model5["ss_dist"] * (wt_svm_bp.ss_dist);
					Prediction.model5_scr += Model5["svm_scr"] * wt_svm_bp.svm_scr;
					Prediction.model5_scr += Model5["svm_bp_ct"] * wt_svm_bp.svm_bp_ct;

					if (wt_svm_bp.ss_dist == mt_svm_bp.ss_dist) {
						Prediction.model2_scr += Model2["same_bpTRUE"]; if (debug) printf("same_bpTRUE\n");
					}
					if (Prediction.bNovelAG && (new3ss - 100) >= 0)
					{
						novel_svm_bp = Run_SVM_BP_Finder_Docker(mt_x3_dist, full_mt_seq.substr(new3ss - 100, 100));
						Prediction.model3_scr += Model3["novel.svm_scr"] * novel_svm_bp.svm_scr;
					}
					//RNAcofold
					if ((wt_x3ss - mt_x3_dist - bp_info.bp - 5 >= 0))
					{
							wt_subseq = full_wt_seq.substr(wt_x3ss - mt_x3_dist - bp_info.bp - 5, 9); u2_mfe_c = Run_RNAcofold(wt_subseq);
							Prediction.model1_scr += Model1["U2.mfe.c"] * u2_mfe_c.mfe;
							Prediction.model5_scr += Model5["U2.mfe.c"] * u2_mfe_c.mfe;
					}
					if (!Prediction.bNovelAG)
					{
						mt_subseq = full_mt_seq.substr(wt_x3ss - mt_x3_dist - bp_info.bp - 5, 9); u2_mfe_cmt = Run_RNAcofold(mt_subseq);
						Prediction.model2_scr += Model2["U2.mfe.cmt"] * u2_mfe_cmt.mfe; if (debug) printf("U2.mfe.cmt=%f\n", u2_mfe_cmt.mfe);
					}
					// GC_ratio
					wt_gc_info = Get_GC_Info(wt_x3ss, full_wt_seq); mt_gc_info = Get_GC_Info(mt_x3ss, full_mt_seq);
					//printf("mt.wt_delta_gc.percent=%f\n", 100 * (mt_gc_info.gc - wt_gc_info.gc));
					Prediction.model2_scr += Model2["mt.wt_delta_gc.percent."] * (100 * (mt_gc_info.gc - wt_gc_info.gc)); if (debug) printf("mt.wt_delta_gc.percent=%f\n", (100 * (mt_gc_info.gc - wt_gc_info.gc)));
					Prediction.model3_scr += Model3["mt.wt_delta_gc.percent."] * (100 * (mt_gc_info.gc - wt_gc_info.gc));
					Prediction.model4_scr += Model4["mt.ex.in._delta_gc.percent."] * (100 * (mt_gc_info.ex_gc - mt_gc_info.in_gc));
					if (Prediction.bNovelAG)
					{
						novel_gc_info = Get_GC_Info(new3ss, mt_x3ss, full_mt_seq);
						Prediction.model1_scr += Model1["novel.ex.short_in.gc"] * (novel_gc_info.short_ex_gc - novel_gc_info.in_gc);
						Prediction.model3_scr += Model3["novel.ex.short_in.gc"] * (novel_gc_info.short_ex_gc - novel_gc_info.in_gc);
						Prediction.model3_scr += Model3["novel.ex.short.gc"] * novel_gc_info.short_ex_gc;
					}
					Prediction.model4_scr += Model4["wt.ex.in._delta_gc.percent."] * (100 * (wt_gc_info.ex_gc - wt_gc_info.in_gc));
					Prediction.model5_scr += Model5["wt.ex.in._delta_gc.percent."] * (100 * (wt_gc_info.ex_gc - wt_gc_info.in_gc));
					// PhyloP
					if (Prediction.bNovelAG) Prediction.model1_scr += Model1["PhyloP100"] * PhyloP100Map[chr][pos];
					if (!Prediction.bNovelAG)
					{
						Prediction.model2_scr += Model2["PhyloP20"] * PhyloP20Map[chr][pos]; if (debug) printf("PhyloP20=%f\n", PhyloP20Map[chr][pos]);
						Prediction.model2_scr += Model2["PhyloP7"] * PhyloP7Map[chr][pos]; if (debug) printf("PhyloP7=%f\n", PhyloP7Map[chr][pos]);
						Prediction.model2_scr += Model2["PhyloP100"] * PhyloP100Map[chr][pos]; if (debug) printf("PhyloP100=%f\n", PhyloP100Map[chr][pos]);
						Prediction.model4_scr += Model4["PhastCons100"] * PhastCons100Map[chr][pos];
						Prediction.model4_scr += Model4["PhyloP100"] * PhyloP100Map[chr][pos];
						Prediction.model4_scr += Model4["PhyloP20"] * PhyloP20Map[chr][pos];
					}
					// Free energy (full length)
					wt_dG_info = Get_FreeEnergy_Info(wt_x3ss, full_wt_seq);
					mt_dG_info = Get_FreeEnergy_Info(mt_x3ss, full_mt_seq);
					//printf("delta.indG=%.1f, mt.dG=%.1f, mt.exdG=%.1f\n", (mt_dG_info.in_dG - wt_dG_info.in_dG), mt_dG_info.dG, mt_dG_info.ex_dG);
					Prediction.model1_scr += Model1["delta.indG"] * (mt_dG_info.in_dG - wt_dG_info.in_dG);
					Prediction.model2_scr += Model2["delta.indG"] * (mt_dG_info.in_dG - wt_dG_info.in_dG); if (debug) printf("delta.indG=%f\n", (mt_dG_info.in_dG - wt_dG_info.in_dG));
					Prediction.model2_scr += Model2["mt.dG"] * mt_dG_info.dG; if (debug) printf("mt.dG=%f\n", mt_dG_info.dG);
					Prediction.model2_scr += Model2["mt.exdG"] * mt_dG_info.ex_dG; if (debug) printf("mt.exdG=%f\n", mt_dG_info.ex_dG);
					Prediction.model3_scr += Model3["wt.exdG"] * wt_dG_info.ex_dG;
					Prediction.model3_scr += Model3["delta.indG"] * (mt_dG_info.in_dG - wt_dG_info.in_dG);
					Prediction.model4_scr += Model4["mt.dG"] * mt_dG_info.dG;
					Prediction.model4_scr += Model4["mt.exdG"] * mt_dG_info.ex_dG;
					Prediction.model4_scr += Model4["mt.indG"] * mt_dG_info.in_dG;
					Prediction.model4_scr += Model4["wt.exdG"] * wt_dG_info.ex_dG;
					Prediction.model4_scr += Model4["wt.indG"] * wt_dG_info.in_dG;
					Prediction.model5_scr += Model5["wt.dG"] * wt_dG_info.dG;
					Prediction.model5_scr += Model5["wt.exdG"] * wt_dG_info.ex_dG;
					Prediction.model5_scr += Model5["wt.indG"] * wt_dG_info.in_dG;
					// openness
					wt_2ndstr = RNAfold_str(full_wt_seq); mt_2ndstr = RNAfold_str(full_mt_seq);
					tmp = wt_2ndstr.substr(wt_x3ss - mt_x3_dist - 5, 11); wt_openness = CountOpenness(tmp);
					tmp = mt_2ndstr.substr(mt_x3ss - mt_x3_dist - 5, 11); mt_openness = CountOpenness(tmp); //printf("mt.openness=%d\n", mt_openness);
					Prediction.model2_scr += Model2["mt.openness"] * mt_openness; if (debug) printf("mt.openness=%d\n", mt_openness);
					Prediction.model3_scr += Model3["mt.openness"] * mt_openness;
					Prediction.model4_scr += Model4["mt_wt.openness"] * (mt_openness - wt_openness);
					tmp = wt_2ndstr.substr(wt_x3ss - 20, 17); Prediction.model5_scr += Model5["wt17open"] * CountOpenness(tmp);
					tmp = mt_2ndstr.substr(mt_x3ss - 20, 17); Prediction.model3_scr += Model3["mt17open"] * CountOpenness(tmp);
					// RNA fold (23nt)
					wt_subseq = full_wt_seq.substr(wt_x3ss - 20, 23); wt23_dG = RNAfold(wt_subseq);
					mt_subseq = full_mt_seq.substr(mt_x3ss - 20, 23); mt23_dG = RNAfold(mt_subseq);
					Prediction.model3_scr += Model3["abs.delta23.dG"] * abs(mt23_dG - wt23_dG);
					Prediction.model4_scr += Model4["wt23.dG"] * wt23_dG;

					wt_subseq = full_wt_seq.substr(wt_x3ss - mt_x3_dist - 19, 23); wt_snp23_dG = RNAfold(wt_subseq);
					mt_subseq = full_mt_seq.substr(mt_x3ss - mt_x3_dist - 19, 23); mt_snp23_dG = RNAfold(mt_subseq);
					Prediction.model4_scr += Model4["abs.delta.snp23.dG"] * abs(mt_snp23_dG - wt_snp23_dG);
					Prediction.model4_scr += Model4["delta.snp23.dG"] * (mt_snp23_dG - wt_snp23_dG);

					if (HGMD_Query(chr, pos)) Prediction.model1_scr += Model1["databasehgmd"];
					if (ClinVar_Query(chr, pos) || HGMD_Query(chr, pos))
					{
						Prediction.model2_scr += Model2["simp.variant_classPathogenic/Likely pathogenic"]; if (debug) printf("pathogenic\n");
						Prediction.model3_scr += Model3["simp.variant_classPathogenic/Likely pathogenic"];
						Prediction.model4_scr += Model4["simp.variant_classPathogenic/Likely pathogenic"];
					}
					else
					{
						Prediction.model2_scr += Model2["simp.variant_classBenign/Likely benign"]; if (debug) printf("benign\n");
						if (dbSNPmap_Query(chr, pos, intron.strand)) {
							Prediction.model2_scr += Model2["databasedbSNP"]; if (debug) printf("dbSNP\n");
						}
					}
					if (AltSpliceRegion_Query(chr, pos, intron.strand))
					{
						Prediction.model2_scr += Model2["UCSCAltTRUE"]; if (debug) printf("UCSCAlt\n");
						Prediction.model3_scr += Model3["UCSCAltTRUE"];
						Prediction.model4_scr += Model4["UCSCAltTRUE"];
						Prediction.model5_scr += Model5["UCSCAltTRUE"];
					}
					//ATtRACT
					wt_subseq = wt_seq.substr(0, obs_intron_len); wt_MotifMap = MotifScan(wt_subseq); // wt.intron
					Prediction.model1_scr += Model1["DHX58.intron.wt"] * wt_MotifMap["DHX58"];
					Prediction.model3_scr += Model3["HNRNPA1L2.intron.wt"] * wt_MotifMap["HNRNPA1L2"];
					Prediction.model4_scr += Model4["ENOX1.intron.wt"] * wt_MotifMap["ENOX1"];
					Prediction.model4_scr += Model4["F2.intron.wt"] * wt_MotifMap["F2"];
					Prediction.model4_scr += Model4["NOVA1.intron.wt"] * wt_MotifMap["NOVA1"];
					Prediction.model4_scr += Model4["PHAX.intron.wt"] * wt_MotifMap["PHAX"];
					Prediction.model4_scr += Model4["RBM42.intron.wt"] * wt_MotifMap["RBM42"];
					Prediction.model4_scr += Model4["RBM46.intron.wt"] * wt_MotifMap["RBM46"];
					Prediction.model4_scr += Model4["YBX2.intron.wt"] * wt_MotifMap["YBX2"];
					Prediction.model4_scr += Model4["ZC3H10.intron.wt"] * wt_MotifMap["ZC3H10"];
					Prediction.model5_scr += Model5["AKAP1.intron.wt"] * wt_MotifMap["AKAP1"];
					Prediction.model5_scr += Model5["CELF4.intron.wt"] * wt_MotifMap["CELF4"];
					Prediction.model5_scr += Model5["CELF5.intron.wt"] * wt_MotifMap["CELF5"];
					Prediction.model5_scr += Model5["DHX58.intron.wt"] * wt_MotifMap["DHX58"];
					Prediction.model5_scr += Model5["EIF4B.intron.wt"] * wt_MotifMap["EIF4B"];
					Prediction.model5_scr += Model5["ESRP1.intron.wt"] * wt_MotifMap["ESRP1"];
					Prediction.model5_scr += Model5["F2.intron.wt"] * wt_MotifMap["F2"];
					Prediction.model5_scr += Model5["GRSF1.intron.wt"] * wt_MotifMap["GRSF1"];
					Prediction.model5_scr += Model5["HNRNPCL1.intron.wt"] * wt_MotifMap["HNRNPCL1"];
					Prediction.model5_scr += Model5["PABPN1.intron.wt"] * wt_MotifMap["PABPN1"];
					Prediction.model5_scr += Model5["PTBP2.intron.wt"] * wt_MotifMap["PTBP2"];
					Prediction.model5_scr += Model5["RBM28.intron.wt"] * wt_MotifMap["RBM28"];
					Prediction.model5_scr += Model5["RBM46.intron.wt"] * wt_MotifMap["RBM46"];
					Prediction.model5_scr += Model5["RNASEL.intron.wt"] * wt_MotifMap["RNASEL"];
					Prediction.model5_scr += Model5["SAMD4A.intron.wt"] * wt_MotifMap["SAMD4A"];
					Prediction.model5_scr += Model5["ZC3H10.intron.wt"] * wt_MotifMap["ZC3H10"];

					mt_subseq = mt_seq.substr(0, (int)mt_seq.length() - obs_exon_len); mt_MotifMap = MotifScan(mt_subseq); // mt.intron
					Prediction.model1_scr += Model1["DHX58.intron.mt"] * mt_MotifMap["DHX58"];
					Prediction.model3_scr += Model3["CELF4.intron.mt"] * mt_MotifMap["CELF4"];
					Prediction.model3_scr += Model3["CELF5.intron.mt"] * mt_MotifMap["CELF5"];
					Prediction.model3_scr += Model3["CMTR1.intron.mt"] * mt_MotifMap["CMTR1"];
					Prediction.model3_scr += Model3["PABPC5.intron.mt"] * mt_MotifMap["PABPC5"];
					Prediction.model4_scr += Model4["AKAP1.intron.mt"] * mt_MotifMap["AKAP1"];
					Prediction.model4_scr += Model4["GRSF1.intron.mt"] * mt_MotifMap["GRSF1"];
					Prediction.model4_scr += Model4["PABPC3.intron.mt"] * mt_MotifMap["PABPC3"];
					Prediction.model4_scr += Model4["PABPN1.intron.mt"] * mt_MotifMap["PABPN1"];
					Prediction.model4_scr += Model4["RBM24.intron.mt"] * mt_MotifMap["RBM24"];
					Prediction.model4_scr += Model4["RBMS1.intron.mt"] * mt_MotifMap["RBMS1"];
					Prediction.model4_scr += Model4["SAMD4A.intron.mt"] * mt_MotifMap["SAMD4A"];

					Prediction.model4_scr += Model4["delta.CELF4.intron"] * (mt_MotifMap["CELF4"] - wt_MotifMap["CELF4"]);
					Prediction.model4_scr += Model4["delta.CELF5.intron"] * (mt_MotifMap["CELF5"] - wt_MotifMap["CELF5"]);
					Prediction.model4_scr += Model4["delta.CELF6.intron"] * (mt_MotifMap["CELF6"] - wt_MotifMap["CELF6"]);

					wt_subseq = wt_seq.substr(obs_intron_len); wt_MotifMap = MotifScan(wt_subseq); // wt.exon
					Prediction.model1_scr += Model1["CMTR1.exon.fl"] * wt_MotifMap["CMTR1"];
					Prediction.model1_scr += Model1["SAMD4A.exon.fl"] * wt_MotifMap["SAMD4A"];
					Prediction.model3_scr += Model3["CMTR1.exon.fl"] * wt_MotifMap["CMTR1"];
					Prediction.model3_scr += Model3["RNASEL.exon.fl"] * wt_MotifMap["RNASEL"];
					Prediction.model3_scr += Model3["SAMD4A.exon.fl"] * wt_MotifMap["SAMD4A"];
					Prediction.model3_scr += Model3["SART3.exon.fl"] * wt_MotifMap["SART3"];
					Prediction.model4_scr += Model4["ACO1.exon.fl"] * wt_MotifMap["ACO1"];
					Prediction.model4_scr += Model4["CMTR1.exon.fl"] * wt_MotifMap["CMTR1"];
					Prediction.model4_scr += Model4["CNOT4.exon.fl"] * wt_MotifMap["CNOT4"];
					Prediction.model4_scr += Model4["CPEB1.exon.fl"] * wt_MotifMap["CPEB1"];
					Prediction.model4_scr += Model4["CPEB2.exon.fl"] * wt_MotifMap["CPEB2"];
					Prediction.model4_scr += Model4["DHX58.exon.fl"] * wt_MotifMap["DHX58"];
					Prediction.model4_scr += Model4["GRSF1.exon.fl"] * wt_MotifMap["GRSF1"];
					Prediction.model4_scr += Model4["HNRNPA1L2.exon.fl"] * wt_MotifMap["HNRNPA1L2"];
					Prediction.model4_scr += Model4["PHAX.exon.fl"] * wt_MotifMap["PHAX"];
					Prediction.model4_scr += Model4["PTBP2.exon.fl"] * wt_MotifMap["PTBP2"];
					Prediction.model4_scr += Model4["RBM28.exon.fl"] * wt_MotifMap["RBM28"];
					Prediction.model4_scr += Model4["RBM46.exon.fl"] * wt_MotifMap["RBM46"];
					Prediction.model4_scr += Model4["RNASEL.exon.fl"] * wt_MotifMap["RNASEL"];
					Prediction.model4_scr += Model4["SART3.exon.fl"] * wt_MotifMap["SART3"];
					Prediction.model4_scr += Model4["ZNF638.exon.fl"] * wt_MotifMap["ZNF638"];
					Prediction.model5_scr += Model5["ACO1.exon.fl"] * wt_MotifMap["ACO1"];
					Prediction.model5_scr += Model5["CELF4.exon.fl"] * wt_MotifMap["CELF4"];
					Prediction.model5_scr += Model5["CMTR1.exon.fl"] * wt_MotifMap["CMTR1"];
					Prediction.model5_scr += Model5["CNOT4.exon.fl"] * wt_MotifMap["CNOT4"];
					Prediction.model5_scr += Model5["CPEB1.exon.fl"] * wt_MotifMap["CPEB1"];
					Prediction.model5_scr += Model5["ESRP2.exon.fl"] * wt_MotifMap["ESRP2"];
					Prediction.model5_scr += Model5["FXR2.exon.fl"] * wt_MotifMap["FXR2"];
					Prediction.model5_scr += Model5["GRSF1.exon.fl"] * wt_MotifMap["GRSF1"];
					Prediction.model5_scr += Model5["HNRNPA1L2.exon.fl"] * wt_MotifMap["HNRNPA1L2"];
					Prediction.model5_scr += Model5["PABPC4.exon.fl"] * wt_MotifMap["PABPC4"];
					Prediction.model5_scr += Model5["PABPC5.exon.fl"] * wt_MotifMap["PABPC5"];
					Prediction.model5_scr += Model5["PHAX.exon.fl"] * wt_MotifMap["PHAX"];
					Prediction.model5_scr += Model5["PIWIL1.exon.fl"] * wt_MotifMap["PIWIL1"];
					Prediction.model5_scr += Model5["PTBP2.exon.fl"] * wt_MotifMap["PTBP2"];
					Prediction.model5_scr += Model5["RBM46.exon.fl"] * wt_MotifMap["RBM46"];
					Prediction.model5_scr += Model5["RNASEL.exon.fl"] * wt_MotifMap["RNASEL"];
					Prediction.model5_scr += Model5["SART3.exon.fl"] * wt_MotifMap["SART3"];

					if (Prediction.bNovelAG)
					{
						mt_subseq = mt_seq.substr(0, obs_intron_len - (mt_x3ss - new3ss)); novel_MotifMap = MotifScan(mt_subseq); // novel.intron
						Prediction.model1_scr += Model1["SF1.intron"] * novel_MotifMap["SF1"];
						Prediction.model1_scr += Model1["AGO1.intron"] * novel_MotifMap["AGO1"];
						Prediction.model1_scr += Model1["OAS1.intron"] * novel_MotifMap["OAS1"];
						Prediction.model1_scr += Model1["SRP68.intron"] * novel_MotifMap["SRP68"];

						mt_subseq = mt_seq.substr(obs_intron_len - (mt_x3ss - new3ss)); novel_MotifMap = MotifScan(mt_subseq); // novel.exon
						Prediction.model1_scr += Model1["PUM1.exon"] * novel_MotifMap["PUM1"];
						Prediction.model1_scr += Model1["HNRNPLL.exon"] * novel_MotifMap["HNRNPLL"];
						Prediction.model1_scr += Model1["SRSF5.exon"] * novel_MotifMap["SRSF5"];
						Prediction.model1_scr += Model1["SRSF7.exon"] * novel_MotifMap["SRSF7"];
					}
					AddIntercept(Prediction);
					if (debug) printf("yhat1: %f\nyhat2: %f\nyhat3: %f\nyhat4: %f\nyhat5: %f\n\n", Prediction.model1_scr, Prediction.model2_scr, Prediction.model3_scr, Prediction.model4_scr, Prediction.model5_scr);
					f = exp(Prediction.model1_scr); Prediction.model1_exp = f / (1 + f);
					f = exp(Prediction.model2_scr); Prediction.model2_exp = f / (1 + f);
					f = exp(Prediction.model3_scr); Prediction.model3_exp = f / (1 + f);
					f = exp(Prediction.model4_scr); Prediction.model4_exp = f / (1 + f);
					f = exp(Prediction.model5_scr); Prediction.model5_exp = f / (1 + f);
					if (debug) printf("model1: %f\nmodel2: %f\nmodel3: %f\nmodel4: %f\nmodel5: %f\n\n", Prediction.model1_exp, Prediction.model2_exp, Prediction.model3_exp, Prediction.model4_exp, Prediction.model5_exp);

					if (Prediction.bNovelAG) Prediction.effect = Prediction.model1_exp >= AG_Threshold ? true : false;
					else Prediction.effect = Prediction.model2_exp >= NonAG_Threshold ? true : false;
					PredictionVec.push_back(Prediction);
				}
				else if (intron.beg <= pos && intron.end >= pos)//out of scope
				{
					Prediction.bValid = false;
					Prediction.gene = geneID; Prediction.transcript = TrantIter->first; Prediction.intron = intron;
					PredictionVec.push_back(Prediction);
				}
			}
		}
	}
	return PredictionVec;
}

void FindIntronSeqFromTestDataSet()
{
	//char strand;
	//fstream file;
	//int i, n=0, pos;
	//stringstream ss;
	//vector<string> vec(5);
	//string str, chr, ref, alt, tmp;
	//vector<Prediction_Info_t> PredictionVec;

	//file.open("/home/bbsc/lab_data/N420/20210037/kernel/test_data.csv", ios_base::in); getline(file, str);
	//while (!file.eof())
	//{
	//	getline(file, str); if (str == "") break;
	//	n++; ss.clear(); ss.str(str);
	//	//if (n+1 != 14) continue;
	//	for (i = 0; i < 5; i++) getline(ss, vec[i], ',');
	//	i = vec[0].find_first_of(':'); chr = vec[0].substr(0, i);
	//	i++; tmp = vec[0].substr(i, vec[0].find_first_of('-', i) - i); pos = atoi(tmp.c_str());
	//	ref = vec[1]; alt = vec[2];
	//	if ((strand = vec[3][0]) == '-') {
	//		ref = GetReverseSeq(ref); 
	//		alt = GetReverseSeq(alt);
	//	}
	//	//if (vec[4] != "BRCA1_c.4987-20A>G") continue;
	//	if (ref.length() == 1 && alt.length() == 1)
	//	{
	//		fprintf(stderr, "%d (%s): pos = %d, ref = %s, alt = %s\n", n + 1, vec[4].c_str(), pos, ref.c_str(), alt.c_str());
	//		PredictionVec = SplicingEffectAnalysis(chr, pos, ref, alt);
	//		
	//		for (i = 0; i < (int)PredictionVec.size(); i++)
	//		{
	//			if (!PredictionVec[i].bValid)
	//			{
	//				fprintf(stderr, "Warning! prediction failed!\n");
	//				exit(1);
	//			}
	//			printf("%s\t%f\t%s\n", vec[4].c_str(), PredictionVec[i].bNovelAG ? PredictionVec[i].model1_scr : PredictionVec[i].model2_scr, PredictionVec[i].effect ? "sig" : "insig");
	//		}
	//	}
	//	//if (n == 11) break;
	//}
	//file.close();
	//chr = ; pos = ; ref = "T"; alt = "C"; //3'
	//chr = "chr1"; pos = 113618627; ref = "T"; alt = "A"; //5'forward strand
	//chr = "chr1"; pos = 119042244; ref = "G"; alt = "C"; //5'reverse strand
	vector<Prediction_Info_t> PredictionVec;
	PredictionVec = SplicingEffectAnalysis("chr1", 70418035, "A", "T");
	PredictionVec = SplicingEffectAnalysis("chr19", 43990112, "G", "C");

}
