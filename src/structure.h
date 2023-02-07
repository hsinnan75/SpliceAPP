#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <map>
#include <ctime>
#include <algorithm>
#include <arpa/inet.h>
#include <assert.h>
#include <errno.h>
#include <netinet/in.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <netdb.h>
#include <unistd.h>
#include <ctime>

#define obs_exon_len 37
#define obs_intron_len 78
#define left_primer_len 84
#define right_primer_len 60

#define LeftPrimer "aagttcagcgtgtccggcgagggcgaggTGAGTATGGATCCtgaggaatgcagggagtggggtgggctggggcctgaagacctg"
#define RightPrimer "CCCCGGGGTCATGTGCGCCTTGCCCGCTGCCTTGAGGAACTACAGAGACAGGAGCCTTCG"

#define MAXLINE 40960

using namespace std;

typedef struct
{
	int beg;
	int end;
	bool strand;
	int trant_idx;
} exon_t;

typedef exon_t intron_t;

typedef struct
{
	int beg;
	int end;
	//vector<string> TrantVec;
	//vector<intron_t> IntronVec;
	map<string, vector<exon_t> > TranscriptExonMap;
	map<string, vector<intron_t> > TranscriptIntronMap;
} Gene_t;

typedef struct
{
	int bp;
	bool bp0;
	bool bp2;
	bool bp3;
} bp_info_t;

typedef struct
{
	bool bValid;
	string bp_seq;
	int svm_bp_ct;
	int ss_dist;
	float bp_scr;
	int ppt_scr;
	float svm_scr;
	int mut_bp_diff; //0,-2,-3,-4(NA)
	bool bp0, bp2, bp3;
} svm_bp_info_t;

typedef struct
{
	bool bValid;
	float maxent; //maxent
	float mm; //mm
	float wmm; // wmm
	float mdd;
} Hollywood_info_t;

typedef struct
{
	float mfe;
	string str;
} RNAcofold_info_t;

typedef struct
{
	float phyloP100;
	float phyloP30;
	float phyloP20;
	float phyloP7;
} PhyloP_info_t;

typedef struct
{
	float gc;
	float ex_gc;
	float in_gc;
	float short_ex_gc;
} GC_info_t;

typedef struct
{
	float dG;
	float ex_dG;
	float in_dG;
} FreeEnergy_Info_t;

typedef struct
{
	int beg;
	int end;
	bool strand;
} AltSpliceInfo_t;

typedef struct
{
	char pred_type; // 5=5' or 3=3'
	int mt_dist;
	bool bValid;
	string gene;
	string transcript;
	bool effect;
	bool bNovelAG;
	intron_t intron;
	float model0_scr;
	float model1_scr;
	float model2_scr;
	float model3_scr;
	float model4_scr;
	float model5_scr;
	float model0_exp;
	float model1_exp;
	float model2_exp;
	float model3_exp;
	float model4_exp;
	float model5_exp;
} Prediction_Info_t;

//extern map<string, int> TrantIDmap;
extern map<string, string> RefSeqMap;
extern map<string, string> AttractMap;
extern map<string, map<int, int> > ClinVarMap;
extern map<pair<string, int>, bp_info_t> BPmap;
extern map<string, map<string, Gene_t> > GeneMap;
extern map<string, map<int, string> > GeneRegionMap;
extern map<pair<string, int>, bool> HGMDmap, dbSNPmap;
extern map<string, map<int, AltSpliceInfo_t> > AltSpliceMap;
extern map<string, float> Model0, Model1, Model2, Model3, Model4, Model5;
extern map<string, vector<float> > PhyloP100Map, PhastCons100Map, PhyloP20Map, PhyloP7Map;

extern void GetIntronRegions();
extern float GC_count(string& seq);
extern string GetReverseSeq(string seq);
extern int Count_AG_Occurrences(string str);
extern void InitWebKitFormBoundary(char* boundary);
extern int FindNewAG(int mt_dist, string wt3seq, string mt3seq);

extern void LoadClinVarMap();
extern string ToUpperStr(string seq);
extern int CountOpenness(string str);
extern int Count_Adenine(string seq);
extern void LocalMode(char* filename);
extern float RNAfold(string query_seq);
extern string RNAfold_str(string query_seq);
extern bool HGMD_Query(string chr, int pos);
extern float Get_FreeEnergy_Info(string seq);
extern map<string, int> MotifScan(string seq);
extern bool ClinVar_Query(string chr, int pos);
extern Hollywood_info_t Run_Hollywood3(string seq);
extern Hollywood_info_t Run_Hollywood5(string seq);
extern RNAcofold_info_t Run_RNAcofold(string bp_seq);
extern GC_info_t Get_GC_Info(int wt_x3ss, string seq);
extern bool shellCmd(const string& cmd, string& result);
extern float AvgPhyloP100_5(string chr, intron_t intron);
extern float AvgPhyloP100_5(string chr, int beg, int end);
extern bool dbSNPmap_Query(string chr, int pos, bool strand);
extern float Get_Phylop(string chr, int pos, string database);
extern GC_info_t Get_GC_Info(int new3ss, int mt_x3ss, string seq);
extern FreeEnergy_Info_t Get_FreeEnergy_Info(int x3ss, string seq);
extern bool AltSpliceRegion_Query(string chr, int pos, bool strand);
extern svm_bp_info_t Run_SVM_BP_Finder(int mt_dist, string intron_seq);
extern svm_bp_info_t Run_SVM_BP_Finder_Docker(int mt_dist, string intron_seq);
extern vector<Prediction_Info_t> SplicingEffectAnalysis(string chr, int pos, string ref, string alt);

