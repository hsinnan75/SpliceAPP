#include "structure.h"

map<string, string> RefSeqMap;
map<string, string> AttractMap;
map<string, map<int, int> > ClinVarMap;
map<pair<string, int>, bp_info_t> BPmap;
map<pair<string, int>, bool> HGMDmap, dbSNPmap;
map<string, map<int, AltSpliceInfo_t> > AltSpliceMap;

map<string, float> Model0, Model1, Model2, Model3, Model4, Model5;
map<string, vector<float> > PhyloP100Map, PhastCons100Map, PhyloP20Map, PhyloP7Map;

void LoadRefSeq()
{
	int p;
	fstream file;
	string str, chr;

	file.open("data/hg38_ref.fa", ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		if (str[0] == '>')
		{
			if ((p = str.find_first_of(' ')) > 0) str.resize(p);
			chr = str.substr(1); RefSeqMap[chr] = "N"; //0-base --> 1-base
			fprintf(stderr, "\33[2K\rLoad chromosome %s...", chr.c_str());
		}
		else RefSeqMap[chr] += str;
	}
	file.close();
	fprintf(stderr, "Done!\n");
}

void LoadExpBPs(const char* filename)
{
	fstream file;
	stringstream ss;
	int pos1, pos2, bp;
	string str, chr, tmp;
	map<pair<string, int>, bp_info_t>::iterator iter;

	fprintf(stderr, "Load experimental branchsites...");
	file.open(filename, ios_base::in);
	if (!file.is_open())
	{
		fprintf(stderr, "Error! %s is not accessible!\n", filename);
		exit(1);
	}
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		ss.clear(); ss.str(str); ss >> chr >> pos1 >> pos2 >> tmp >> bp;
		if ((iter = BPmap.find(make_pair(chr, pos2))) == BPmap.end()) BPmap[make_pair(chr, pos2)].bp = bp;
		else if (iter->second.bp < bp) iter->second.bp = bp;

		switch (bp)
		{
		case 0: BPmap[make_pair(chr, pos2)].bp0 = true; break;
		case -2: BPmap[make_pair(chr, pos2)].bp2 = true; break;
		case -3: BPmap[make_pair(chr, pos2)].bp3 = true; break;
		}
	}
	file.close();
	fprintf(stderr, "Get %d observed branch sites\n", (int)BPmap.size());
}

void LoadHGMD(const char* filename)
{
	fstream file;
	int i, p, pos;
	stringstream ss;
	string str, chr, tmp;
	vector<string> vec(9);

	fprintf(stderr, "Load HGMD database...");
	file.open(filename, ios_base::in); getline(file, str);
	if (!file.is_open())
	{
		fprintf(stderr, "Error! %s is not accessible!\n", filename);
		exit(1);
	}
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		ss.clear(); ss.str(str);
		for (i = 0; i < 9; i++) getline(ss, vec[i], '\t');
		if (vec[5].find("DM") != string::npos)
		{
			p = vec[8].find_first_of(':'); chr = vec[8].substr(0, p); tmp = vec[8].substr(p + 1); pos = atoi(tmp.c_str());
			HGMDmap[make_pair(chr, pos)] = true;
		}
	}
	file.close();
	fprintf(stderr, "Done!\n");
}

void LoadSNPdatabase(const char* filename)
{
	char strand;
	fstream file;
	int pos1, pos2;
	stringstream ss;
	string str, chr, tmp;
	map<pair<string, int>, bp_info_t>::iterator iter;

	fprintf(stderr, "Load SNP database...");
	file.open(filename, ios_base::in);
	if (!file.is_open())
	{
		fprintf(stderr, "Error! %s is not accessible!\n", filename);
		exit(1);
	}
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		ss.clear(); ss.str(str); ss >> chr >> pos1 >> pos2 >> tmp >> tmp >> strand;
		dbSNPmap[make_pair(chr, pos2)] = (strand == '+' ? true : false);
	}
	file.close();
	fprintf(stderr, "Done!\n");
}

void LoadClinVarMap(const char* filename)
{
	fstream file;
	int i, p1, p2;
	stringstream ss;
	vector<string> vec;
	string str, chr, tmp;

	fprintf(stderr, "Load ClinVar database...");
	file.open(filename, ios_base::in); getline(file, str);
	if (!file.is_open())
	{
		fprintf(stderr, "Error! %s is not accessible!\n", filename);
		exit(1);
	}
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		ss.clear(); ss.str(str); vec.clear(); vec.resize(24);
		for (i = 0; i < 24; i++) getline(ss, vec[i], '\t');
		if (vec[7].find("Pathogenic") != string::npos || vec[7].find("Likely pathogenic") != string::npos)
		{
			p1 = atoi(vec[20].c_str()); p2 = atoi(vec[21].c_str());
			if (p2 - p1 < 15)
			{
				if (vec[19] == "MT") vec[19] = "M"; 
				chr = "chr" + vec[19];
				//printf("%s: %d-%d\n", chr.c_str(), p1, p2);
				ClinVarMap[chr].insert(make_pair(p2, p1));
			}
		}
	}
	fprintf(stderr, "Done!\n");
}

void LoadAltSplicedRegion(const char* filename)
{
	fstream file;
	stringstream ss;
	string str, chr, tmp;
	AltSpliceInfo_t AltSpliceInfo;

	fprintf(stderr, "Load UCSC alternative spliced regions...");
	file.open(filename, ios_base::in); getline(file, str);
	if (!file.is_open())
	{
		fprintf(stderr, "Error! %s is not accessible!\n", filename);
		exit(1);
	}
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		ss.clear(); ss.str(str); ss >> chr >> AltSpliceInfo.beg >> AltSpliceInfo.end;
		AltSpliceInfo.strand = str[str.length() - 1] == '+' ? true : false;
		AltSpliceMap[chr].insert(make_pair(AltSpliceInfo.end, AltSpliceInfo));
	}
	fprintf(stderr, "Done!\n");
}

void LoadAttract_Database(const char* filename)
{
	int i;
	fstream file;
	stringstream ss;
	map<string, string> m;
	vector<string> vec(12);
	string id, gene, str, tmp;

	fprintf(stderr, "Load Attract database...");
	file.open(filename, ios_base::in); getline(file, str);
	if (!file.is_open())
	{
		fprintf(stderr, "Error! %s is not accessible!\n", filename);
		exit(1);
	}
	while (!file.eof())
	{
		getline(file, str); if (str == "") break;
		ss.clear(); ss.str(str); vec.clear();
		for (i = 0; i < 12; i++) getline(ss, vec[i], '\t');
		AttractMap.insert(make_pair(vec[11], vec[0]));
	}
	file.close();
	fprintf(stderr, "Done!\n");
}

void LoadModelParameters(const char* filename, map<string, float>& Model)
{
	int i;
	fstream file;
	stringstream ss;
	vector<string> vec(3);
	string str, tmp, factor;

	file.open(filename, ios_base::in);
	if (!file.is_open())
	{
		fprintf(stderr, "Error! %s is not accessible!\n", filename);
		exit(1);
	}
	getline(file, str); // header
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		ss.clear(); ss.str(str); for (i = 0; i < 3; i++) getline(ss, vec[i], ',');
		Model.insert(make_pair(vec[1], atof(vec[2].c_str())));
	}
	file.close();

	//for (map<string, float>::iterator iter = Model.begin(); iter != Model.end(); iter++) printf("%s: %f\n", iter->first.c_str(), iter->second);
}

void ReadPhyloP(map<string, vector<float> >& map, const char* filename) // wig file
{
	float f;
	fstream file;
	int p1, p2, pos;
	stringstream ss;
	string str, chr, tmp;
	bool bValideChr = false;

	fprintf(stderr, "Read %s...\n", filename);
	file.open(filename, ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		if (str[0] == 'f')
		{
			p1 = str.find("chrom=") + 6; p2 = str.find_first_of(' ', p1); chr = str.substr(p1, p2 - p1);
			bValideChr = chr.find_last_of('_') == string::npos ? true : false;
			p1 = str.find("start=", p2) + 6; p2 = str.find_first_of(' ', p1); pos = atoi(str.substr(p1, p2 - p1).c_str());
		}
		if (bValideChr)
		{
			f = atof(str.c_str());
			map[chr][pos++] = f;
		}
	}
	file.close();
}

void LoadBinaryPhyloP(map<string, vector<float> >& PhyloPmap, const char* filename) // binary phylop
{
	fstream file;
	fprintf(stderr, "Load %s...", filename);
	file.open(filename, ios_base::in | ios_base::binary);
	if (!file.is_open())
	{
		fprintf(stderr, "Error! %s is not accessible!\n", filename);
		exit(1);
	}
	for (map<string, vector<float> >::iterator iter = PhyloPmap.begin(); iter != PhyloPmap.end(); iter++) file.read((char*)&iter->second[0], sizeof(float) * iter->second.size());
	file.close();
	fprintf(stderr, "Done!\n");
}

void LoadPhyloPData()
{
	int chr_size;
	fstream file;
	stringstream ss;
	string str, chr;

	file.open("data/hg38_ref.fa.fai", ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		ss.clear(); ss.str(str); ss >> chr >> chr_size; chr_size += 1;
		PhyloP100Map[chr].resize(chr_size); 
		PhyloP20Map[chr].resize(chr_size);
		PhyloP7Map[chr].resize(chr_size);
		PhastCons100Map[chr].resize(chr_size);
	}
	file.close();

	LoadBinaryPhyloP(PhyloP7Map, "data/phylop7.bin");
	LoadBinaryPhyloP(PhyloP20Map, "data/phylop20.bin");
	LoadBinaryPhyloP(PhyloP100Map, "data/phylop100.bin");
	LoadBinaryPhyloP(PhastCons100Map, "data/phastcons100.bin");
}

int main(int argc, char* argv[])
{
	LoadExpBPs("data/bs_exp.bed"); // experimental break points
	LoadSNPdatabase("data/snp150Flagged_intron.bed");
	LoadHGMD("data/HGMD_intron-ss-hg38-strand.txt");
	LoadClinVarMap("data/ClinVar.txt");
	LoadAltSplicedRegion("data/GENCODEv24_intron-alt.txt");
	LoadAttract_Database("data/ATtRACT_db.txt");

	LoadPhyloPData();
	LoadModelParameters("data/5ss_noWT_coef_nostd_959_corefeature_wmm.csv", Model0);
	LoadModelParameters("data/AG_coef_nostd_footnote.csv", Model1);
	LoadModelParameters("data/nonAG_coef_nostd_footnote.csv", Model2);
	LoadModelParameters("data/AG_efficiency_coef_nostd_footnote.csv", Model3);
	LoadModelParameters("data/nonAG_efficiency_coef_nostd_footnote.csv", Model4);
	LoadModelParameters("data/wt_efficiency_coef_nostd_footnote.csv", Model5);

	LoadRefSeq();
	GetIntronRegions();
	
	LocalMode(argv[1]);

	return 0;
}
