#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cstring>
#include <map>
#include <unistd.h>

using namespace std;

map<string, vector<float> > PhyloP100Map, PhastCons100Map;
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

	file.open("/home/bbsc/ReferenceGenomes/hg38_lncRNA/ref.fa.fai", ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		ss.clear(); ss.str(str); ss >> chr >> chr_size; chr_size += 1;
		PhyloP100Map[chr].resize(chr_size);
		PhastCons100Map[chr].resize(chr_size);
	}
	file.close();

	LoadBinaryPhyloP(PhyloP100Map, "/home/bbsc/lab_data/N420/20210037/kernel/data/phylop100.bin");
	LoadBinaryPhyloP(PhastCons100Map, "/home/bbsc/lab_data/N420/20210037/kernel/data/phastcons100.bin");
}

int main(int argc, char* argv[])
{
	float f, sum;
	fstream file;
	stringstream ss;
	int p, beg, end, n;
	string str, chr, tmp;

	LoadPhyloPData();

	file.open(argv[1], ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		ss.clear(); ss.str(str); ss >> chr >> beg >> end;
		if (end - beg == 1) printf("%s\t%.3f\t%.3f\n", str.c_str(), PhyloP100Map[chr][end], PhastCons100Map[chr][end]);
		else
		{
			sum = 0;
			for (p = beg + 1; p <= end; p++) sum += PhyloP100Map[chr][p];
			printf("%s\t%.3f\n", str.c_str(), sum / (end - beg));
		}
	}
	file.close();

	return 0;
}
