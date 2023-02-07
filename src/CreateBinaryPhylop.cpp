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


void process_mem_usage(double& vm_usage, double& resident_set)
{
	vm_usage = 0.0;
	resident_set = 0.0;

	// the two fields we want
	unsigned long vsize;
	long rss;
	{
		std::string ignore;
		std::ifstream ifs("/proc/self/stat", std::ios_base::in);
		ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
			>> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
			>> ignore >> ignore >> vsize >> rss;
	}

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}

void GetFactorID(const char* filename, map<string, int>& m)
{
	fstream file;
	stringstream ss;
	string str, tmp, factor;

	file.open(filename, ios_base::in); getline(file, str);
	while (!file.eof())
	{
		getline(file, str); if (str == "") break;
		ss.clear(); ss.str(str); getline(ss, tmp, ','); getline(ss, factor, ',');
		m[factor]++;
	}
	file.close();
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
		else if (bValideChr)
		{
			f = atof(str.c_str());
			map[chr][pos++] = f;
		}
	}
	file.close();
}

void CreateBinaryPhyloP(const char* wig_filename, const char* bin_filename)
{
	int chr_size;
	fstream file;
	stringstream ss;
	string str, chr;
	map<string, vector<float> > m;

	file.open("/home/bbsc/ReferenceGenomes/hg38_lncRNA/ref.fa.fai", ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		ss.clear(); ss.str(str); ss >> chr >> chr_size; chr_size += 1;
		m[chr].resize(chr_size);
	}
	file.close();

	ReadPhyloP(m, wig_filename);
	file.clear(); file.open(bin_filename, ios_base::out | ios_base::binary);
	for (map<string, vector<float> >::iterator iter = m.begin(); iter != m.end(); iter++) file.write((char*)&iter->second[0], sizeof(float) * iter->second.size());
	file.close();
}

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "usage: %s wig_file out_bin_file\n", argv[0]);
		exit(0);
	}
	CreateBinaryPhyloP(argv[1], argv[2]);

	return 0;
}
