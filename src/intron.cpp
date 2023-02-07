#include "structure.h"

map<string, map<string, Gene_t> > GeneMap;
map<string, map<int, string> > GeneRegionMap;

bool CompByPos(const exon_t& p1, const exon_t& p2)
{
	if (p1.beg == p2.beg) return p1.end < p2.end;
	else return p1.beg < p2.beg;
}

bool CompByIdentity(const intron_t& p1, const intron_t& p2)
{
	return (p1.beg == p2.beg) && (p1.end == p2.end);
}

string GenRevStrand(string seq)
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

void GetIntronRegions()
{
	int i, j;
	char strand;
	exon_t exon;
	fstream file;
	intron_t intron;
	stringstream ss;
	string str, chr, database, region, info1, info2, tmp, trantID, geneID;

	file.open("/data/Lab_data/N420/20210037/kernel/data/gencode.v37.annotation.gtf", ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "" || str[0] == '#') continue;
		ss.clear(); ss.str(str); ss >> chr >> database >> region;

		if (region == "exon") 
		{
			ss >> exon.beg >> exon.end >> tmp >> strand >> tmp >> tmp >> tmp >> tmp >> info2 >> tmp >> tmp >> tmp >> info1;
			geneID = info1.substr(1, info1.length() - 3);
			trantID = info2.substr(1, info2.length() - 3);
			exon.strand = (strand == '+' ? true : false);
			GeneMap[chr][geneID].TranscriptExonMap[trantID].push_back(exon);

			fprintf(stderr, "\33[2K\rGet gene %s...", geneID.c_str());
			if (GeneMap[chr][geneID].beg == 0 || GeneMap[chr][geneID].beg > exon.beg) GeneMap[chr][geneID].beg = exon.beg;
			if (GeneMap[chr][geneID].end == 0 || GeneMap[chr][geneID].end < exon.end) GeneMap[chr][geneID].end = exon.end;
		}
	}
	file.close();

	for (map<string, map<string, Gene_t> >::iterator ChrIter = GeneMap.begin(); ChrIter != GeneMap.end(); ChrIter++)
	{
		for (map<string, Gene_t>::iterator GeneIter = ChrIter->second.begin(); GeneIter != ChrIter->second.end(); GeneIter++)
		{
			for (map<string, vector<exon_t> >::iterator TrantIter = GeneIter->second.TranscriptExonMap.begin(); TrantIter != GeneIter->second.TranscriptExonMap.end(); TrantIter++)
			{
				sort(TrantIter->second.begin(), TrantIter->second.end(), CompByPos);
				for (i = 0, j = 1; j < (int)TrantIter->second.size(); i++, j++)
				{
					intron.beg = TrantIter->second[i].end + 1; intron.end = TrantIter->second[j].beg - 1; intron.strand = TrantIter->second[i].strand;
					GeneIter->second.TranscriptIntronMap[TrantIter->first].push_back(intron);
				}
			}
		}
	}
	for (map<string, map<string, Gene_t> >::iterator ChrIter = GeneMap.begin(); ChrIter != GeneMap.end(); ChrIter++)
	{
		for (map<string, Gene_t>::iterator GeneIter = ChrIter->second.begin(); GeneIter != ChrIter->second.end(); GeneIter++)
		{
			GeneRegionMap[ChrIter->first].insert(make_pair(GeneIter->second.beg, GeneIter->first));
		}
	}
	fprintf(stderr, "Done!\n");
}
