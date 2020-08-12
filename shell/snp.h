// by J@:)
#ifndef SHELL_SNP_H
#define SHELL_SNP_H

#include "taldriver.h"
#include <fstream>

using namespace std;


class SnpData
{
public:
    string chr;
    int pos;
    string ref;
    string alt;
};

vector<SnpData> LoadSnpSet(const char* filename)
{
    vector<SnpData> SNPs;

    ifstream ifile(filename); string line;
    while(getline(ifile, line, '\n'))
    {
        auto fields = split(line, '\t');
        SnpData data;
        data.chr = fields[0];
        data.pos = atoi(fields[1].c_str());
        data.ref = fields[2];
        data.alt = fields[3];
        SNPs.push_back(data);
    }
    sort(SNPs.begin(), SNPs.end(), [](const SnpData& a, const SnpData& b)
        -> bool { return a.pos < b.pos; });
    sort(SNPs.begin(), SNPs.end(), [](const SnpData& a, const SnpData& b)
        -> bool { return a.chr < b.chr; });
    return SNPs;
}

void CastCoordinates(vector<SnpData>& SNPs, int offset)
{
    for(auto& snp : SNPs)
    {
        snp.pos -= offset;
    }
}
#endif
