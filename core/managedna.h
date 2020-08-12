// by J@:)
#ifndef MANAGEDNA_H
#define MANAGEDNA_H
#include <algorithm>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

struct Rindex
{
    int p;
    int i;
    Rindex() : p(0), i(-1) {}
};

inline string UpperDNA(const string& DNA)
{
    string upperDNA = DNA;
    transform(upperDNA.begin(), upperDNA.end(), upperDNA.begin(), ::toupper);
    return upperDNA;
}

inline string ReverseAndComplementDNA(const string& DNA)
{
    string revComDNA = DNA;
    int j = 0;
    for(int i = DNA.size()-1; i >= 0; i--)
    {
        switch(DNA[i]) {
            case 'A': revComDNA[j] = 'T'; break;
            case 'G': revComDNA[j] = 'C'; break;
            case 'T': revComDNA[j] = 'A'; break;
            case 'C': revComDNA[j] = 'G'; break;
            default:  revComDNA[j] = 'N'; break; // T-Aligner will not handle Ns!
        }
        j++;
    }
    return revComDNA;
}

inline bool DoReadSequenceSanityCheck(const string& Read)
{
    for(int i = 0; i < Read.size(); i++)
    {
        char letter = Read[i];
        if(letter != 'A' && letter != 'C' && letter != 'G' && letter != 'T')
        {
            return false;
        }
    }
    return true;
}

inline vector<string> split(string str, char dlm)
{
    stringstream ss(str);
    string item;
    vector<string> splitted;
    while (std::getline(ss, item, dlm))
    {
       splitted.push_back(item);
    }
    return splitted;
}

// auxilary function to get genomic interval
// using 's' and 'e' 1-based coordinates, inclusive
static inline string subseq(const string& sequence, int s, int e)
{
    return sequence.substr(s-1, e-s+1);
}

static inline vector<vector<string> > FastaReader(const char* filename, bool split_space = true)
{
    vector<vector<string> > entries;

    string line;
    ifstream ifile(filename);
    while(getline(ifile, line, '\n'))
    {
        if(line[line.size()-1] == '\r')
        {
            line.pop_back();
        }
        if(line[0] == '>')
        {
            string name = line.substr(1);
            if(split_space)
            {
                name = split(line.substr(1), ' ')[0];
            }
            entries.push_back({name, ""});
        }
        else
        {
            entries.back()[1] += UpperDNA(line);
        }
    }
    return entries;
}

class GffReader
{

public:

    struct GffData
    {
        map<string, string> tags;
        int s;
        int e;
        char strand;
        string auth;
        string type;
        string chr;
    };

    void Load(const char* filename)
    {
        string line;
        gffdata.clear();

        ifstream ifile(filename);
        while(getline(ifile, line, '\n'))
        {
            if(line.substr(0,1) == "#") continue;

            auto fields = split(line, '\t');
            GffData gffitem;
            gffitem.s = stoi(fields[3].c_str());
            gffitem.e = stoi(fields[4].c_str());
            gffitem.strand = fields[6][0];
            gffitem.type   = fields[2];
            gffitem.chr    = fields[0];
            gffitem.auth   = fields[1];

            auto tags = split(fields[8],';');
            for(auto& tag : tags)
            {
                auto tag_p = split(tag, '=');
                gffitem.tags[tag_p[0]] = tag_p[1];
            }

            gffdata.push_back(gffitem);
        }
    }

    GffData Get(int index)
    {
        return gffdata[index];
    }

    GffData Get(string id)
    {
        for(auto& item : gffdata)
        {
            for(auto& tag : item.tags)
            {
                if(tag.first == "ID" && tag.second == id) return item;
            }
        }
        throw runtime_error("\nNo such item " + id + " in gff dataset!\n\n");
    }

    string GetSeq(int index, const string& ref_upper)
    {
        string seq = subseq(ref_upper, gffdata[index].s, gffdata[index].e);
        if(gffdata[index].strand == '-')
        {
            return ReverseAndComplementDNA(seq);
        }
        else
        {
            return seq;
        }
    }

private:
    vector<GffData> gffdata;

};

#endif
