// by J@:)
#ifndef INDEXREF_H
#define INDEXREF_H

#include "tless.h"
#include "taoptions.h"
#include "managedna.h"
#include <unordered_map>

using namespace std;

int ReferenceLoader(const char* filename, vector<vector<string> >& ReferenceHolder)
{
    ifstream ifile(filename); string line;
    int fasta_counter = 0;
    while(getline(ifile, line, '\n'))
    {
        if(line[line.size() - 1] == '\r')
            line = line.substr(0, line.size() - 1);
        if(line[0] == '>')
        {
            ReferenceHolder.push_back({line.substr(1), ""});
            fasta_counter++;
        }
        else
        {
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            for(int i = 0; i < line.size(); i++)
            {
                char letter = line[i];
                if(letter != 'A' && letter != 'C' && letter != 'G' && letter != 'T')
                {
                    throw runtime_error("\nInput error: reference " +
                    ReferenceHolder.back()[0] + " from file " +
                    string(filename) + " contains " + string(letter) + " in position " +
                    std::to_string(ReferenceHolder.back()[1].size()+i+1) +
                    " which is not A/T/G/C!\n Exiting...\n\n");
                    return -1;
                }
            }
            ReferenceHolder.back()[1] += line;
        }
    }
    return fasta_counter;
}

void ProcessReferenceFile(const char* filename, vector<vector<string> >& ReferenceHolder,
    vector<TlessDNA>& ReferenceTlpHolder)
{

    TAlignerOptions::Options().Log("Loading references from file: " + string(filename));
    int operation_result = ReferenceLoader(filename, ReferenceHolder);
    TAlignerOptions::Options().Log("Total reference sequences processed: " + std::to_string(operation_result));
    if(operation_result > 0)
    {
        for(int i = ReferenceHolder.size() - operation_result; i < ReferenceHolder.size(); i++)
        {
            ReferenceTlpHolder.push_back(MakeTless(ReferenceHolder[i][1]));
        }
    }
    else
    {
        TAlignerOptions::Options().Log("WARNING! ReferenceHolder warns that there are 0 sequences in file: " + string(filename));
    }
    TAlignerOptions::Options().Log("Finished with reference file: " + string(filename));
}


// Build hash from collection of Texts
static void Build_Fasta_Index(vector<TlessDNA>& Texts, unordered_map<string, vector<Rindex> >& Index)
{
    int seedLen  = TAlignerOptions::Options().getOption("Aligner|Seed");
    for(size_t T = 0; T < Texts.size(); T++)
    {
        for(size_t p = 0; p < Texts[T].T.size() - seedLen; p++) // extract seed with minimal possible step 1
        {
            Rindex ri; ri.p = p; ri.i = T;
            Index[Texts[T].T.substr(p, seedLen)].push_back(ri);
        }
    }
}

// Look for a Seed match in Index
// S should be size of seed (8)
static vector<Rindex> Match_Index(const string& S, unordered_map<string, vector<Rindex> >& Index)
{
    auto hit = Index.find(S);
    if(hit != Index.end())
    {
        return hit->second;
    }
    return vector<Rindex>(1);
}

#endif
