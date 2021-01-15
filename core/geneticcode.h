// by J@:)
#ifndef GENETICCODE_H
#define GENETICCODE_H
#include <map>
#include <string>

using namespace std;

static inline map<string, string> setGeneticCode(string geneticCodeId)
{
	map<string, string> geneticCodeTable;
	if(geneticCodeId == "T-Aligner-default")
	{

		// sense codon assignments
		geneticCodeTable["TTT"]="F";	geneticCodeTable["TCT"]="S";	geneticCodeTable["TAT"]="Y";	geneticCodeTable["TGT"]="C";
		geneticCodeTable["TTC"]="F";	geneticCodeTable["TCC"]="S";	geneticCodeTable["TAC"]="Y";	geneticCodeTable["TGC"]="C";
		geneticCodeTable["TTA"]="L";	geneticCodeTable["TCA"]="S";	geneticCodeTable["TAA"]="*";	geneticCodeTable["TGA"]="W";
		geneticCodeTable["TTG"]="L";	geneticCodeTable["TCG"]="S";	geneticCodeTable["TAG"]="*";	geneticCodeTable["TGG"]="W";
		geneticCodeTable["CTT"]="L";	geneticCodeTable["CCT"]="P";	geneticCodeTable["CAT"]="H";	geneticCodeTable["CGT"]="R";
		geneticCodeTable["CTC"]="L";	geneticCodeTable["CCC"]="P";	geneticCodeTable["CAC"]="H";	geneticCodeTable["CGC"]="R";
		geneticCodeTable["CTA"]="L";	geneticCodeTable["CCA"]="P";	geneticCodeTable["CAA"]="Q";	geneticCodeTable["CGA"]="R";
		geneticCodeTable["CTG"]="L";	geneticCodeTable["CCG"]="P";	geneticCodeTable["CAG"]="Q";	geneticCodeTable["CGG"]="R";
		geneticCodeTable["ATT"]="I";	geneticCodeTable["ACT"]="T";	geneticCodeTable["AAT"]="N";	geneticCodeTable["AGT"]="S";
		geneticCodeTable["ATC"]="I";	geneticCodeTable["ACC"]="T";	geneticCodeTable["AAC"]="N";	geneticCodeTable["AGC"]="S";
		geneticCodeTable["ATA"]="I";	geneticCodeTable["ACA"]="T";	geneticCodeTable["AAA"]="K";	geneticCodeTable["AGA"]="R";
		geneticCodeTable["ATG"]="M";	geneticCodeTable["ACG"]="T";	geneticCodeTable["AAG"]="K";	geneticCodeTable["AGG"]="R";
		geneticCodeTable["GTT"]="V";	geneticCodeTable["GCT"]="A";	geneticCodeTable["GAT"]="D";	geneticCodeTable["GGT"]="G";
		geneticCodeTable["GTC"]="V";	geneticCodeTable["GCC"]="A";	geneticCodeTable["GAC"]="D";	geneticCodeTable["GGC"]="G";
		geneticCodeTable["GTA"]="V";	geneticCodeTable["GCA"]="A";	geneticCodeTable["GAA"]="E";	geneticCodeTable["GGA"]="G";
		geneticCodeTable["GTG"]="V";	geneticCodeTable["GCG"]="A";	geneticCodeTable["GAG"]="E";	geneticCodeTable["GGG"]="G";

		// nonsense codon assignments
		geneticCodeTable["->ATG"] = "start";
		geneticCodeTable["->TTG"] = "start";
        geneticCodeTable["->GTG"] = "start";
        geneticCodeTable["->CTG"] = "start";
        geneticCodeTable["->ATA"] = "start";
		geneticCodeTable["<-TAA"] = "stop";
		geneticCodeTable["<-TAG"] = "stop";

	}

    if(geneticCodeId == "T-Aligner-default-ATG")
    {

        // sense codon assignments
        geneticCodeTable["TTT"]="F";	geneticCodeTable["TCT"]="S";	geneticCodeTable["TAT"]="Y";	geneticCodeTable["TGT"]="C";
        geneticCodeTable["TTC"]="F";	geneticCodeTable["TCC"]="S";	geneticCodeTable["TAC"]="Y";	geneticCodeTable["TGC"]="C";
        geneticCodeTable["TTA"]="L";	geneticCodeTable["TCA"]="S";	geneticCodeTable["TAA"]="*";	geneticCodeTable["TGA"]="W";
        geneticCodeTable["TTG"]="L";	geneticCodeTable["TCG"]="S";	geneticCodeTable["TAG"]="*";	geneticCodeTable["TGG"]="W";
        geneticCodeTable["CTT"]="L";	geneticCodeTable["CCT"]="P";	geneticCodeTable["CAT"]="H";	geneticCodeTable["CGT"]="R";
        geneticCodeTable["CTC"]="L";	geneticCodeTable["CCC"]="P";	geneticCodeTable["CAC"]="H";	geneticCodeTable["CGC"]="R";
        geneticCodeTable["CTA"]="L";	geneticCodeTable["CCA"]="P";	geneticCodeTable["CAA"]="Q";	geneticCodeTable["CGA"]="R";
        geneticCodeTable["CTG"]="L";	geneticCodeTable["CCG"]="P";	geneticCodeTable["CAG"]="Q";	geneticCodeTable["CGG"]="R";
        geneticCodeTable["ATT"]="I";	geneticCodeTable["ACT"]="T";	geneticCodeTable["AAT"]="N";	geneticCodeTable["AGT"]="S";
        geneticCodeTable["ATC"]="I";	geneticCodeTable["ACC"]="T";	geneticCodeTable["AAC"]="N";	geneticCodeTable["AGC"]="S";
        geneticCodeTable["ATA"]="I";	geneticCodeTable["ACA"]="T";	geneticCodeTable["AAA"]="K";	geneticCodeTable["AGA"]="R";
        geneticCodeTable["ATG"]="M";	geneticCodeTable["ACG"]="T";	geneticCodeTable["AAG"]="K";	geneticCodeTable["AGG"]="R";
        geneticCodeTable["GTT"]="V";	geneticCodeTable["GCT"]="A";	geneticCodeTable["GAT"]="D";	geneticCodeTable["GGT"]="G";
        geneticCodeTable["GTC"]="V";	geneticCodeTable["GCC"]="A";	geneticCodeTable["GAC"]="D";	geneticCodeTable["GGC"]="G";
        geneticCodeTable["GTA"]="V";	geneticCodeTable["GCA"]="A";	geneticCodeTable["GAA"]="E";	geneticCodeTable["GGA"]="G";
        geneticCodeTable["GTG"]="V";	geneticCodeTable["GCG"]="A";	geneticCodeTable["GAG"]="E";	geneticCodeTable["GGG"]="G";

        // nonsense codon assignments
        geneticCodeTable["->ATG"] = "start";
        //geneticCodeTable["->TTG"] = "start";
        //geneticCodeTable["->GTG"] = "start";
        //geneticCodeTable["->ATA"] = "start";
        geneticCodeTable["<-TAA"] = "stop";
        geneticCodeTable["<-TAG"] = "stop";

    }
	return geneticCodeTable;
}

static inline map<string, double> setCodonNs(map<string, string>& geneticCode)
{
    map<string, double> codonNs;
    static const string dnaAlphabet = "AGTC";
    for(auto& codon_aa : geneticCode)
    {
        string codon = codon_aa.first;
        double sN = 0;
        for(int i = 0; i < 3; i++)
        {
            double sn = 0;
            for(const char& letter : dnaAlphabet)
            {
                if(letter == codon[i]) continue;
                string mut = codon;
                mut[i] = letter;
                if(geneticCode[mut] != geneticCode[codon])
                {
                    sn += 1.0;
                }
            }
            sn = sn/3;
            sN += sn;
        }
        codonNs[codon] = sN;
    }
    return codonNs;
}

static inline string Translate(const std::string& sequence, map<string, string>& geneticCodeTable,
	int frame=0)
{
    string translation = "";

    for(int i = frame; i < sequence.size()-3; i+=3)
        translation += geneticCodeTable[sequence.substr(i,3)];

    return translation;
}

#endif
