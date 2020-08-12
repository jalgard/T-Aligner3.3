#include "taldriver.h"
#include <fstream>

int AlignLibOnRef(string ref, string mrna)
{

    string ref_file_name = split(ref, '/').back();
    cerr << "Starting T-Aligner protocol:\t\"Alignment of mRNA to single reference\"\n";

	vector<vector<string> > refsBank, mrnaBank;
	vector<TlessDNA> refTless, mrnaTless;
	unordered_map<string, vector<Rindex> > refIndex;

	ProcessReferenceFile(ref.c_str(), refsBank, refTless);
    ProcessReferenceFile(mrna.c_str(), mrnaBank, mrnaTless);
	Build_Fasta_Index(refTless, refIndex);

    cerr << "Checking if reference fasta is good:\t";
    if(refsBank.size() == 1) cerr << "OK!\n";
    else if(refsBank.size() > 1) cerr << "WARNING! Two or more references found!\n";
    else { cerr << "ERROR! Reference .fasta is not properly formatted.\n"; return 1; }

	ReadAlignmentTask gRAT;
	gRAT.refIndex = &refIndex;
	gRAT.refTless = &refTless;

    cerr << "Starting alignment process...\n";

    for(int mId = 0; mId < mrnaBank.size(); mId++)
    {
        ReadAlignmentTask fwRAT = gRAT;
        fwRAT.read = mrnaTless[mId];
        TAlignment alignFw = TAlignRead(fwRAT);
        string alignResult = "Unmapped";
        if(alignFw.isValid())
        {
            alignResult = "Mapped";
            vector<TAlignment> alignVec = {alignFw};
            string taf_output_file   = ref_file_name + "_" + mrnaBank[mId][0] + "_TA3.taf";
            string fasta_output_file = ref_file_name + "_" + mrnaBank[mId][0] + "_TA3.fasta";

        	WriteAlignedReadsTAF(alignVec, refTless, taf_output_file.c_str());
            WriteAlignedReadFasta(alignFw, refTless, fasta_output_file.c_str());
        }
        cerr << "@TAL:\t" << refsBank[0][0] << "\t" << mrnaBank[mId][0] << "\t" << alignResult << "\n";
    }

    cerr << "Protocol finished!\n";

    return 0;
}


int main(int argc, char** argv)
{
    return AlignLibOnRef(argv[1], argv[2]);
}
