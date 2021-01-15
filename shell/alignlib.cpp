#include "taldriver.h"
#include <fstream>

int AlignLibOnRef()
{
    cerr << "Starting T-Aligner protocol:\t\"Alignment of single library to single reference\"\n";

	vector<vector<string> > refsBank;
	vector<TlessDNA> refTless;
	unordered_map<string, vector<Rindex> > refIndex;

    string reference_filename = TAlignerOptions::Options().getOption("Input|Reference|Filename");
    string fastqInputFile = TAlignerOptions::Options().getOption("Input|Library1|Filename");

	ProcessReferenceFile(reference_filename.c_str(), refsBank, refTless);
	Build_Fasta_Index(refTless, refIndex);

    cerr << "Checking if reference fasta is good:\t";
    if(refsBank.size() == 1) cerr << "OK!\n";
    else if(refsBank.size() > 1) cerr << "WARNING! Two or more references found!\n";
    else { cerr << "ERROR! Reference .fasta is not properly formatted.\n"; return 1; }

	ReadAlignmentTask gRAT;
	gRAT.refIndex = &refIndex;
	gRAT.refTless = &refTless;
	vector<vector<TAlignment> > alignedBank;

    cerr << "Starting alignment process...\n";

    auto Alignments = SingleEndAlignFastqMT(fastqInputFile.c_str(), gRAT);

    string output_prefix = TAlignerOptions::Options().getOption("Output|Prefix");
    string output_matrix = TAlignerOptions::Options().getOption("Output|Files|ReadAlignmentsMatrix");
    string matrix_output_file = output_prefix + output_matrix;

    cerr << "Saving alignment matrix to:\t" << matrix_output_file << "\n";
	auto selectedRefAlignments = AlignmentSplitterByRefIdx(Alignments, 0);
	vector<vector<int> > Mx = CalcTalMatrix(selectedRefAlignments[0], gRAT, 0);

    ofstream ofile(matrix_output_file.c_str());

    for(int i = 0; i < 40; i++)
	{
		for(int j = 0; j < Mx.size(); j++)
		{
			ofile << "\t" << Mx[j][i];
		}
		ofile << "\n";
	}
	ofile.close();

    string output_fastq = TAlignerOptions::Options().getOption("Output|Files|MappedReadsFastq");
    string output_taf = TAlignerOptions::Options().getOption("Output|Files|MappedReadsTaf");
    string taf_output_file   = output_prefix + output_taf;
    string fastq_output_file = output_prefix + output_fastq;

    cerr << "Saving read alignments in TAF format:\t" << taf_output_file << "\n";
	WriteAlignedReadsTAF(selectedRefAlignments[0], refTless, refsBank,
            taf_output_file.c_str());

    cerr << "Saving read alignments in FASTQ format:\t" << fastq_output_file << "\n";
	WriteAlignedReadsFastq(selectedRefAlignments[0], refTless,
			fastq_output_file.c_str());

    cerr << "Protocol finished!\n";
    return 0;
}


int main(int argc, char** argv)
{
    TAlignerOptions::Options().StartLogging();
    ParseFromCommandLine(argc, argv);
    int result = AlignLibOnRef();
    return result;
}
