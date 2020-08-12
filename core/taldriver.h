// by J@:)
#ifndef TALDRIVER_H
#define TALDRIVER_H

#include "aligner.h"
#include "talmatrix.h"
#include<future>
#include<thread>

using namespace std;

vector<TAlignment> AlignDriverSE(vector<string>& readBufferChain, ReadAlignmentTask& preRAT)
{
	vector<TAlignment> Alignments;
	int minReadLength = TAlignerOptions::Options().getOption("Aligner|MinMappedSegment");
	int seedLength    = TAlignerOptions::Options().getOption("Aligner|Seed");
	minReadLength += (seedLength * 1);
	for(int i = 0; i < readBufferChain.size(); i++)
	{
		ReadAlignmentTask fwRAT = preRAT;
		ReadAlignmentTask rvRAT = preRAT;
		fwRAT.read = MakeTless(readBufferChain[i]);
		rvRAT.read = MakeTless(ReverseAndComplementDNA(readBufferChain[i]));
		TAlignment alignFw, alignRv;
		if(fwRAT.read.T.size() > minReadLength)
		{
			alignFw = TAlignRead(fwRAT);
		}
		if(rvRAT.read.T.size() > minReadLength)
		{
			alignRv = TAlignRead(rvRAT);
			alignRv.setReversed();
		}
		if(alignFw.isValid())
		{
			if(alignRv.isValid())
			{
				if(alignFw.length() >= alignRv.length())
				{
					Alignments.push_back(alignFw);
				}
				else
				{
					Alignments.push_back(alignRv);
				}
			}
			else
			{
				Alignments.push_back(alignFw);
			}
		}
		else if(alignRv.isValid())
		{
			Alignments.push_back(alignRv);
		}
	}
	return Alignments;
}

vector<TAlignment> SingleEndAlignFastqMT(const char* libFile, ReadAlignmentTask& preRAT)
{
	vector<TAlignment> Alignments;
	TAlignerOptions::Options().Log("Staring to map reads [SE mode]: " + string(libFile));
	ifstream ifile(libFile);
	string line1, line2, line3, line4;

	int totalReadReadCounter  = 0;
	int totalReadEntryIgnored = 0;

	int nThreads = thread::hardware_concurrency() - 1;
	int userTheads = TAlignerOptions::Options().getOption("T-Aligner|Threads");
	if(userTheads > 0) nThreads = userTheads;
	TAlignerOptions::Options().Log("Using " + std::to_string(nThreads) + " threads for that.");
	string memPreset = TAlignerOptions::Options().getOption("T-Aligner|Reads/job");
	int jobSizeInReads = TAlignerOptions::Options().getOption("T-Aligner|MemoryUsePreset|" + memPreset);
	TAlignerOptions::Options().Log("Job size " + std::to_string(jobSizeInReads) + " as suggested by memory preset \"" + memPreset + "\"");

	vector<vector<string> > readsBuffer(nThreads);
	vector<future<vector<TAlignment> > > alignFutures;
	const int maximalJobSize = jobSizeInReads * nThreads;
	int totalJobSizeCounter = 0;

	auto DoMapping = [&]()
    {
		// MT version

		for(int i = 0; i < readsBuffer.size(); i++)
		{
			alignFutures.push_back(async(launch::async,
				AlignDriverSE, std::ref(readsBuffer[i]), std::ref(preRAT)));
		}
		for(auto &alignResult : alignFutures)
		{
			auto result = alignResult.get();
			Alignments.insert(Alignments.end(), result.begin(), result.end());
		}
		alignFutures.clear();

		// single-threaded version
		/*
		for(int i = 0; i < readsBuffer.size(); i++)
		{
			auto result = AlignDriverSE(std::ref(readsBuffer[i]), std::ref(preRAT));
			Alignments.insert(Alignments.end(), result.begin(), result.end());
		}
		*/
		readsBuffer = vector<vector<string> >(nThreads);
	};

	while(  std::getline(ifile, line1, '\n') &&
			std::getline(ifile, line2, '\n') &&
			std::getline(ifile, line3, '\n') &&
			std::getline(ifile, line4, '\n'))
	{
		if(line1[0] == '@' && line2.size() > 1 && line2.size() == line4.size())
		{
			// check if the read's sequence is composed of only
			// legal charaters A/T/C/G.
			// T-Aligner will skip read if it contains any other charater
			if(DoReadSequenceSanityCheck(line2)) {

				totalReadReadCounter++;
				readsBuffer[totalJobSizeCounter % nThreads].push_back(line2);
				totalJobSizeCounter++;
				if(totalJobSizeCounter > maximalJobSize)
				{
					totalJobSizeCounter = 0;
					DoMapping();
				}
			}
		}
		else
		{
			totalReadEntryIgnored++;
		}
	}
	DoMapping();
	return Alignments;
}

string PrintReadTaf(TAlignment& read, vector<TlessDNA>& refTless, const string& read_id = "@tal")
{
    string taf_line_for_read = "";
    // ref id
    taf_line_for_read += std::to_string(read.ref_Idx);
    taf_line_for_read += "\t";
    // read if
    taf_line_for_read += read_id;
    taf_line_for_read += "\t";
    // ref start
    taf_line_for_read += std::to_string(read.ref_S);
    taf_line_for_read += "\t";
    // ref end
    taf_line_for_read += std::to_string(read.ref_E);
    taf_line_for_read += "\t";
    // dR
    for(auto& dr : read.dR)
    {
        taf_line_for_read += std::to_string(dr);
        taf_line_for_read += ";";
    }
    taf_line_for_read += "\t";
    // sequence
    taf_line_for_read += Alignment2Seq(read, refTless);

    return taf_line_for_read;
}

void WriteAlignedReadsTAF(vector<TAlignment>& reads,
	vector<TlessDNA>& refTless, const char* file)
{
	ofstream ofile(file);
	int readId = 0;
	for(int i = 0; i < reads.size(); i++)
	{
		ofile << reads[i].ref_Idx << "\t";
		ofile << "@tal" << ++readId << "\t";
		ofile << reads[i].ref_S << "\t";
		ofile << reads[i].ref_E << "\t";
		for(int j = 0; j < reads[i].dR.size(); j++)
		{
			ofile << reads[i].dR[j] << ";";
		}
        ofile << "\t";
        auto read_seq = Alignment2Seq(reads[i], refTless);
        ofile << read_seq;
		ofile << "\n";
	}
	ofile.close();
}

void WriteAlignedReadsFastq(vector<TAlignment>& reads,
	vector<TlessDNA>& refTless, const char* file)
{
	ofstream ofile(file);
	int readId = 0;
	for(int i = 0; i < reads.size(); i++)
	{
		ofile << "@tal" << ++readId << "\n";
		auto rSeq = Alignment2Seq(reads[i], refTless);
		int  rSz  = rSeq.size();

		ofile << rSeq << "\n";
		ofile << "+\n";
		for(int j = 0; j < rSz; j++) ofile << "J";

		ofile << "\n";

	}
	ofile.close();
}


string PrintAlignedReadFasta(TAlignment& read, vector<TlessDNA>& refTless,
    const string& reference_name = "ref", const string& read_name = "tal")
{
    string result = "";

    const string& refT = refTless[read.ref_Idx].T;
    vector<difft_t>& refDt = refTless[read.ref_Idx].dT;

    string alignedRefSeq = "";
    string alignedReadSeq = "";
    int curRefPos = 0;
    while(curRefPos < read.ref_S)
    {
        alignedReadSeq += "-";
        alignedRefSeq  += refT.substr(curRefPos, 1);
        difft_t tacc = refDt[curRefPos+1];
        for(difft_t t = 0; t < tacc; t++)
        {
            alignedReadSeq += "-";
            alignedRefSeq  += "T";
        }
        curRefPos++;
    }
    int curReadPos = 0;
    while(curRefPos < read.ref_E - 1)
    {
        alignedReadSeq += refT.substr(curRefPos, 1);
        alignedRefSeq  += refT.substr(curRefPos, 1);
        difft_t readTs = (read.dR[curReadPos] + refDt[curRefPos+1]);
        difft_t refTs  = refDt[curRefPos+1];
        if(refTs > readTs)
        {
            for(int t = 0; t < readTs; t++)
            {
                alignedReadSeq += "T";
                alignedRefSeq  += "T";
            }
            for(int t = 0; t < refTs-readTs; t++)
            {
                alignedReadSeq += "-";
                alignedRefSeq  += "T";
            }
        }
        else
        {
            for(int t = 0; t < refTs; t++)
            {
                alignedReadSeq += "T";
                alignedRefSeq  += "T";
            }
            for(int t = 0; t < readTs-refTs; t++)
            {
                alignedReadSeq += "t";
                alignedRefSeq  += "-";
            }
        }

        curRefPos++;
        curReadPos++;
    }


    while(curRefPos < refT.size())
    {
        alignedReadSeq += "-";
        alignedRefSeq  += refT.substr(curRefPos, 1);
        difft_t tacc = refDt[curRefPos+1];
        for(difft_t t = 0; t < tacc; t++)
        {
            alignedReadSeq += "-";
            alignedRefSeq  += "T";
        }
        curRefPos++;
    }
    result += (">" + reference_name + "\n");
    result += alignedRefSeq;
    result += "\n";
    result += (">" + read_name + "\n");
    result += alignedReadSeq;
    result += "\n";

    return result;
}


void WriteAlignedReadFasta(TAlignment& read, vector<TlessDNA>& refTless,
    const char* file)
{
    ofstream ofile(file);

    string result =
        PrintAlignedReadFasta(read, refTless);

    ofile << result;

    ofile.close();
}

#ifdef DOUNITTEST

void Test3()
{
	//const char* referenceFile = "/Users/jalgard/rps12.fa";
    //const char* fastqFile     = "/Users/jalgard/Lsey23_2.connect.fq";

	//const char* referenceFile = "../tst/rps12.fa";
    //const char* fastqFile     = "../tst/test.fastq";

	const char* referenceFile = "/Users/jalgard/nd8.fa";
    const char* fastqFile     = "/Users/jalgard/polyAmix_merged_filtered.fastq";

	vector<vector<string> > refs;
    vector<TlessDNA> refTless;
    unordered_map<string, vector<Rindex> > refIndex;

    ProcessReferenceFile(referenceFile, refs, refTless);
    Build_Fasta_Index(refTless, refIndex);

	ReadAlignmentTask RAT;
	RAT.refIndex = &refIndex;
	RAT.refTless = &refTless;

	auto All = SingleEndAlignFastqMT(fastqFile, RAT);
	return;
	int matchingRef = 0;
	for(int i = 0; i < All.size(); i++)
	{
		if(All[i].doesMatchRef()) matchingRef++;
	}
	//cout << "Aligned " << All.size() << ", match ref " << matchingRef << "\n";

	/*
	auto dedupAl = FilterMappingDuplicatesReads(All);
	for(int i = 0; i < 10; i++)
	{
		cout << dedupAl[i].ref_S << ", " << dedupAl[i].ref_E << "\n";
		//for(int j = 0; j < dedupAl[i].dR.size(); j++)
		//	cout << static_cast<char>(dedupAl[i].dR[j]) << " ";
		//cout << "\n";
		cout << Alignment2Seq(dedupAl[i], refTless) << "\n";
	}
	*/
	//cout << "Total " << All.size() << ", dedup " << dedupAl.size() << "\n";
	vector<vector<TAlignment> > SplittedRC = AlignmentSplitterByRevcom(All);
	//cout << All.size() << "\n\n";
	vector<vector<int> > rMfw = CalcTalMatrix(SplittedRC[0], RAT, 0);
	vector<vector<int> > rMrc = CalcTalMatrix(SplittedRC[1], RAT, 0);

	for(int i = 0; i < 40; i++)
	{
		for(int j = 0; j < rMfw.size(); j++)
		{
			cout << "\t" << rMfw[j][i];
		}
		cout << "\n";
	}

	for(int i = 0; i < 40; i++)
	{
		for(int j = 0; j < rMrc.size(); j++)
		{
			cout << "\t" << rMrc[j][i];
		}
		cout << "\n";
	}

}

#endif
#endif
