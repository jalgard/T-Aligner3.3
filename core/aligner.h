// by J@:)
#ifndef ALIGNER_H
#define ALIGNER_H

#include "indexref.h"

using namespace std;

// namespace TA for all global constants

namespace TA {
    enum AlignmentFlag {
        AF_Valid             = 1 << 0,
        AF_Paired            = 1 << 1,
        AF_MateAligned       = 1 << 2,
        AF_ExactMatch        = 1 << 3,
        AF_MatchReference    = 1 << 4,
        AF_Multimapper       = 1 << 5,
        AF_WasReversed = 1 << 6,
        AF_MatesOverlap = 1 << 7,
        AF_Reserved4 = 1 << 8,
        AF_Reserved5 = 1 << 9,
        AF_Reserved6 = 1 << 10,
        AF_Reserved7 = 1 << 11,
        AF_Reserved8 = 1 << 12,
        /*
            Simple rules:
            A & B           check
            A |= B          set ON
            A &= ~B         set OFF
            A ^= B          toggle
            A |= (B1 | B2)  set ON for two opts at once
        */
    };
};



struct TAlignment
{
    vector<difft_t> dR;
    int ref_S;
    int ref_E;
    int ref_Idx;

    int mmcount;
    int ref_ESD;
    int rd_id;
    unsigned int flags;

    TAlignment() : ref_Idx(-1), flags(0) {}

    bool asRef()    { return flags & TA::AF_ExactMatch; }
    void setAsRef() { flags |=  TA::AF_ExactMatch; }

    bool isValid()  { return flags & TA::AF_Valid; }
    void setVaild() { flags |= TA::AF_Valid; }

    bool isReversed() { return flags & TA::AF_WasReversed; }
    void setReversed() { flags |= TA::AF_WasReversed; }

    bool doesMatchRef() { return flags & TA::AF_MatchReference; }

    int length() { return ref_E - ref_S; }
};

static inline string Alignment2Seq(TAlignment& alignment, vector<TlessDNA>& refTless)
{
    string alignedSeq = "";
    const string& refT = refTless[alignment.ref_Idx].T;
    vector<difft_t>& refDt = refTless[alignment.ref_Idx].dT;
    int j = 0;
    for(int i = alignment.ref_S; i < alignment.ref_E-1; i++) // correction from "i < alignment.ref_E"
    {
        alignedSeq += refT.substr(i,1);
        difft_t tacc = (alignment.dR[j] + refDt[i+1]);

        for(difft_t t = 0; t < tacc; t++)
        {
            alignedSeq += "T";
        }

        j++;
    }
    alignedSeq += refT.substr(alignment.ref_E, 1);
    return alignedSeq;
}

struct Subalignment
{
    int rd_leftmost;
    int rd_rightmost;
    int rf_start;
    int ref;
};

class ReadAlignmentTask
{
public:

    ReadAlignmentTask()
    {
        mismatch_max = TAlignerOptions::Options().getOption("Aligner|Mismatch");
        seed_length  = TAlignerOptions::Options().getOption("Aligner|Seed");
        aligned_segment_min = TAlignerOptions::Options().getOption("Aligner|MinMappedSegment");
        aligned_fraction_min = TAlignerOptions::Options().getOption("Aligner|MinMapped%");
    }
    bool init;
    int mismatch_max;
    int seed_length;
    int aligned_segment_min;
    int aligned_fraction_min;
    unordered_map<string, vector<Rindex> >* refIndex;
    vector<TlessDNA>* refTless;
    TlessDNA read;
    int rd_id;
};

TAlignment TAlignRead(ReadAlignmentTask& T)
{
    vector<Subalignment> Alignments;
    const string RdTlp = T.read.T;
    int seedLen  = T.seed_length;
    int seedPerRead = RdTlp.size() / seedLen;

    for(int seedIter = 0; seedIter < seedPerRead - 1; seedIter++)
    {
        int pos = seedIter * seedLen;
        auto hits = Match_Index(RdTlp.substr(pos, seedLen), *(T.refIndex));
        for(auto& hit : hits)
        {
            if(hit.i == -1) continue;
            bool extendsExistingHit = false;
            for(auto& align : Alignments)
            {
                //if(hit.i == align.ref && hit.p == align.rd_rightmost + seedLen)

                // 19oct2019: CRITICAL FIX on this code!
                if(hit.i == align.ref &&
                   hit.p == align.rf_start + (align.rd_rightmost - align.rd_leftmost))
                {
                    extendsExistingHit = true;
                    align.rd_rightmost += seedLen;
                    //break;
                }
            }
            if(!extendsExistingHit)
            {
                Subalignment newAlignment;
                newAlignment.rd_leftmost = pos;
                newAlignment.rd_rightmost = pos + seedLen;
                newAlignment.rf_start = hit.p;
                newAlignment.ref = hit.i;
                Alignments.push_back(newAlignment);
            }
        }
    }

    TAlignment bestAlignment;
    bestAlignment.flags &= ~TA::AF_Valid;
    bestAlignment.rd_id = T.rd_id;

    if(Alignments.size() > 0)
    {
        int alignment_probe = 0;
        RepickAlignment:

        sort(Alignments.begin(), Alignments.end(), [&](Subalignment a, Subalignment b)
    		{ return (a.rd_rightmost-a.rd_leftmost) > (b.rd_rightmost-b.rd_leftmost); });

        Subalignment& longestMatch = Alignments[alignment_probe];

        longestMatch.rd_rightmost--;
        const string RfTlp = (*(T.refTless))[longestMatch.ref].T;
        // Extend exact match
        // Extend right till first mismatch found

        while(
              longestMatch.rf_start +
              longestMatch.rd_rightmost -
              longestMatch.rd_leftmost < RfTlp.size() - 1
              &&
              longestMatch.rd_rightmost < RdTlp.size() - 1
             )
        {
            if(RfTlp[longestMatch.rf_start + longestMatch.rd_rightmost -
               longestMatch.rd_leftmost+1] != RdTlp[longestMatch.rd_rightmost+1])
            {
                break;
            }
            longestMatch.rd_rightmost++;
        }

        int mmmax = T.mismatch_max;
        if(mmmax > 2) mmmax = 2;
        int mm = 0;
        // Simple extension with up-to-2-algorithm
        while(longestMatch.rf_start > 0 &&
              longestMatch.rd_leftmost > 0 &&
              mm < mmmax
             )
        {
            if(RfTlp[longestMatch.rf_start] !=
                RdTlp[longestMatch.rd_leftmost])
                {
                    mm++;
                }
            longestMatch.rf_start--;
            longestMatch.rd_leftmost--;
        }
        while(
             longestMatch.rf_start
             + longestMatch.rd_rightmost
             - longestMatch.rd_leftmost < RfTlp.size() - 1
             &&
             longestMatch.rd_rightmost  < RdTlp.size() - 1
             &&
             mm < mmmax
             )
        {
            if(RfTlp[longestMatch.rf_start
            + longestMatch.rd_rightmost
            - longestMatch.rd_leftmost] !=
            RdTlp[longestMatch.rd_rightmost])
            {
                mm++;
            }
            longestMatch.rd_rightmost++;
        }

        // Filters
        bool alignPassFilters = true;
        int alignmentLen = longestMatch.rd_rightmost - longestMatch.rd_leftmost; // remove '-1'
        const int alignSegmentMin = T.aligned_segment_min;
        if(alignmentLen <= alignSegmentMin)
        { // DROP if minimal mapped segment is less then some absolute value
            alignPassFilters = false;
        }
        const int alignFractionMin = T.aligned_fraction_min;
        if(alignmentLen <= (1.0*alignFractionMin / 100.0) * RdTlp.size())
        { // DROP if minimal mapped segment do not exceed % of read's length
            alignPassFilters = false;
        }
        if(alignPassFilters)
        {
        // Filling
        //
        // -------------------------------
        //      S              E
        //       \/\/\/\/\/\/\/
        //
        //   TTATTA+++G++CTTA++G++CATTG   Rf
        //   XXA--ATTTGTTCTTATTGTTCAXXX   Rd
        //
        // RfTlp:                   *********
        //       T: A A G C A G C A G
        //      dT: 2 2 0 0 2 0 0 0 2 0
        //          0 1 2 3 4 5 6 7 8 9
        //            ^
        // RdTlp:
        //       T: X X A A G C A G C A X X X
        //      dT: 0 0 0 0 3 2 2 2 2 0 0 0 0
        //          0 1 2 3 4 5 6 7 8 9 0 1 2
        //                ^
        //  Alignment: rf_start     = 0
        //             rd_leftmost  = 2
        //             rd_rightmost = 9 (alignmentLen = 9-2 = 7)
        //  ref_E = 0+7 = 7
        //  dR = (0-2), (3-0), (2-0), (2-2), (2-0), (2-0), (0-0)
        //         -2     3      2      0      2      2      0
        //
        // for(int i = 0; i < alignmentLen-1; i++): 0123456
        // dR = RdDt[rd_leftmost + i + 1] - RfDt[rf_start + i + 1];

            bestAlignment.ref_Idx = longestMatch.ref;
            bestAlignment.ref_S   = longestMatch.rf_start;
            bestAlignment.ref_E   = longestMatch.rf_start + alignmentLen;
            bestAlignment.mmcount = mm;

            if(bestAlignment.ref_E > RfTlp.size() ||
                alignmentLen + longestMatch.rd_leftmost > RdTlp.size())
            {
                if(Alignments.size() > alignment_probe)
                {
                    alignment_probe++;
                    bestAlignment.flags &= ~TA::AF_Valid;
                    goto RepickAlignment;
                }
                string err = "Mapping error:\n";
                err += "\nbestAlignment.ref_S="; err += std::to_string(bestAlignment.ref_S);
                err += "\nbestAlignment.ref_E="; err += std::to_string(bestAlignment.ref_E);
                err += "\nlongestMatch.rd_leftmost="; err += std::to_string(longestMatch.rd_leftmost);
                err += "\nlongestMatch.rd_rightmost="; err += std::to_string(longestMatch.rd_rightmost);
                err += "\nalignmentLen="; err += std::to_string(alignmentLen);
                err += "\nRfTlp.size()="; err += std::to_string(RfTlp.size());
                err += "\nRdTlp.size()="; err += std::to_string(RdTlp.size());
                err += "\n\n";
                throw runtime_error(err);
            }

            const vector<difft_t>& RfDt = (*(T.refTless))[longestMatch.ref].dT;
            const vector<difft_t>& RdDt = T.read.dT;

            int distESD = 0;
            for(int i = 0; i < alignmentLen - 1; i++)
            // changed from '< alignmentLen' to preserve
            // [start-\/\/\/\/-end] alignment notation, where
            // the dRs count is alignmentLen-1
            {
                difft_t edR = RdDt[longestMatch.rd_leftmost + i + 1] -
                    RfDt[longestMatch.rf_start + i + 1];

                if(edR != 0) distESD++;
                bestAlignment.dR.push_back(edR);
            }

            bestAlignment.setVaild();
            bestAlignment.ref_ESD = distESD;
            // Reporting
            // set EXACT MATCH flag
            if(mm == 0) bestAlignment.flags |= TA::AF_ExactMatch;
            else bestAlignment.flags &= ~TA::AF_ExactMatch;
            // set REF MATCH flag
            if(distESD == 0) bestAlignment.flags |= TA::AF_MatchReference;
            else bestAlignment.flags &= ~TA::AF_MatchReference;
            // set MUPLIMAPPER flag
            //if(seedMatchCounter == 0) bestAlignment.flags |= TA::AF_Multimapper;
            //else bestAlignment.flags &= ~TA::AF_Multimapper;
        }
    }
    return bestAlignment;
}


#ifdef DOUNITTEST

void Test2()
{
    const char* referenceFile = "/Users/jalgard/Yandex.Disk/jalgard/projects/H10mini/assembly.fa";
    const char* fastqFile     = "/Users/jalgard/Yandex.Disk/jalgard/projects/H10mini/h10.win";

    vector<vector<string> > refs;
    vector<TlessDNA> refTless;
    unordered_map<string, vector<Rindex> > refIndex;

    ProcessReferenceFile(referenceFile, refs, refTless);
    Build_Fasta_Index(refTless, refIndex);

    int numReadsAligned = 0;
    int numReadsSeen    = 0;

    std::ifstream ifile(fastqFile); std::string line1, line2, line3, line4;
	while(  std::getline(ifile, line1, '\n') &&
			std::getline(ifile, line2, '\n') &&
			std::getline(ifile, line3, '\n') &&
			std::getline(ifile, line4, '\n'))
	{
		if(line1[0] == '@' && line2.size() > 1 && line2.size() == line4.size())
		{
			cout << "Aligning read " << line1 << "\n";
            ReadAlignmentTask RAT;
            RAT.refIndex = &refIndex;
            RAT.refTless = &refTless;
            RAT.read     = MakeTless(line2);
            auto A = TAlignRead(RAT);
            if(A.isValid()) cout << "Valid(true)\n";
            else cout << "Valid(false)\n";
            cout << A.ref_S << "\n";

            RAT.read     = MakeTless(ReverseAndComplementDNA(line2));
            auto B = TAlignRead(RAT);
            if(B.isValid()) cout << "Valid(true)\n";
            else cout << "Valid(false)\n";
            cout << B.ref_S << "\n";
		}
	}
}

#endif

#endif
