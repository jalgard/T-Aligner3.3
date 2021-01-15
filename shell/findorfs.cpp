
#include "orffinder.h"
#include <fstream>


int FindOrfs()
{

    vector<vector<string> > refsBank;
    vector<TlessDNA> refTless;
    unordered_map<string, vector<Rindex> > refIndex;

    cerr << "Starting protocol: ORF finding\n";

    string inFastaFile = TAlignerOptions::Options().getOption("Input|Reference|Filename");
    ProcessReferenceFile(inFastaFile.c_str(), refsBank, refTless);
    Build_Fasta_Index(refTless, refIndex);

    ReadAlignmentTask RAT;
    RAT.refIndex = &refIndex;
    RAT.refTless = &refTless;

    int nThreads = thread::hardware_concurrency() - 1;
	int userThreads = TAlignerOptions::Options().getOption("T-Aligner|Threads");
    if(userThreads > 0) nThreads = userThreads;

    cerr << "Alignment stage: using " << nThreads << " threads\n";
    string fastqInputFile = TAlignerOptions::Options().getOption("Input|Library1|Filename");
    auto AlignmentsBank = SingleEndAlignFastqMT(fastqInputFile.c_str(), RAT);

    cerr << "Alignment stage: counting reads matching reference...\n";

    int matchingRef = 0;
    for(int i = 0; i < AlignmentsBank.size(); i++)
    {
        if(AlignmentsBank[i].doesMatchRef()) matchingRef++;
    }

    cerr << "Alignment stage:\tAligned:" << AlignmentsBank.size() << "\tMatching ref:" << matchingRef << "\n";
    cerr << "@TAL:\tAligned:\t" << AlignmentsBank.size() << "\tMatching ref:\t" << matchingRef << "\n";
    cerr << "Alignment stage: filtering duplicates...\n";
    auto dedupAlignBank = FilterMappingDuplicatesReads(AlignmentsBank);
    cerr << "@TAL:\tAlignments after dedup:\t" << dedupAlignBank.size() << "\n";


    cerr << "Assembly stage: building overlap graph... ";
    auto OGr = BuildOverlapGraph(dedupAlignBank);
    cerr << "final graph has " << OGr[-1][-1] << " edges.\n";

    string prefMode = TAlignerOptions::Options().getOption("ORFinder|PreferredMode");
    cerr << "Assembly stage: searching ORFs with " << prefMode << " algorithm.\n";

    vector<vector<int> > ORFs;
    for(int i = 0; i < dedupAlignBank.size(); i++)
    {
        TraceORF(i, OGr, dedupAlignBank, ORFs);
    }

    cerr << "@TAL:\tTotal traced attempts completed:\t" << ORFs.size() << "\n";

    cerr << "Assembly stage: processing ORFs\n";

//    sort(ORFs.begin(), ORFs.end(), [](const vector<int> & a, const vector<int> & b)
//        {  return a.size() > b.size(); });

    vector<tuple<string, string, TAlignment> > orfs_vector;

    map<string, int> orfs_found;

    cerr << "Decoding and translating...\n";
    int orf_length_min_aa = TAlignerOptions::Options().getOption("ORFinder|MinOrfLength|aa");

    string genetic_code_name = TAlignerOptions::Options().getOption("ORFinder|GeneticCode|Table");

    auto genetic_code = setGeneticCode(genetic_code_name);


    for(int i = 0; i < ORFs.size(); i++)
    {
        auto result = DecodeORF(ORFs[i], dedupAlignBank, refTless, genetic_code);
        auto cds_seq  = get<0>(result);
        auto pep_seq   = get<1>(result);
        auto alignment = get<2>(result);
        bool complete  = get<3>(result);
        if(orfs_found.count(cds_seq) == 0 &&
           orfs_found.count(pep_seq)  == 0 &&
           complete)
        {
            orfs_vector.push_back({cds_seq, pep_seq, alignment});
            orfs_found[cds_seq] = 1;
            orfs_found[pep_seq]  = 1;
        }
    }

    sort(orfs_vector.begin(), orfs_vector.end(),
        [](const tuple<string, string, TAlignment>& a,
            const tuple<string, string, TAlignment>& b)
        {  return get<1>(a).size() > get<1>(b).size(); });


    cerr << "Writing output files...\n";
    string output_prefix = TAlignerOptions::Options().getOption("Output|Prefix");
    string output_cds_fasta = TAlignerOptions::Options().getOption("Output|Files|OrfFinderCdsFasta");
    string output_pep_fasta = TAlignerOptions::Options().getOption("Output|Files|OrfFinderPepFasta");
    string output_mrna_fasta = TAlignerOptions::Options().getOption("Output|Files|OrfFinderMrnaFasta");
    string output_mrna_taf = TAlignerOptions::Options().getOption("Output|Files|OrfFinderMrnaTaf");

    ofstream ofile_cds_fasta((output_prefix + output_cds_fasta).c_str());
    ofstream ofile_pep_fasta((output_prefix + output_pep_fasta).c_str());
    ofstream ofile_mrna_fasta((output_prefix + output_mrna_fasta).c_str());
    ofstream ofile_mrna_taf((output_prefix + output_mrna_taf).c_str());

    int orf_name_suffix = 1;
    string reference_name_prefix = refsBank[0][0];

    string print_each_orf_alignment =
        TAlignerOptions::Options().getOption("ORFinder|OutputEachOrfAlignment");


    string orf_mode = TAlignerOptions::Options().getOption("ORFinder|PreferredMode");

    for(auto&& orf : orfs_vector)
    {
        if(get<1>(orf).size() < orf_length_min_aa) continue;

        //ofile_cds_fasta << ">" << reference_name_prefix << " assembled ORF " << orf_name_suffix << "\n";
        //ofile_pep_fasta << ">" << reference_name_prefix << " assembled ORF (translated, " << genetic_code_name  <<  ") " << orf_name_suffix << "\n";
        //ofile_mrna_fasta << ">" << reference_name_prefix << " assembled mRNA " << orf_name_suffix << "\n";

        ofile_cds_fasta << ">" << reference_name_prefix << "_" << orf_mode << "_ORF_" << orf_name_suffix << "\n";
        ofile_pep_fasta << ">" << reference_name_prefix << "_" << orf_mode << "_ORF_" << orf_name_suffix << "\n";
        ofile_mrna_fasta << ">" << reference_name_prefix << "_" << orf_mode << "_mRNA_" << orf_name_suffix << "\n";

        ofile_cds_fasta << get<0>(orf) << "\n";
        ofile_pep_fasta << get<1>(orf) << "\n";
        ofile_mrna_fasta << Alignment2Seq(get<2>(orf), refTless) << "\n";

        ofile_mrna_taf << PrintReadTaf(get<2>(orf), refTless) << "\n";

        if(print_each_orf_alignment == "true")
        {

            ofstream ofile_alignment((output_prefix + "orf_" +
                std::to_string(orf_name_suffix) + "_alignment.fasta").c_str());
            string orf_name = "ORF" + std::to_string(orf_name_suffix);
            ofile_alignment << PrintAlignedReadFasta(get<2>(orf), refTless, "ref", orf_name);
            ofile_alignment.close();
        }

        orf_name_suffix++;
    }

    ofile_cds_fasta.close();
    ofile_pep_fasta.close();
    ofile_mrna_fasta.close();
    ofile_mrna_taf.close();

    return 0;

}


int main(int argc, char** argv)
{
    TAlignerOptions::Options().StartLogging();
    TAlignerOptions::Options().Log(ParseFromCommandLine(argc, argv));
    int result = FindOrfs();
    return result;
}


/*
for i in $(ls ../WSD_REF/*.fasta);
do
name=$(echo $i | awk '{split($1, aa, ".fa"); split(aa[1], bb, "F/"); print bb[2];}');
echo Working with $name;
/home/jalgard/ta330/bin/findorfs --in_ref $i --in_lib ../RNA_DATASETS/WSD.paired.merged.strict.fq --orf_filter_esd 0 --out_prefix $name
done
*/
