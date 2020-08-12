#include "gmatcher.h"
#include "managedna.h"


static map<string, float> InitScoringSystem()
{
    map<string, float> scoring_matrix;
    // [GRNA : MRNA] = SCORE
    //  TCTTTAAGTGTAACAGAAAATTAGG mRNA
    //  +++++++::::++++++++++++#:
    //  TCTTTAAACACAACAGAAAATTACA gRNA
    int score_em = 1.0;
    int score_gu = 0.0;
    int score_mm = -2.0;

    scoring_matrix["AA"] = score_em;
    scoring_matrix["TT"] = score_em;
    scoring_matrix["GG"] = score_em;
    scoring_matrix["CC"] = score_em;

    scoring_matrix["AG"] = score_gu;
    scoring_matrix["CT"] = score_gu;

    scoring_matrix["CA"] = score_mm;
    scoring_matrix["CG"] = score_mm;

    scoring_matrix["AC"] = score_mm;
    scoring_matrix["AT"] = score_mm;

    scoring_matrix["GT"] = score_mm;
    scoring_matrix["GA"] = score_mm;
    scoring_matrix["GC"] = score_mm;

    scoring_matrix["TA"] = score_mm;
    scoring_matrix["TC"] = score_mm;
    scoring_matrix["TG"] = score_mm;

    return scoring_matrix;
}

tuple<int, float, int> Align(int pGrna, string& gRNA, int pMrna, string& mRNA,
    float& init_score, float& min_score, map<string, float>& scoring_system)
{
    float score = 0; int run = 1; int last = pMrna;
    while(pGrna < gRNA.size() && pMrna < mRNA.size())
    {
        float bonus = scoring_system[gRNA.substr(pGrna,1) + mRNA.substr(pMrna,1)];
        if(bonus >= 0)
        {
            last = pMrna;
        }
        score += bonus;
        if(float(init_score + score) / run < min_score)
        {
            break;
        }
        else
        {
            pGrna++; pMrna++; run++;
        }
    }
    return make_tuple(last, score, run);
}


vector<tuple<string, string, int, int, int, int> > FindLongestGrna(string& gRNA, string& mRNA,
    float& init_score, float& min_score, map<string, float>& scoring_system,
    int& min_grna_length, int& mm_max, int& gu_max)
{
    vector<tuple<string, string, int, int, int, int> > results;
    int pMrna = 0;
    while(pMrna < mRNA.size())
    {
        int pGrna = 0;
        while(pGrna < gRNA.size())
        {
            auto align = Align(pGrna, gRNA, pMrna, mRNA, init_score, min_score, scoring_system);

            int aln_length = get<0>(align) - pMrna;
            if(aln_length > min_grna_length)
            {
                string MrnaSeq = mRNA.substr(pMrna, aln_length);
                string GrnaSeq = gRNA.substr(pGrna, aln_length);
                int GU = 0; int MM = 0;
                for(int p = 0; p < aln_length; p++)
                {
                    if(MrnaSeq[p] != GrnaSeq[p])
                    {
                        if( (MrnaSeq[p] == 'G' && GrnaSeq[p] == 'A') ||
                            (MrnaSeq[p] == 'T' && GrnaSeq[p] == 'C') )
                        {
                            GU++;
                        }
                        else
                        {
                            MM++;
                        }
                    }
                }
                if(MM <= mm_max && GU <= gu_max)
                {
                    results.push_back(make_tuple(GrnaSeq, MrnaSeq, pGrna, pMrna, MM, GU));
                    pMrna += aln_length;
                }
            }
            pGrna++;
        }
        pMrna++;
    }
    return results;
}

int main(int argc, char** argv)
{
    auto grna_fasta_data = FastaReader(argv[1], false);
    auto mrna_fasta_data = FastaReader(argv[2], false);

    int mm_max = 10;
    int gu_max = 15;
    float init_score = 4.0;
    float min_score  = 0.58;
    int min_grna_length = 20;

    map<string, float> scoring_system = InitScoringSystem();

    cerr << "Staring process for potential gRNA source:\n" << argv[1];
    cerr << "\nand mRNA dataset from:\n" << argv[2] << "\n";
    cerr << "Output file is set to: STDOUT\n\n";

    for(auto& mRNA : mrna_fasta_data)
    {
        cerr << "Inspecting " << mRNA[0] << " sequence of " << mRNA[1].size() << " bp.\n";
        for(auto& gRNA : grna_fasta_data)
        {
            cerr << "...... processing " << gRNA[0] << " of " << gRNA[1].size() << " bp.";
            auto results = FindLongestGrna(gRNA[1], mRNA[1], init_score, min_score,
                scoring_system, min_grna_length, mm_max, gu_max);
            for(auto& result: results)
            {
                cout << mRNA[0] << "\t" << get<3>(result) << "\t";
                cout << gRNA[0] << "\t" << get<2>(result) << "\t";
                cout << get<1>(result) << "\t" << get<0>(result) << "\t";
                cout << get<5>(result) << "\t" << get<4>(result) << "\n";
            }
            if(results.size() == 0)
            {
                cerr <<  "..... No hits found.\n";
            }
            else
            {
                cerr << "..... " << results.size() << " hits spotted.\n";
            }
        }
    }
    cerr << "\nAll done!\n\n";
}
