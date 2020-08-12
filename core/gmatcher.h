#include <iostream>
#include <map>
#include <fstream>
#include <unordered_map>
#include <random>

#include "managedna.h"



using namespace std;

class AlignmentScoringScheme
{

public:

    AlignmentScoringScheme() { setDefaultScoring(); }

    void setDefaultScoring()
    {
        edit_bonus     = 4;
        gu_bonus       = 4;
        match_bonus    = 4;
        mismatch_bonus = -7; // -10
        mismatch_max   = 5;
        gu_max         = 10; //10
        score_min      = 60; //50
        // use dynamic scoring system
        per_step       = 2.5;
    }

    void setDefaultScoring_v1()
    {
        edit_bonus     = 2;
        gu_bonus       = 1;
        match_bonus    = 4;
        mismatch_bonus = -10;
        mismatch_max   = 2;
        gu_max         = 9;
        score_min      = 70;
    }

    int edit_bonus;
    int gu_bonus;
    int match_bonus;
    int mismatch_bonus;
    int mismatch_max;
    int gu_max;
    int score_min;
    float per_step;

};

class GrnaAlignment
{
public:

    GrnaAlignment() : aligned_rd(""), aligned_rf(""), align_action(""),
        rf_begin(0), rd_begin(0), mm(0), gu(0), score(0),
        state(0), pos_pointer(0) {}
    int score;
    int mm;
    int gu;
    int rf_begin;
    int rd_begin;
    int state;
    int pos_pointer;
    string aligned_rd;
    string aligned_rf;
    string align_action;


};

class GrnaMatcher
{
private:

    struct GrnaNode
    {
        GrnaNode() : N('N'), insertion_node(false) {}
        GrnaNode(char n, bool as = false) : N(n), insertion_node(as) {}
        vector<vector<int > > transitions;
        char N;
        bool insertion_node;

    };

    map<int, GrnaNode> reference_graph_E;
    map<int, int> ref_coordinate_caster;
    string t_less_ref;
    vector<int> tlr_t_holder;
    AlignmentScoringScheme scoring_scheme;

    bool GU(char G, char R)
    {
        if( (G == 'C' && R == 'T') || (G == 'A' && R == 'G') ) return true;
        return false;
    }

    void tless(const string& ref, string& tl_ref, vector<int>& t_holder, map<int, int>& rcc)
    {
        tl_ref = "";
        int p = 0;

        for(int i = 0; i < ref.length(); i++)
        {
            int t = 0;
            while(ref[i] == 'T') { i++; t++; }
            t_holder.push_back(t);
            tl_ref = tl_ref + ref[i];
            rcc[p] = i;
            p++;
        }
    }


public:

    int refSizeTless()
    {
        return t_less_ref.size();
    }

    int getScoreMin()
    {
        return scoring_scheme.score_min;
    }

    vector<int> getAlignmentPath(GrnaAlignment& al)
    {
        // TODO!!!
        vector<int> path;
        int t_less_pos = al.rf_begin+1;
        int t_holder = 0;
        for(int pp = 0; pp < al.align_action.size(); pp++)
        {
            int dT = tlr_t_holder[t_less_pos];
            if(al.align_action[pp] == 'M' || al.align_action[pp] == '*')
            {
                dT = t_holder - dT;
                t_holder = 0;
                t_less_pos++;
                path.push_back(dT);
            }
            else
            {
                t_holder++;
            }
        }
        return path;
    }

    void setScoringScheme(AlignmentScoringScheme user_scoring_scheme)
    {
        scoring_scheme = user_scoring_scheme;
    }

    void init(const string& ref)
    {

        tless(ref, t_less_ref, tlr_t_holder, ref_coordinate_caster);
        auto& tl = t_less_ref;
        int nc = 0;
        GrnaNode N0 = GrnaNode(tl[0], false);
        reference_graph_E[nc] = N0;

        for(int p = 1; p < tl.size(); p++)
        {
            char nuc_next = tl[p];

            GrnaNode NT = GrnaNode('T', true);
            GrnaNode NC = GrnaNode('C', true);

            reference_graph_E[nc].transitions.push_back({nc, nc+1}); // T
            reference_graph_E[nc+1] = NT;
            reference_graph_E[nc].transitions.push_back({nc, nc+2}); // C
            reference_graph_E[nc+2] = NC;

            GrnaNode NN = GrnaNode(nuc_next, false);

            reference_graph_E[nc].transitions.push_back({nc, nc+3}); // Nn
            reference_graph_E[nc+3] = NN;

            reference_graph_E[nc+1].transitions.push_back({nc+1, nc+1}); // T->T
            reference_graph_E[nc+1].transitions.push_back({nc+1, nc+2}); // T->C
            reference_graph_E[nc+2].transitions.push_back({nc+2, nc+2}); // C->C
            reference_graph_E[nc+2].transitions.push_back({nc+2, nc+1}); // C->T
            reference_graph_E[nc+1].transitions.push_back({nc+1, nc+3}); // T->Nn
            reference_graph_E[nc+2].transitions.push_back({nc+2, nc+3}); // C->Nn

            nc = nc + 3;
        }
    }

    vector<GrnaAlignment> matchNext(const string& grna, GrnaAlignment align)
    {

        vector<GrnaAlignment> results;
        if(align.pos_pointer >= grna.size()-1)
        {
            return vector<GrnaAlignment>({align});
        }
        char next_n   = grna[align.pos_pointer];
        int  state    = align.state;
        auto cur_node = reference_graph_E[state];

        for(auto& transition : cur_node.transitions)
        {

            auto move = reference_graph_E[transition[1]];
            GrnaAlignment next_align = align;
            next_align.pos_pointer++;
            next_align.state = transition[1];

            // TODO
            //bool ref_contains_T = ?;

            if(move.N == next_n)  // match transition
            {
                if(move.insertion_node)
                {
                    if(move.N == 'T')
                    {
                        next_align.score += scoring_scheme.edit_bonus;
                    }
                    else
                    {
                        next_align.score += scoring_scheme.edit_bonus;
                        next_align.gu++;
                    }

                    next_align.aligned_rd   += next_n;
                    next_align.aligned_rf   += "T";
                    next_align.align_action += "+";
                }
                else
                {
                    next_align.score += scoring_scheme.match_bonus;
                    next_align.aligned_rd   += next_n;
                    next_align.aligned_rf   += move.N;
                    next_align.align_action += "M";
                }
            }
            // GU transition
            else if(!move.insertion_node && GU(next_n, move.N))
            {
                next_align.score += scoring_scheme.gu_bonus;
                next_align.gu++;
                next_align.aligned_rd += next_n;
                next_align.aligned_rf += move.N;
                next_align.align_action += "M";
            }
            // MM transition
            else if(!move.insertion_node) // move.N can't be T by definition 
            {
                next_align.score += scoring_scheme.mismatch_bonus;
                next_align.mm++;
                next_align.aligned_rd += next_n;
                next_align.aligned_rf += move.N;
                next_align.align_action += "*";
            }
            // break other invalid alignments
            else
            {
                next_align.mm = 1 + scoring_scheme.mismatch_max;
            }
            if(next_align.mm <= scoring_scheme.mismatch_max && next_align.gu <= scoring_scheme.gu_max)
            {
                auto r = matchNext(grna, next_align);
                if(r.size() > 0)
                    results.insert(results.end(), r.begin(), r.end());
            }
            else
            {
                    if(next_align.score >= scoring_scheme.score_min)
                    {
                        results.push_back(next_align);
                    }
            }
        }

        vector<GrnaAlignment> valid_alignments;
        for(auto& result : results)
        {
            // heuristic filtration rules
            // 1. should start with at least 2 matches
            // 2. should have at least 4 nuc-long anchor
            // 3. any two mismathces should not go together
            // 4. mismatch can't be between two editing sites
            // DYNAMIC scoring system
            if(result.score >= scoring_scheme.score_min &&
               float(result.score) >= scoring_scheme.per_step *
               float(result.align_action.size()) &&
                //result.align_action.substr(0,2) == "MM" &&
                //result.align_action.substr(result.align_action.size()-4) == "MMMM" &&
                result.align_action.find("**") == string::npos
                //result.align_action.find("+*+") == string::npos
            )
            {
                valid_alignments.push_back(result);
            }
        }

        return valid_alignments;
    }


    vector<GrnaAlignment> MatchGrna(const string& grna, int refs = 0, int minlen = 25, int grnas = 0)
    {
        if(grna.size() < minlen) return vector<GrnaAlignment>({});

        vector<GrnaAlignment> finalResults;
        GrnaAlignment initialAlignment = GrnaAlignment();
        initialAlignment.state = 3 * refs;
        vector<GrnaAlignment> subResults = matchNext(grna, initialAlignment);
        for(int i = 0; i < subResults.size(); i++)
        {
            if(subResults[i].aligned_rd.size() >= minlen)
            {
                subResults[i].rd_begin = grnas;
                subResults[i].rf_begin = refs; //ref_coordinate_caster[refs];
                finalResults.push_back(subResults[i]);
            }
        }
        return finalResults;
    }

    vector<GrnaAlignment> MatchGrnaOverSequence(const string& grna, int minlen = 25, int grnas = 0)
    {
        vector<GrnaAlignment> finalResults;
        for(int pp = 0; pp < refSizeTless(); pp++)
        {
            auto subResults = MatchGrna(grna, pp, minlen, grnas);
            for(int i = 0; i < subResults.size(); i++)
            {
                finalResults.push_back(subResults[i]);
            }
        }
        return finalResults;
    }

};


void FindGrnaCandidatesWrapper(
    const char* i_crypto_file,
    const char* i_guiding_file,
    const char* o_duplex_file,
    const int gmer = 25,
    const int gwin = 19)
{
    auto crypto_fasta = FastaReader(i_crypto_file, false);
    auto guiding_fasta = FastaReader(i_guiding_file, false);
    ofstream o_duplex_stream(o_duplex_file);
    unordered_map<string, int> filtered_results;

    for(auto& cryptogene : crypto_fasta)
    {
        auto& crypto_name = cryptogene[0];
        auto& crypto_seq  = cryptogene[1];
        GrnaMatcher matcher = GrnaMatcher();
        matcher.init(crypto_seq);

        for(auto& guiderna : guiding_fasta)
        {
            auto grna_name = guiderna[0];
            auto grna_seq  = guiderna[1];

            for(int p = 0; p < static_cast<int>(grna_seq.size()) - gmer; p++)
            {
                string grna = grna_seq.substr(p, gmer);
                auto results = matcher.MatchGrnaOverSequence(grna, gwin, p);

                for(auto& duplex : results)
                {
                    if(duplex.score > gmer * 2)
                    {
                        auto align_dR = matcher.getAlignmentPath(duplex);
                        string align_str = "";
                        for(int pp = 0; pp < align_dR.size(); pp++)
                        {
                            align_str += ",";
                            align_str += std::to_string(align_dR[pp]);
                        }

                        string data_row = crypto_name + "\t" + std::to_string(duplex.rf_begin) +
                            "\t" + grna_name + "\t" + std::to_string(p) + "\t" + std::to_string(duplex.score) +
                            "\t" + duplex.aligned_rf + "\t" + duplex.aligned_rd + "\t" + duplex.align_action +
                            "\t" + align_str.substr(1);
                        filtered_results[data_row]++;
                        /*
                            writes output as
                            mRNA name | mRNA start | gRNA source name | gRNA source start | score | mRNA seq | gRNA seq | alignment action
                        */
                    }
                }
            }
        }

    }

    /*
        sort results by mRNA name, then by mRNA start
    */
    vector<tuple<string, string, int> > sorted_results;
    for(auto& r : filtered_results)
    {
        sorted_results.push_back({r.first, split(r.first, '\t')[0], std::stoi(split(r.first, '\t')[1])});
    }

    sort(sorted_results.begin(), sorted_results.end(),
        [](const tuple<string, string, int> & a, const tuple<string, string, int> & b) -> bool
        {
            return get<2>(a) > get<2>(b);
        });
    sort(sorted_results.begin(), sorted_results.end(),
        [](const tuple<string, string, int> & a, const tuple<string, string, int> & b) -> bool
        {
            return get<1>(a) > get<1>(b);
        });

    for(auto& r : sorted_results)
    {
        o_duplex_stream << get<0>(r) << "\n";
    }
    o_duplex_stream.close();

}
