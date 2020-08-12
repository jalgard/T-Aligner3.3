
#ifndef TA_OPTIONS_H
#define TA_OPTIONS_H

#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include <ctime>

using namespace std;

struct TAOptionValue
{

    string ta_value;

    template<typename T> operator T() const
    {
       stringstream ss(ta_value);
       T convertedValue;
       if (ss >> convertedValue) return convertedValue;
       else throw runtime_error("Can't convert option value \"" + ta_value + "\"");
    }

};

class TAlignerOptions final
{
// functions
public:

    static TAlignerOptions& Options() {
        static TAlignerOptions TAinstance;
        return TAinstance;
    }

    void setOption(string option, string value)
    {
        OptionsHolder[option] = value;
    }

    void setOption(string option, double value)
    {
        setOption(option, to_string(value));
    }

    void setOption(string option, int value)
    {
        setOption(option, to_string(value));
    }

    TAOptionValue getOption(string option)
    {
        return GetValue(option);
    }

    TAlignerOptions(TAlignerOptions const&) = delete;
    //void operator=(TAlignerOptions const&)  = delete;
    TAlignerOptions(TAlignerOptions&&) = delete;                  // Move construct
    TAlignerOptions& operator=(TAlignerOptions const&) = delete;  // Copy assign
    TAlignerOptions& operator=(TAlignerOptions &&) = delete;      // Move assign

    inline string TimestampFull()
    {
        time_t rawtime; struct tm * timeinfo; char buffer [80]; time (&rawtime);
        timeinfo = localtime (&rawtime);
        strftime (buffer,80,"[ %A   %d %b %G   %H:%M:%S ] ",timeinfo);
        return string(buffer);
    }
    inline string Timestamp()
    {
        time_t rawtime; struct tm * timeinfo; char buffer [80]; time (&rawtime);
        timeinfo = localtime (&rawtime);
        strftime (buffer,80,"[ %H:%M:%S ] ",timeinfo);
        return string(buffer);
    }
    void StartLogging()
    {
        Logfile.open("taligner.log");
        Logfile << TimestampFull() << "Starting run.\n";
    }
    void Log(const string& logdata)
    {
        Logfile << Timestamp() << logdata << "\n";
    }
    void StopLogging()
    {
        Logfile.close();
    }
// functions
private:

    ofstream Logfile;

    map<string, string> OptionsHolder

    /*
        Default options initialization block

        values of boolean type should be set as 0/1
    */
    {
        { "T-Aligner|Version", "4.0.SA" },
        { "T-Aligner|Threads", "-1" },
        { "T-Aligner|Reads/job", "Hi" },
        { "T-Aligner|MemoryUsePreset|Hi", "400000" },
        { "T-Aligner|MemoryUsePreset|Mid", "200000" },
        { "T-Aligner|MemoryUsePreset|Low", "5000" },

        // legacy options 'Workflow'
        // not shell python scripts and plotters are available
        { "Workflow|MapReads", "1" },
        { "Workflow|FindORFs", "1" },
        { "Workflow|Draw", "1" },

        { "Aligner|Seed", "10" },
        { "Aligner|Mismatch", "2" },
        { "Aligner|MinMappedSegment", "10" },
        { "Aligner|MinMapped%", "70" },

        { "ORFinder|OverlapMin", "10" },
        { "ORFinder|OverlapGraph|FilterEsd", "5" },
        { "ORFinder|PreferredMode", "overlap" },   // can be 'overlap',
                                                   // 'all', 'coverage',
                                                   // 'editing' or 'extension'
        { "ORFinder|Tracing|SearchDepth", "1" },  // how many best positions to trace
        { "ORFinder|MinOrfLength|aa", "50" },
        { "ORFinder|GeneticCode|Table", "T-Aligner-default" },
        { "ORFinder|OutputEachOrfAlignment", "false" },

        //{ "ORFinder|Homologs|FastaInput", "-" }, // should be set(!) if orf tracing mode is 'homology'

        // suggested names for output files by default
        { "Output|Prefix", "TA_run_" },
        { "Output|Files|MappedReadsFastq", "mapped_reads.fastq" },
        { "Output|Files|MappedReadsTaf", "mapped_reads.taf" },
        { "Output|Files|ReadAlignmentsMatrix", "alignments.tam" },
        { "Output|Files|OrfFinderCdsFasta", "assembled_CDS.fasta" },
        { "Output|Files|OrfFinderPepFasta", "assembled_pep.fasta" },
        { "Output|Files|OrfFinderMrnaFasta", "assembled_mrna.fasta" },
        { "Output|Files|OrfFinderMrnaTaf", "assembled_mrna.taf" },



        { "Output|Files|Pdf", "1" },


        { "Input|Libraries|Count", "0" },

    };

    //TAlignerOptions* TAinstance;
    TAlignerOptions() {};

    TAOptionValue GetValue(string key)
    {
        if(OptionsHolder.count(key))
        {
            string value = OptionsHolder[key];
            return { value };
        }
        else throw runtime_error("Invalid option \"" + key + "\"");
    }

};

string ParseFromCommandLine(int argc, char** argv)
{
    auto& TOP = TAlignerOptions::Options();
    // iterate through options
    for(int oId = 0; oId < argc; oId++)
    {
        string ckey = string(argv[oId]);

        // KEY: library with reads
        if(ckey == "--in_lib")
        {
            int clib = TOP.getOption("Input|Libraries|Count");
            TOP.setOption("Input|Libraries|Count", ++clib);
            TOP.setOption("Input|Library" + to_string(clib) + "|Filename", string(argv[++oId]));
            TOP.setOption("Input|Library" + to_string(clib) + "|Type", "SE");
            if(oId+1 < argc)
            {
                string ltype = string(argv[oId+1]);
                if(ltype == "PE")
                {
                    TOP.setOption("Input|Library" + to_string(clib) + "|Type", ltype);
                }
            }
        }   // end library fastq

        // KEY: reference fasta file
        else if(ckey == "--in_ref")
        {
            TOP.setOption("Input|Reference|Filename", string(argv[++oId]));
        }   // end library fasta

        // KEY: prefix for output files
        else if(ckey == "--out_prefix")
        {
            TOP.setOption("Output|Prefix", string(argv[++oId]));
        }

        // KEY: set genetic code table name
        else if(ckey == "--orf_codon_table")
        {
            TOP.setOption("ORFinder|GeneticCode|Table", string(argv[++oId]));
        }

        // KEY: change preferred mode DEFAULT: overlap
        else if(ckey == "--orf_tracing_mode")
        {
            string mode = string(argv[++oId]);
            if(mode == "overlap" || mode == "editing" || mode == "extension"
                || mode == "all" || mode == "coverage")
            {
                TOP.setOption("ORFinder|PreferredMode", string(mode));
            }
        }

        // KEY: set treshold for ESD when building overlap graph
        else if(ckey == "--orf_search_depth")
        {
            TOP.setOption("ORFinder|Tracing|SearchDepth", string(argv[++oId]));
        }

        // KEY: set treshold for ESD when building overlap graph
        else if(ckey == "--orf_filter_esd")
        {
            TOP.setOption("ORFinder|OverlapGraph|FilterEsd", string(argv[++oId]));
        }

        // KEY: set treshold for overlap when building overlap graph and ORF tracing
        else if(ckey == "--orf_min_overlap")
        {
            TOP.setOption("ORFinder|OverlapMin", string(argv[++oId]));
        }

        // KEY: determine if orf alignments will be printed in .fasta
        // files 'true'/'false'    DEFAULT: false
        else if(ckey == "--orf_output_orf_alignments")
        {
            TOP.setOption("ORFinder|OutputEachOrfAlignment", string(argv[++oId]));
        }

        else if(ckey == "--orf_min_orf_aa")
        {
            TOP.setOption("ORFinder|MinOrfLength|aa", string(argv[++oId]));
        }

        else if(ckey == "--aln_mismatch_max")
        {
            TOP.setOption("Aligner|Mismatch", string(argv[++oId]));
        }

        else if(ckey == "--aln_min_segment")
        {
            TOP.setOption("Aligner|MinMappedSegment", string(argv[++oId]));
        }

        else if(ckey == "--aln_min_mapped_percent")
        {
            TOP.setOption("Aligner|MinMapped%", string(argv[++oId]));
        }


    }

    string input_command_line = "T-Aligner command line:";
    for(int i = 0; i < argc; i++)
    {
        input_command_line += " ";
        input_command_line += string(argv[i]);
    }
    return input_command_line;
}


#endif
