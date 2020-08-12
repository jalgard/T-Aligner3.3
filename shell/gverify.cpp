#include "gmatcher.h"
#include "managedna.h"

bool FileExists(const char *fileName)
{
    ifstream infile(fileName);
    return infile.good();
}

void CreateRnaLibraryKmerHash(const char* rnalib_fastq,
    unordered_map<string, int>& kmer_lib, int kmer_size = 16)
{
    cerr << "Start building kmer index for " << rnalib_fastq << "\n";
    cerr << "with kmer size " << kmer_size << "\n";
    ifstream ifile(rnalib_fastq);
	string line1, line2, line3, line4;
    while(  std::getline(ifile, line1, '\n') &&
			std::getline(ifile, line2, '\n') &&
			std::getline(ifile, line3, '\n') &&
			std::getline(ifile, line4, '\n'))
	{
		if(line1[0] == '@' && line2.size() > kmer_size && line2.size() == line4.size())
		{
            auto rc = ReverseAndComplementDNA(line2);
            for(int p = 0; p < line2.size() - kmer_size; p++)
            {
                kmer_lib[line2.substr(p, kmer_size)]++;
                kmer_lib[rc.substr(p, kmer_size)]++;
            }
        }
    }
    ifile.close();
    cerr << "Kmer indexing finished successfully!\n";
    cerr << "Totally " << kmer_lib.size() << " kmers extracted\n";

    cerr << "Writing to file " << rnalib_fastq << "." << kmer_size << "hash\n";
    ofstream ofile((string(rnalib_fastq) + std::to_string(kmer_size) + "hash").c_str());

    for(auto& key : kmer_lib)
    {
        ofile << key.first << "\t" << key.second << "\n";
    }
    ofile.close();
}

void LoadKmerLib(const char* rnalib_fastq,
    unordered_map<string, int>& kmer_lib, int kmer_size = 16)
{
    cerr << "Loading kmer library file\n";
    ifstream ifile((string(rnalib_fastq) + std::to_string(kmer_size) + "hash").c_str());

    string line1;
    while(std::getline(ifile, line1, '\n'))
    {
        auto toks = split(line1, '\t');
        kmer_lib[toks[0]] = stoi(toks[1]);
    }
    cerr << "Done!\n";
}

void CheckGrnaProtocol(const char* grna_predictions_file,
    const char* grna_verified_output_file,
    unordered_map<string, int>& kmer_lib,
    int which_row_to_check = 5,
    int kmer_size = 16, int minimal_support = 2)
{
    ifstream ifile(grna_predictions_file);
    ofstream ofile(grna_verified_output_file);
    string line1;
    cerr << "Start checking gRNA patterns from " << grna_predictions_file << "\n";
    cerr << "Kmer size is " << kmer_size << ", minimal support value is " << minimal_support << "\n";
    cerr << "Output is written to " << grna_verified_output_file << "\n";
    int total_seen = 0;
    int total_good = 0;
    while( std::getline(ifile, line1, '\n') )
    {
        auto toks = split(line1, '\t');
        total_seen++;
        int max_support_value = 0;
        int kmers_with_support = 0;
        for(size_t p = 1; p < toks[which_row_to_check].size() - kmer_size - 1; p++)
        {
            string pat = toks[which_row_to_check].substr(p, kmer_size);
            if(kmer_lib.find(pat) != kmer_lib.end())
            {
                int supp = kmer_lib[pat];
                if(supp >= minimal_support)
                {
                    if(supp > max_support_value)
                    {
                        max_support_value = supp;
                    }
                    kmers_with_support++;
                }
            }
        }
        if(kmers_with_support >= toks[which_row_to_check].size() - kmer_size - 3)
        {
            ofile << line1 << "\t" << max_support_value << "\n";
            total_good++;
        }
    }
    cerr << "Verification finished successfully!\n";
    cerr << "From total " << total_seen << " gRNAs seen " << total_good << " were found to be good!\n";
}

int main(int argc, char** argv)
{
    unordered_map<string, int> kmer_lib;

    if(argc < 4)
    {
        cerr << "Running in kmer library generation mode\n";
        // params: kmer_lib_fastq,  kmer_size

        CreateRnaLibraryKmerHash(argv[1], kmer_lib, atoi(argv[2]));
    }
    else
    {
        cerr << "Running in verification mode\n";
        // params: kmer_lib_fastq,  kmer_size,  grna_input_file, output_file,  min_supp,  row_to_check

        LoadKmerLib(argv[1], kmer_lib, atoi(argv[2]));

        CheckGrnaProtocol(argv[3], argv[4], kmer_lib, atoi(argv[6]), atoi(argv[2]), atoi(argv[5]));
    }

    return 0;
}
