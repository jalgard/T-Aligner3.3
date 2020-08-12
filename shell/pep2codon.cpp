#include "taldriver.h"
#include <fstream>


void TransferProteinAlignmentToMrna(int argc, char** argv, bool split_space = false)
{
    // read aligned protein fasta file
    prot_alignment = FastaReader(argv[1], split_space);
    // read correspoding CDS fasta files
    cds_dataset = FastaReader(argv[2], split_space);
    // write output
    cds_output = ofstream(argv[3]);

    // check if aligned protein names are also in CDS dataset
    for(auto& prot : prot_alignment)
    {
        bool name_found = false;
        bool length_check_passed = false;
        int prot_len = 0;
        for(int i = 0; i < prot[1].size(); i++)
        {
            if(prot[1][i] != '-')
            {
                prot_len++;
            }
        }
        // loop over CDSs
        for(auto& cds : cds_dataset)
        {
            if(cds[0] == prot[0])
            {
                name_found = true;
                int cds_len = 0;
                for(int i = 0; i < cds[1].size(); i++)
                {
                    if(cds[1][i] != '-')
                    {
                        cds_len++;
                    }
                }
                if(cds_len / 3 == prot_len)
                {
                    length_check_passed = true;
                }
                break;
            }
        }
        // report consistency
        if(!name_found)
        {
            cerr << "ERROR (In dataset)! No matching CDS sequence for " << prot[0] << " found in dataset " << argv[2] << "!\n\n";
            return 1;
        }
        if(!length_check_passed)
        {
            cerr << "ERROR (In dataset)! Length of protein sequence " << prot[0] <<
                    " and length of the corresponding CDS sequence does not match!\n\n";
            return 2;
        }
    }

    

    cds_output.close();
}

int main(int argc, char** argv)
{


}
