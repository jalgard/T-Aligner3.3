#include "snp.h"
#include "geneticcode.h"




int main(int argc, char** argv)
{
    auto snps = LoadSnpSet(argv[1]);
    GffReader gffReader;
    gffReader.Load(argv[2]);
    auto rps12_gene = gffReader.Get("RPS12");
    CastCoordinates(snps, rps12_gene.s);

    for(int i = 0; i < snps.size(); i++)
    {
        cout << snps[i].pos << "\t" << snps[i].ref << "\t" << snps[i].alt << "\n";
    }

    auto geneticCode = setGeneticCode("T-Aligner-default");
    auto codonNs     = setCodonNs(geneticCode);
    cout << codonNs["CGG"] << "\n";
    cout << codonNs["TGG"] << "\n";
    cout << codonNs["AAA"] << "\n";
    cout << codonNs["TTT"] << "\n";
    cout << codonNs["GTG"] << "\n";
    cout << codonNs["CGC"] << "\n";
    cout << codonNs["AGG"] << "\n";
    cout << codonNs["TAC"] << "\n";
    cout << codonNs["CAT"] << "\n";
    cout << codonNs["GGA"] << "\n";
    return 0;
}
