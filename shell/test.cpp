#include "taldriver.h"
#include <fstream>



int main(int argc, char** argv)
{
    GffReader g;
    g.Load(argv[1]);
    auto item = g.Get(atoi(argv[3]));
    cout << item.chr << "; " << item.s << ", " << item.e << "; ";
    for(auto& x : item.tags) cout << " " << x.first << "=" << x.second << ", ";
    cout << "\n";

    auto f = FastaReader(argv[2]);

    cout << g.GetSeq(atoi(argv[3]), f[0][1]) << "\n\n";
    return 0;
}
