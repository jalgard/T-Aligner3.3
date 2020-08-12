// by J@:)
#ifndef TLESS_H
#define TLESS_H

#include <vector>
#include <string>
#include <iostream>

using difft_t = int;

using namespace std;

/*
    Notation of Tless DNA string for this version:
    i-th elemnt of T/dT arrays encodes 0-based
    i-th A/G/C letter of DNA string and preceeding
    number of Ts.

    DNA:    TTATTGCTGA
    TLP:     2 2 0 1 0 0
        T:   A G C G A
*/

class TlessDNA
{
public:
    string T;
    vector<difft_t> dT;
    TlessDNA() = default;
    TlessDNA(const string& P): T("") { dT.reserve(P.size()); }
};

static TlessDNA MakeTless(const string& P)
{
    TlessDNA tlp(P);
    difft_t dTaccum = 0;
    for(size_t i = 0; i < P.size(); i++)
    {
    	while(i < P.size() && P[i] == 'T')
        {
            dTaccum++;
            i++;
        }
        tlp.dT.push_back(dTaccum);
        dTaccum = 0;
    	if(i < P.size())
        {
    	   tlp.T += P.substr(i,1);
        }
    }

    if(tlp.dT.size() == tlp.T.size())
    {
        tlp.dT.push_back(0);
    }

    return tlp;
}

static string PrintTless(const TlessDNA& Tl)
{
    string R = "";
    size_t i = 0;
    for(; i < Tl.T.size(); i++)
    {
        for(difft_t t = 0; t < Tl.dT[i]; t++)
        {
            R += "T";
        }
        R += Tl.T.substr(i,1);
    }
    if(i == Tl.T.size())
    {
        for(difft_t t = 0; t < Tl.dT[i]; t++)
        {
            R += "T";
        }
    }

    return R;
}

static string DumpTless(const TlessDNA& Tl)
{
    string R = "[T]";
    for(size_t i = 0; i < Tl.T.size(); i++)
    {
        R += ("\t" + Tl.T.substr(i,1));
    }
    R += "\n[dT]";
    for(size_t i = 0; i < Tl.dT.size(); i++)
    {
        R += ("\t" + std::to_string(Tl.dT[i]));
    }
    return R;
}

// translate from T-less DNA coordinate to sequence coordinate
static int Translate(const TlessDNA& Tl, int pos)
{
    int seqp = 0;
    int p    = 0;

    // account for 5' Ts
    seqp += Tl.dT[0];

    while(p != pos)
    {
        seqp += 1; seqp += Tl.dT[p+1];
        p++;
    }

    return seqp;
}

// translate from sequence coordinate to T-less DNA coordinate
// returns the position of first non-T symbol if
// pos points to 'T'
static int Translate(const string& Dn, int pos)
{
    int tlp = -1;
    int p   = 0;

    while(p != pos + 1)
    {
        if(Dn[p] != 'T' &&
           Dn[p] != 't') tlp++;
        p++;
    }

    return tlp;
}

/*
    TESTS FOR TLESS.H
*/

#ifdef DOUNITTEST
void Test1()
{
    string DNA = "AATGTTCGTTTATAGAT";
    TlessDNA TLP = MakeTless(DNA);
    cout << "Doing TlessDNA tests:\n";

    cout << "DNA:\n" << DNA << "\n";
    cout << "PrintTless test:\n" << PrintTless(TLP) << "\n";
    cout << "DumpTless test:\n" << DumpTless(TLP) << "\nAll done!\n\n";
}
#endif
#endif
