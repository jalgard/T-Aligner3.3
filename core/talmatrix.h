// by J@:)
#ifndef TALMAT_H
#define TALMAT_H
#include<functional>

vector<vector<TAlignment> > AlignmentSplitterGeneric(vector<TAlignment>& alignments, std::function<bool(TAlignment&)> filterFunc)
{
	vector<TAlignment> rAlignedReads, lAlignedReads;

	for(size_t i = 0; i < alignments.size(); i++)
	{
		if(filterFunc(alignments[i]))
		{
			rAlignedReads.push_back(alignments[i]);
		}
		else
		{
			lAlignedReads.push_back(alignments[i]);
		}
	}

	vector<vector<TAlignment> > Split;
	Split.push_back(rAlignedReads);
	Split.push_back(lAlignedReads);
	return Split;
}

vector<vector<TAlignment> > AlignmentSplitterByRevcom(vector<TAlignment>& alignments)
{
	std::function<bool (TAlignment&)> filter = [](TAlignment& tal)
	{
		return tal.isReversed();
	};
	return AlignmentSplitterGeneric(alignments, filter);
}

vector<vector<TAlignment> > AlignmentSplitterByMatchRef(vector<TAlignment>& alignments)
{
	std::function<bool (TAlignment&)> filter = [](TAlignment& tal)
	{
		return tal.doesMatchRef();
	};
	return AlignmentSplitterGeneric(alignments, filter);
}

vector<vector<TAlignment> > AlignmentSplitterByRefIdx(
	vector<TAlignment>& alignments, int refIdx)
{
	std::function<bool (TAlignment&)> filter = [=](TAlignment& tal)
	{
		return tal.ref_Idx == refIdx;
	};
	return AlignmentSplitterGeneric(alignments, filter);
}

// overloaded function: call with int to count reads with at least ESD_cutoff
// edited non-reference positions
// or call with float to count reads with at least ESD_cutoff_p % of
// edited sites (w.r.o.f alignment length)
vector<vector<TAlignment> > AlignmentSplitterByEditDist(
	vector<TAlignment>& alignments, double ESD_cutoff_p)
{
	std::function<bool (TAlignment&)> filter = [=](TAlignment& tal)
	{
		return (static_cast<double>(tal.ref_ESD) / static_cast<double>(tal.dR.size()) >= ESD_cutoff_p);
	};
	return AlignmentSplitterGeneric(alignments, filter);
}

vector<vector<TAlignment> > AlignmentSplitterByEditDist(
	vector<TAlignment>& alignments, int ESD_cutoff)
{
	std::function<bool (TAlignment&)> filter = [=](TAlignment& tal)
	{
		return (tal.ref_ESD >= ESD_cutoff);
	};
	return AlignmentSplitterGeneric(alignments, filter);
}

vector<vector<int> > CalcTalMatrix(vector<TAlignment>& alignments,
	ReadAlignmentTask& preRAT, int refID)
{
	const string& RfTlp = (*(preRAT.refTless))[refID].T;
	vector<vector<int> > rMatrix(RfTlp.size(), vector<int>(40, 0));

	for(size_t a = 0; a < alignments.size(); a++)
	{
		auto& A = alignments[a];
		for(size_t i = 0; i < A.dR.size(); i++)
		{
			if(A.dR[i] > -17 && A.dR[i] < 17 && A.ref_S + i < rMatrix.size())
				rMatrix[A.ref_S + i][20-A.dR[i]] += 1;
		}
	}

	return rMatrix;
}


#endif
