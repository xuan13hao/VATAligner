

#ifndef TRANSLATE_H_
#define TRANSLATE_H_

class Translator
{

public:

	static const DNA reverseNucleotide[5];
	static const Protein reverseProtein[5];
	static const Protein lookup[5][5][5];
	static const Protein lookupReverse[5][5][5];
	static const Protein STOP;

	static DNA getReverseComplement(DNA nucleotide)
	{ return reverseNucleotide[(int) nucleotide]; }
	static Protein getReverseComplement(Protein nucleotide)
	{ return reverseProtein[(int) nucleotide]; }
	static Protein getAminoAcid(vector<DNA> const &dnaSequence, size_t pos)
	{ return lookup[(int) dnaSequence[pos]][(int)dnaSequence[pos+1]][(int)dnaSequence[pos+2]]; }

	static Protein getAminoAcidReverse(vector<DNA> const &dnaSequence, size_t pos)
	{ return lookupReverse[(int) dnaSequence[pos + 2]][(int)dnaSequence[pos + 1]][(int)dnaSequence[pos]]; }

	static vector<DNA> reverse(const vector<DNA> &seq)
	{
		vector<DNA> r;
		for(vector<DNA>::const_reverse_iterator i=seq.rbegin(); i!=seq.rend(); ++i)
			r.push_back(getReverseComplement(*i));
		// std::reverse(r.begin(),r.end());
		return r;
	}
	static vector<Protein> reverse(const vector<Protein> &seq)
	{
		vector<Protein> r;
		for(vector<Protein>::const_reverse_iterator i=seq.rbegin(); i!=seq.rend(); ++i)
			r.push_back(getReverseComplement(*i));
		return r;
	}
	static size_t translate(vector<DNA> const &dnaSequence, vector<Protein> *proteins)
	{
		size_t length_ = dnaSequence.size(), d, n;
		proteins[0].resize(d = length_ / 3);
		proteins[3].resize(d);
		n = 2*d;
		proteins[1].resize(d = (length_-1) / 3);
		proteins[4].resize(d);
		n += 2*d;
		proteins[2].resize(d = (length_-2) / 3);
		proteins[5].resize(d);
		n += 2*d;

		size_t r = length_ - 2;
		unsigned pos = 0;
		unsigned i = 0;
		while(r > 2) {
			proteins[0][i] = getAminoAcid(dnaSequence, pos++);
			proteins[3][i] = getAminoAcidReverse(dnaSequence, --r);
			proteins[1][i] = getAminoAcid(dnaSequence, pos++);
			proteins[4][i] = getAminoAcidReverse(dnaSequence, --r);
			proteins[2][i] = getAminoAcid(dnaSequence, pos++);
			proteins[5][i] = getAminoAcidReverse(dnaSequence, --r);
			++i;
		}
		if(r) {
			proteins[0][i] = getAminoAcid(dnaSequence, pos++);
			proteins[3][i] = getAminoAcidReverse(dnaSequence, --r);
		}
		if(r) {
			proteins[1][i] = getAminoAcid(dnaSequence, pos);
			proteins[4][i] = getAminoAcidReverse(dnaSequence, r);
		}
		return n;
	}

	static Protein const* nextChar(Protein const*p, Protein const*end)
	{
		while(*(p) != STOP && p < end)
			++p;
		return p;
	}

	static void mask_runs(vector<Protein> &query, unsigned run_len)
	{
		Protein *last = &query[0]-1, *i = &query[0], *end = &query.back();
		while (i <= end) {
			if(*i == STOP) {
				if(last != 0 && i - last - 1 < run_len) {
					for(Protein *j = last+1; j < i; ++j)
						*j = 23;
				}
				last = i;
			}
			++i;
		}
		if(last != 0 && i - last - 1 < run_len) {
			for(Protein *j = last+1; j < i; ++j)
				*j = 23;
		}
	}

	static unsigned computeGoodFrames(vector<Protein>  const *queries, unsigned runLen)
	{
		unsigned set = 0;

		for (unsigned i = 0; i < 6; ++i) {
			if (queries[i].size() > 0) {
				unsigned run = 0;
				Protein const*p =  &(queries[i][0]);
				Protein const*q;
				Protein const*end = p + queries[i].size();
				while((q = nextChar(p, end)) < end) {
					run = q-p;
					if (run >= runLen)
						set |= 1 << i;
					p=q+1;
				}
				run = q-p;
				if (run >= runLen)
					set |= 1 << i;
			}
		}
		return set;
	}

	static void mask_runs(vector<Protein> *queries, unsigned run_len)
	{
		for (unsigned i = 0; i < 6; ++i)
			mask_runs(queries[i], run_len);
	}

};

const DNA Translator::reverseNucleotide[5] = { 3, 2, 1, 0, 4 };
const Protein Translator::reverseProtein[5] = { 3, 2, 1, 0, 4 };
const Protein Translator::lookup[5][5][5] = {
{ { 11,2,11,2,23 },
{ 16,16,16,16,16 },
{ 1,15,1,15,23 },
{ 9,9,12,9,23 },
{ 23,23,23,23,23 },
 },
{ { 5,8,5,8,23 },
{ 14,14,14,14,14 },
{ 1,1,1,1,1 },
{ 10,10,10,10,10 },
{ 23,23,23,23,23 },
 },
{ { 6,3,6,3,23 },
{ 0,0,0,0,0 },
{ 7,7,7,7,7 },
{ 19,19,19,19,19 },
{ 23,23,23,23,23 },
 },
{ { 23,18,23,18,23 },
{ 15,15,15,15,15 },
{ 23,4,17,4,23 },
{ 10,13,10,13,23 },
{ 23,23,23,23,23 },
 },
{ { 23,23,23,23,23 },
{ 23,23,23,23,23 },
{ 23,23,23,23,23 },
{ 23,23,23,23,23 },
{ 23,23,23,23,23 },
} };

const Protein Translator::lookupReverse[5][5][5] = {
{ { 13,10,13,10,23 },
{ 4,17,4,23,23 },
{ 15,15,15,15,15 },
{ 18,23,18,23,23 },
{ 23,23,23,23,23 },
 },
{ { 19,19,19,19,19 },
{ 7,7,7,7,7 },
{ 0,0,0,0,0 },
{ 3,6,3,6,23 },
{ 23,23,23,23,23 },
 },
{ { 10,10,10,10,10 },
{ 1,1,1,1,1 },
{ 14,14,14,14,14 },
{ 8,5,8,5,23 },
{ 23,23,23,23,23 },
 },
{ { 9,12,9,9,23 },
{ 15,1,15,1,23 },
{ 16,16,16,16,16 },
{ 2,11,2,11,23 },
{ 23,23,23,23,23 },
 },
{ { 23,23,23,23,23 },
{ 23,23,23,23,23 },
{ 23,23,23,23,23 },
{ 23,23,23,23,23 },
{ 23,23,23,23,23 },
}};

const Protein Translator::STOP (AlphabetAttributes<Protein>::from_char('*'));

#endif /* TRANSLATE_H_ */
