

#ifndef VALUE_TYPE_H_
#define VALUE_TYPE_H_

class AlphatProtein { };
class AlphatDNA { };

template<typename _tag>
class AlphabetType
{
	public:
	AlphabetType():
		value_ ()
	{ }

	AlphabetType(int v):
		value_ ((char)v)
	{ }

	AlphabetType(char v):
		value_ (v)
	{ }

	AlphabetType(unsigned v):
		value_ ((char)v)
	{ }

	operator char() const
	{ return value_; }

	operator int() const
	{ return (int)value_; }

	operator unsigned() const
	{ return (unsigned)value_; }

	operator long() const
	{ return (long)value_; }

	bool operator==(const AlphabetType &rhs) const
	{ return value_ == rhs.value_; }

	bool operator!=(const AlphabetType &rhs) const
	{ return value_ != rhs.value_; }

private:

	char value_;

};

typedef AlphabetType<AlphatProtein> Protein;
typedef AlphabetType<AlphatDNA> DNA;

#endif /* VALUE_TYPE_H_ */
