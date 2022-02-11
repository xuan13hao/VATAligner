

#ifndef VALUE_TYPE_H_
#define VALUE_TYPE_H_

struct Letter_prot { };
struct Letter_nucl { };

template<typename _tag>
struct Value_type
{

	Value_type():
		value_ ()
	{ }

	Value_type(int v):
		value_ ((char)v)
	{ }

	Value_type(char v):
		value_ (v)
	{ }

	Value_type(unsigned v):
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

	bool operator==(const Value_type &rhs) const
	{ return value_ == rhs.value_; }

	bool operator!=(const Value_type &rhs) const
	{ return value_ != rhs.value_; }

private:

	char value_;

};

typedef Value_type<Letter_prot> Protein;
typedef Value_type<Letter_nucl> DNA;

#endif /* VALUE_TYPE_H_ */
