#ifndef __SEQUENCETYPE_H__
#define __SEQUENCETYPE_H__


class SequenceProtein
{

};

class SequenceDNA
{

};

class SequenceRNA
{

};

template<typename _t>
class SequenceTypes
{
    public:
    SequenceTypes():value_()
    {

    }

    SequenceTypes(int i)
    {
        value_ = ((char)i);
    }

    SequenceTypes(char str)
    {
        value_ = str;
    }

    SequenceTypes(unsigned i)
    {
        value_ = ((char)i);
    }

    operator char() const
    {
        return value_;
    }

    operator int() const
    {
        return (int)value_;
    }

    operator unsigned() const
    {
        return (unsigned)value_;
    }

    operator long() const
    {
        return (long)value_;
    }

    bool operator ==(const SequenceTypes &vt) const
    {
        return value_ == vt.value_;
    }

    bool operator!=(const SequenceTypes &vt) const
    {
        return value_ != vt.value_;
    }

    private:
        char value_;
};

typedef SequenceTypes<SequenceDNA> DNA;
typedef SequenceTypes<SequenceProtein> Protein;
typedef SequenceTypes<SequenceRNA> RNA;

#endif // __SEQUENCETYPE_H__