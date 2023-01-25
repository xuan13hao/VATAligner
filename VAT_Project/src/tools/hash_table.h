
#ifndef HASH_TABLE_H_
#define HASH_TABLE_H_

template<typename T, T value> struct value_compare
{
	bool operator()(T x) const
	{ return x == value; }
};

template<typename _K, typename _V, typename _E, typename _H> class hash_table
{

public:

	struct Tuple
	{
		_K	key;
		_V	value;
	} __attribute__((packed));

	hash_table(size_t size):
		table (new Tuple[size]),
		size_ (size)
	{ memset(table, 0, size_ * sizeof(Tuple)); }

	~hash_table()
	{ delete[] table; }

	Tuple* operator[](_K key) const
	{
		Tuple *Tuple = get_entry(key);
		if(_E()(Tuple->value))
			return NULL;
		return Tuple;
	}

	void insert(_K key, _V value)
	{
		Tuple *Tuple = get_entry(key);
		if(_E()(Tuple->value))
			Tuple->key = key;
		Tuple->value = value;
	}

	size_t size() const
	{
		return size_;
	}

	size_t count() const
	{
		size_t n (0);
		for(size_t i=0;i<size_;++i)
			if(!_E()(table[i].value))
				++n;
		return n;
	}

private:

	Tuple* get_entry(_K key) const
	{
		Tuple *p = &table[_H()(key) % size_];
		bool wrapped = false;
		while(p->key != key && !_E()(p->value)) {
			++p;
			if(p == &table[size_]) {
				if(wrapped)
					throw hash_table_overflow_exception();
				p = &table[0];
				wrapped = true;
			}
		}
		return p;
	}

	Tuple	*table;
	size_t	 size_;

};

#endif /* HASH_TABLE_H_ */
