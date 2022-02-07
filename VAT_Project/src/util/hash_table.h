

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

	struct entry
	{
		_K	key;
		_V	value;
	} __attribute__((packed));

	hash_table(size_t size):
		table (new entry[size]),
		size_ (size)
	{ memset(table, 0, size_ * sizeof(entry)); }

	~hash_table()
	{ delete[] table; }

	entry* operator[](_K key) const
	{
		entry *entry = get_entry(key);
		if(_E()(entry->value))
			return NULL;
		return entry;
	}

	void insert(_K key, _V value)
	{
		entry *entry = get_entry(key);
		if(_E()(entry->value))
			entry->key = key;
		entry->value = value;
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

	entry* get_entry(_K key) const
	{
		entry *p = &table[_H()(key) % size_];
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

	entry	*table;
	size_t	 size_;

};

#endif /* HASH_TABLE_H_ */
