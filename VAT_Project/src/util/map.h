
#ifndef MAP_H_
#define MAP_H_

template<typename _it, typename _key>
struct Map
{

	struct Iterator
	{
		Iterator(const _it& begin, const _it& parent_end):
			begin_ (begin),
			parent_end_ (parent_end),
			end_ (get_end())
		{ }
		void operator++()
		{ begin_ = end_; end_ = get_end(); }
		bool valid() const
		{ return begin_ < parent_end_; }
		_it& begin()
		{ return begin_; }
		_it& end()
		{ return end_; }
	private:
		_it get_end() const
		{
			_it i = begin_;
			while(i < parent_end_ && _key()(*i) == _key()(*begin_)) ++i;
			return i;
		}
		_it begin_, parent_end_, end_;
	};

	Map(const _it &begin, const _it &end):
		begin_ (begin),
		end_ (end)
	{ }

	Iterator begin()
	{ return Iterator(begin_, end_); }

private:
	_it begin_, end_;

};

#endif /* MAP_H_ */
