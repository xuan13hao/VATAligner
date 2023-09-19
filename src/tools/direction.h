

#ifndef DIRECTION_H_
#define DIRECTION_H_

struct Right { };
struct Left { };

template<typename _val>
_val get_dir(const _val* x, int i, const Right&)
{ return *(x+i); }

template<typename _val>
_val get_dir(const _val* x, int i, const Left&)
{ return *(x-i); }

template<typename _val>
const _val* get_dir_ptr(const _val* x, int i, const Right&)
{ return x+i; }

template<typename _val>
const _val* get_dir_ptr(const _val* x, int i, const Left&)
{ return x-i; }

template<typename _val>
const _val* inc_dir(const _val* x, const Right&)
{ return x+1; }

template<typename _val>
const _val* inc_dir(const _val* x, const Left&)
{ return x-1; }

#endif /* DIRECTION_H_ */
