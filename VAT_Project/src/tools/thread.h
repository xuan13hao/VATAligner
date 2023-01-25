
#ifndef THREAD_H_
#define THREAD_H_

#include <vector>
#include "fast_mutex.h"
#include "tinythread.h"

using tthread::thread;
using std::vector;

template<typename _t>
struct Atomic
{
	Atomic(const _t &v):
		v_ (v)
	{ }
	volatile _t operator++(int)
	{
		tthread::lock_guard<tthread::mutex> lock (mtx_);
		_t r = v_++;
		return r;
	}
private:
	volatile _t v_;
	tthread::mutex mtx_;
};

template<typename _context>
struct Thread_p
{
	Thread_p(unsigned thread_id, _context &context):
		thread_id (thread_id),
		context (&context)
	{ }
	unsigned thread_id;
	_context *context;
};

template<typename _context>
void pool_worker(void *p)
{
	((Thread_p<_context>*)p)->context->operator()(((Thread_p<_context>*)p)->thread_id);
}

template<typename _context>
void launch_thread_pool(_context &context, unsigned threads)
{
	vector<tthread::thread*> t;
	vector<Thread_p<_context> > p;
	p.reserve(threads);
	unsigned n = 0;
	for(unsigned i=0;i<threads;++i) {
		p.push_back(Thread_p<_context> (i, context));
		t.push_back(new tthread::thread(pool_worker<_context>, (void*)&p.back()));
		n += t.back()->get_id() == tthread::thread::id () ? 0 : 1;
	}
	for(vector<tthread::thread*>::iterator i=t.begin();i!=t.end();++i) {
		(*i)->join();
		delete *i;
	}
	exception_state.sync();
	if(n != threads)
		throw std::runtime_error("Failed to create thread.");
}

template<typename _context>
struct Schedule_context
{
	Schedule_context(_context &context, unsigned count):
		context (context),
		n (0),
		count (count)
	{ }
	void operator()(unsigned thread_id)
	{
		try {
			unsigned idx;
			while(!exception_state() && (idx = n++) < count)
				context(thread_id, idx);
		} catch(std::exception &e) {
			exception_state.set(e);
		}
	}
	_context &context;
	Atomic<unsigned> n;
	const unsigned count;
};

template<typename _context>
void launch_scheduled_thread_pool(_context &context, unsigned count, unsigned threads)
{
	Schedule_context<_context> c (context, count);
	launch_thread_pool(c, threads);
}

template<typename _f, typename _t1, typename _t2, typename _t3, typename _t4>
struct Thread_p4
{
	Thread_p4(_f f, _t1 p1, _t2 p2, _t3 p3, _t4 p4):
		f (f),
		p1 (p1),
		p2 (p2),
		p3 (p3),
		p4 (p4)
	{ }
	_f f;
	_t1 p1;
	_t2 p2;
	_t3 p3;
	_t4 p4;
};

template<typename _f, typename _t1, typename _t2, typename _t3, typename _t4>
void thread_worker(void *p)
{
	Thread_p4<_f,_t1,_t2,_t3,_t4> *q = (Thread_p4<_f,_t1,_t2,_t3,_t4>*)p;
	q->f(q->p1, q->p2, q->p3, q->p4);
	delete q;
}

template<typename _f, typename _t1, typename _t2, typename _t3, typename _t4>
thread* launch_thread(_f f, _t1 p1, _t2 p2, _t3 p3, _t4 p4)
{ return new thread (thread_worker<_f,_t1,_t2,_t3,_t4>, new Thread_p4<_f,_t1,_t2,_t3,_t4> (f, p1, p2, p3, p4)); }

#endif /* THREAD_H_ */
