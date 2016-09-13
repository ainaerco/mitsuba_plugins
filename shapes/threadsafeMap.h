#pragma once
#include <boost/thread.hpp>
#include <boost/atomic.hpp>

template <typename Item, typename Mutex>
class ReaderProxy
{
public:
	ReaderProxy(Item& i, Mutex& m) : lock(m), item(i) {}

	Item* operator->() { return &item; }

private:
	boost::shared_lock<Mutex> lock;
	Item& item;
};

template <typename Item, typename Mutex>
class WriterProxy
{
public:
	WriterProxy(Item& i, Mutex& m) : uplock(m), lock(uplock), item(i) {}

	Item* operator->() { return &item; }

private:
	boost::upgrade_lock<Mutex> uplock;
	boost::upgrade_to_unique_lock<Mutex> lock;
	Item& item;
};

template<class Key, class Data>
class MapSafe
{
	typedef ReaderProxy< std::map<Key, Data>, boost::shared_mutex> Reader;
	typedef WriterProxy< std::map<Key, Data>, boost::shared_mutex> Writer;

public:
	Data find(const Key& k)
	{
		Reader r(map, m);
		auto it = r->find(k);
		if (it == r->end()) { return NULL; }
		return it->second;
	}

	void insert(const Key& k, Data v)
	{
		
		Writer w(map, m);
		w->insert(std::make_pair(k, v));
		
	}
	void clear()
	{
		map.clear();
	}
private:
	boost::shared_mutex m;
	std::map<Key, Data> map;
};