#pragma once
#include <boost/thread.hpp>
#include <boost/asio.hpp>

typedef std::unique_ptr<boost::asio::io_service::work> asio_worker;

struct ThreadPool {
    ThreadPool(size_t threads) :service(), working(new asio_worker::element_type(service)) {
        while(threads--)
        {
            auto worker = boost::bind(&boost::asio::io_service::run, &(this->service));
            g.add_thread(new boost::thread(worker));
        }
    }

    template<class F>
        void enqueue(F f){
            service.post(f);
        }

    ~ThreadPool() {
        working.reset(); //allow run() to exit
        g.join_all();
        service.stop();
    }

    private:
    boost::asio::io_service service; //< the io_service we are wrapping
    asio_worker working;
    boost::thread_group g; //< need to keep track of threads so we can join them
};