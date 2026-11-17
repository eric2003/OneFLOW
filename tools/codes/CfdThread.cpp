#include "CfdThread.h"
#include <QTime>
#include <thread>
#include <iostream>

CfdThread::CfdThread() 
{
    this->paused = false;
    this->stop = false;
}

void CfdThread::beginCFD()
{
    this->paused = false;

}

void CfdThread::pauseCFD()
{
    this->paused = true;
}

void CfdThread::stopCFD()
{
    this->stop = true;
}

void CfdThread::run()
{
    download( "CfdThread::run()" );
}

void CfdThread::download(const std::string &file)
{
    for ( int i = 0; i < 10; ++ i )
    {
        if ( this->stop ) break;
        std::cout << "Downloading " << file
            << " (" << i * 10 << "%)..." << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
    std::cout << "Download complete: " << file << std::endl;
}
