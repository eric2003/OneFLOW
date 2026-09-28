#pragma once
#include <fstream>

void OpenLogFile( int logFileIndex, std::fstream & file );
void CloseLogFile( std::fstream & file );
class LogFile;
extern LogFile logFile;

class Parallel
{
public:
    static int pid;
};

class LogFile
{
public:
    LogFile();
    ~LogFile();
    std::fstream my_fstream;
public:
    void Open();
    void Close();
};

template< typename T >
LogFile & operator << ( LogFile & f, const T & value )
{
    f.Open();
    f.my_fstream << value;
    f.Close();
    return f;
}
