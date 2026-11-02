#include "LogFile.h"
#include "Parallel.h"
#include <sstream>
#include <iostream>

#ifdef _WINDOWS
#include <windows.h>
#include <direct.h>
#include <io.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

LogFile logFile;

class OStream : public std::ostringstream
{
public:
    OStream() {}
    ~OStream() {}
public:
    void ClearAll()
    {
        this->clear();
        this->str("");
    }
};

OStream StrIO;
std::string GetPrjDirName( const std::string & fileName );
bool DirExist( const std::string & dirName );
void MakeDir( const std::string & dirName );

std::string GetPrjDirName( const std::string & fileName )
{
    size_t pos = fileName.find_last_of("\\/");
    if ( std::string::npos == pos )
    {
        return "";
    }
    else
    {
        return fileName.substr( 0, pos );
    }
}

bool DirExist( const std::string & dirName )
{
#ifdef _WINDOWS
    bool flag = ( _access( dirName.c_str(), 0 ) == 0 );
    return flag;
#else
    bool flag = ( access( dirName.c_str(), 0 ) == 0 );
    return flag;
#endif
}

void MakeDir( const std::string & dirName )
{
    int flag;
#ifdef _WINDOWS
    flag = _mkdir( dirName.c_str() );
#else
    flag = mkdir( dirName.c_str(), S_IRWXU );
#endif
    if ( flag == 0 )
    {
        std::cout << dirName << " directory has been created successfully !\n";
    }
}

void CreateDirIfNeeded( std::string & prjFileName )
{
    std::string prj_dir = GetPrjDirName( prjFileName );

    if ( ! DirExist( prj_dir ) )
    {
        MakeDir( prj_dir );
    }
}


void OpenLogFile( int logFileIndex, std::fstream & file )
{
    static int ifReWrite = 0;

    StrIO.ClearAll();
    StrIO << "log/log" << logFileIndex << ".log";
    std::string fileName = StrIO.str();

    std::ios_base::openmode openMode;

    if ( ifReWrite == 0 )
    {
        CreateDirIfNeeded( fileName );

        openMode = std::ios_base::out | std::ios_base::trunc;

        ifReWrite = 1;
    }
    else
    {
        openMode = std::ios_base::out | std::ios_base::app;
    }

    file.open( fileName.c_str(), openMode );
    if ( ! file )
    {
        std::cout << "could not open " << fileName << std::endl;
        exit( 0 );
    }
}

void CloseLogFile( std::fstream & file )
{
    file.close();
    file.clear();
}

LogFile::LogFile()
{
}

LogFile::~LogFile()
{
}

void LogFile::Open()
{
    int pid = Parallel::pid;
    OpenLogFile( pid, this->my_fstream );
}

void LogFile::Close()
{
    CloseLogFile( this->my_fstream );
}
