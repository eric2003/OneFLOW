# -*- coding: UTF-8 -*-
import os

def WriteMainFile():
    file = open("test.cpp", "w")
    tab = "\t"
    line_end = "\r\n"

    file.write( "haha\n" )
    file.write( line_end )
    file.close()
    
def WriteFile1():
    file = open("test1.cpp", "w")
    file.write( "haha" )
    file.close()
    
def WriteFile2():
    file = open("test2.cpp", "w")
    file.write( "haha\n" )
    file.close()
    
def ReadFile1():
    file = open("test1.cpp", "r")
    count = 0
    while True:
        str = file.readline()
        print( "str = " , str )
        count += 1
        print( " count = ", count )
        #if len( str ) == 0:
        if count > 3:
            break
        
    file.close()
    
def ReadFile2():
    file = open("test2.cpp", "r")
    count = 0
    while True:
        str = file.readline()
        print( "str = " , str )
        count += 1
        print( " count = ", count )
        #if len( str ) == 0:
        if count > 3:
            break
        
    file.close()   

def ReadFile3():
    file = open("test3.cpp", "r")
    count = 0
    while True:
        str = file.readline()
        print( "str = " , str )
        count += 1
        print( " count = ", count )
        line_len = len( str )
        print( " line_len = ", line_len )
        if str == '\n' :
            print( " str = change line \n"  )
        if line_len == 0:
            break
        if count > 5:
            break
        
    file.close()
    
def WriteNewLine( filename ):
    file = open( filename, "r+" )
    newline = 0
    last_line = ""
    len_0 = len( last_line )
    while True:
        str = file.readline()
        #print( "str = " , str )
        lenstr = len( str )
        if lenstr == 0:
            len_0 = len( last_line )
            if len_0 == 0:
                file.write('\n')
            else:
                str_end = last_line[-1:]
                if str_end == '\n':
                    break
                else:
                    file.write('\n')
            break
        last_line = str
    file.close() 

def ModifyNewLineForAllFiles( base ):
    numfile = 0
    for root, dirs, files in os.walk( base ):
        for file in files:
            fullname = os.path.join( root, file )
            print( "file = ", fullname, "\n" )
            WriteNewLine( fullname )
            numfile += 1
    print( " total file number = ", numfile )

    

if __name__ == '__main__':
    dir = "/home/kylinlu/work_gitlab/oneflow/OneFLOW_deb/codes"
    ModifyNewLineForAllFiles( dir )
