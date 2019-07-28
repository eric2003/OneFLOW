#-*- coding:utf-8 -*-
'''
---------------------------------------------------------------------------
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

---------------------------------------------------------------------------
'''
import os
import optparse
import shutil
import sys
import datetime
import subprocess
import time
import filecmp

def my_dir_cmd( file_dir ):
    dirs = os.listdir( file_dir )
    for file in dirs:
       print( file )

def RunTest(testprjdir):
    print("testprjdir=",testprjdir)
    testScript = testprjdir + "/autotest/test.txt"
    absTestScript = os.path.normpath(os.path.abspath(testScript))
    print("testScript =", testScript)
    print("absTestScript =", absTestScript)
    fileNameList = []
    with open(absTestScript, 'r', encoding='utf-8-sig') as f:
        for line in f.readlines():
            print('line=', line)
            ll = line.replace("\r", "")
            ll = ll.replace("\n", "")
            fileNameList.append(ll)
    f.close()
    print(fileNameList)
    ffList = []
    for file in fileNameList:
        print(file)
        ff = file.split("/")[-1]
        ffList.append(ff)
    print(ffList)
    testFileList = []
    testFileListPath=[]
    for ff in ffList:
        ffNew = testprjdir + "/autotest/" + ff
        testFileList.append(ffNew)
        ffNewPath = os.path.normpath(os.path.abspath(ffNew))
        testFileListPath.append(ffNewPath)
    print(testFileListPath)
    fileNameListPath = []
    for file in fileNameList:
        fileNew = testprjdir + "/" + file
        fileNewPath = os.path.normpath(os.path.abspath(fileNew))
        fileNameListPath.append(fileNewPath)
    print(fileNameListPath)

    for i in range(0, len(testFileListPath)):
        print("i=", i, " var=", testFileListPath[i], "file=", fileNameListPath[i])

    exedir = '"c:/Program Files (x86)/OneFLOW/bin/"'
    #cmd = exedir +"OneFLOW" + " " + testprjdir
    cmd = "mpiexec -n 1 " + exedir +"OneFLOW" + " " + testprjdir
    print(cmd)
    process = subprocess.Popen(cmd, shell=True)
    while process.poll() is None:
        current = datetime.datetime.now()
        time.sleep(0.5)
    returnCode = process.poll()
    for i in range(0, len(testFileListPath)):
        print("i=", i)
        print(" file1=", testFileListPath[i])
        print(" file2=", fileNameListPath[i])
        flag = filecmp.cmp(testFileListPath[i], fileNameListPath[i])
        print("flag=", flag)
    print("returnCode=",returnCode)

def RunAllTest(filename):
    with open(filename, 'r', encoding='utf-8-sig') as f:
        for line in f.readlines():
            print(os.getcwd())
            print('line=', line)
            dir = "workdir/" + line
            print("dir=",dir)
            RunTest(dir)
    f.close()

def main():
    totalPass = True
    location = os.getcwd()
    print( " location = ", location )
    my_dir_cmd( location )

    RunAllTest("test.txt")
    # print( 'number of param = ', len(sys.argv) )
    # print( 'list of param = ', str(sys.argv) )
    # nargv = len(sys.argv);
    # if nargv > 2 :
    #     name = sys.argv[ 1 ];
    #     ccc  = sys.argv[ 2 ];
    #     print( "name = ", name )
    #     os.chdir( name )
    #     print(" current dir is = ", os.getcwd())
    #     dd = os.getcwd()
    #     my_dir_cmd( dd )
    #     ee = dd + "/"
    #     ff = dd + "/" + ccc
    #     print(ee)
    #     print(ff)
    #     my_dir_cmd(ff)
    # print( " ccc = ", ccc )
    # path = "tmp/" + ccc;
    # #os.makedirs(path)
    #
    # shutil.copytree(ccc, path)

    errorCode = 0
    if totalPass:
        print("All tests passed!")
    else:
        print("ERROR: Some tests failed")
        errorCode = 1
    sys.exit(errorCode)


if __name__ == "__main__":
    main()
