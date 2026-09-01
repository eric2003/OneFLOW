#-*- coding:utf-8 -*-
'''
---------------------------------------------------------------------------
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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

from residual_db import compare_case

def my_dir_cmd(file_dir):
    dirs = os.listdir(file_dir)
    for file in dirs:
       print(file)
def CmpFloatStr(str1, str2):
    ss1 = str1.split()
    ss2 = str2.split()
    passflag = True
    for i in range(0, len(ss1)):
        a = ss1[i]
        b = ss2[i]
        flag = is_number(a)
        if not flag:
            continue
        fa = float(a)
        fb = float(b)
        if abs(fa-fb) > 1.0e-8:
            passflag = False
            return passflag
    return passflag
def CmpFile(fileName1, fileName2):
    f1 = open(fileName1, 'r', encoding='utf-8-sig')
    f2 = open(fileName2, 'r', encoding='utf-8-sig')
    cmpflag = True
    line_id = 0
    while True:
        str1 = f1.readline()
        str2 = f2.readline()
        if not str1 and not str2:
            break
        if not str1 or not str2:
            print(fileName1)
            print(fileName2)
            print("Line", line_id, "has different end-of-file position")
            cmpflag = False
            break
        # aa = str1.split()
        # print(aa)
        if str1 != str2:
            flag = CmpFloatStr(str1, str2)
            if not flag:
                print(fileName1)
                print(fileName2)
                print("Line", line_id)
                print('line1=', str1)
                print('line2=', str2)
                cmpflag = False
                break
        line_id += 1
    f1.close()
    f2.close()
    return cmpflag

def GetFileNameList(filename, filename_list):
    with open(filename, 'r', encoding='utf-8-sig') as f:
        for line in f.readlines():
            fname = line.strip()
            filename_list.append(fname)
    f.close()
def RunTest(testprjdir, residual_db_path=None, residual_tolerance=None):
    print("testprjdir=",testprjdir)
    testScript = testprjdir + "/autotest/test.txt"
    absTestScript = os.path.normpath(os.path.abspath(testScript))
    print("testScript =", testScript)
    print("absTestScript =", absTestScript)
    fileNameList = []
    GetFileNameList(absTestScript,fileNameList)
    #print(fileNameList)
    ffList = []
    for file in fileNameList:
        #print("file=", file)
        ff = file.split("/")[-1]
        ffList.append(ff)
    #print(ffList)
    testFileList = []
    testFileListPath=[]
    for ff in ffList:
        ffNew = testprjdir + "/autotest/" + ff
        testFileList.append(ffNew)
        ffNewPath = os.path.normpath(os.path.abspath(ffNew))
        testFileListPath.append(ffNewPath)
    #print(testFileListPath)
    fileNameListPath = []
    for file in fileNameList:
        fileNew = testprjdir + "/" + file
        fileNewPath = os.path.normpath(os.path.abspath(fileNew))
        fileNameListPath.append(fileNewPath)
    #print(fileNameListPath)

    # for i in range(0, len(testFileListPath)):
    #     print("i=", i, " var=", testFileListPath[i], "file=", fileNameListPath[i])

    #exePath = '"c:/Program Files (x86)/OneFLOW/bin/"'
    #cmd = exedir +"OneFLOW" + " " + testprjdir
    lenth = len(sys.argv)
    print("lenth=", lenth)
    for i in range(0, lenth):
        print("i=", i, "var=", sys.argv[i])
    #mpiPath = '"' + sys.argv[1] + '"'
    #exePath = '"' + sys.argv[2] + '"'
    #opsys = sys.argv[3]
    #mpiCmd = "mpiexec -n 1 "
    #if opsys == "linux":
    #    mpiCmd = "mpirun -np 1 "
    #cmd = mpiPath + mpiCmd + exePath + "OneFLOW" + " 0 " + testprjdir
    mpiCmd = sys.argv[ 1 ]
    exeCmd = sys.argv[ 2 ]
    cmd = mpiCmd + " " + exeCmd + " 0 " + testprjdir
    print(cmd)
    process = subprocess.Popen(cmd, shell=True)
    while process.poll() is None:
        current = datetime.datetime.now()
        time.sleep(0.5)
    returnCode = process.poll()
    totalPass = returnCode == 0
    for i in range(0, len(testFileListPath)):
        #print("i=", i)
        #print(" file1=", testFileListPath[i])
        #print(" file2=", fileNameListPath[i])
        #flag = filecmp.cmp(testFileListPath[i], fileNameListPath[i])
        #print("flag=", flag)
        cmpflag = CmpFile(testFileListPath[i], fileNameListPath[i])
        #print("cmpflag=", cmpflag)
        if not cmpflag:
            totalPass = False
            break;
    if residual_db_path:
        residual_result = compare_case(
            testprjdir, residual_db_path, absolute_tolerance=residual_tolerance
        )
        print("residual_baseline=", residual_result)
        if not residual_result["ok"]:
            totalPass = False
    print("returnCode=",returnCode)
    print("totalPass=", totalPass)
    return totalPass

def RunAllTest(filename, residual_db_path=None, residual_tolerance=None):
    passFlag =[]
    with open(filename, 'r', encoding='utf-8-sig') as f:
        for line in f.readlines():
            print(os.getcwd())
            prjname = line.strip()
            print('prjname=', prjname)
            prjdir = prjname
            print("prjdir=",prjdir)
            flag = RunTest(prjdir, residual_db_path, residual_tolerance)
            passFlag.append( flag )
    f.close()
    return passFlag


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    return False

def main():
    location = os.getcwd()
    print( " location = ", location )
    my_dir_cmd( location )
    errorCode = 0
    suiteFile = "suites/cpu-serial.txt"
    if len(sys.argv) >= 4:
        suiteFile = sys.argv[3]
    residualDb = os.environ.get("ONEFLOW_RESIDUAL_DB")
    if len(sys.argv) >= 5:
        residualDb = sys.argv[4]
    residualTolerance = os.environ.get("ONEFLOW_RESIDUAL_TOLERANCE")
    if len(sys.argv) >= 6:
        residualTolerance = sys.argv[5]
    if residualTolerance is not None:
        residualTolerance = float(residualTolerance)
    print("suiteFile=", suiteFile)
    print("residualDb=", residualDb)
    print("residualTolerance=", residualTolerance)
    print("backend=", os.environ.get("ONEFLOW_ACCEL_BACKEND", "CPU"))
    print("residual_profile=", "strict" if os.environ.get("ONEFLOW_RESIDUAL_TEST_OUTPUT", "").lower() in {"1", "true", "on"} else "normal")
    passFlag = RunAllTest(suiteFile, residualDb, residualTolerance)
    numTest = len(passFlag)
    npass = 0
    nfail = 0
    for i in range(0, numTest):
        if passFlag[i]:
            npass = npass + 1
        else:
            nfail = nfail + 1
    if npass == numTest:
        print("Total tests passed!")
    else:
        print("ERROR: Some tests failed")
    print(npass, " tests passed! ", nfail, " tests failed!")
    return 0 if nfail == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
