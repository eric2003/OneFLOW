Cmd
==================================

Cmd
---------------------------------
#. `CMD相关链接 <https://zhuanlan.zhihu.com/p/662224432/>`_
#. `Windows batch script系列链接（不定期更新） <https://zhuanlan.zhihu.com/p/509954303/>`_
#. `An A-Z Index of Windows CMD commands <https://ss64.com/nt/>`_
#. `ModernCmd <https://github.com/eric2003/ModernCmd/>`_
#. `Windows Commands <https://learn.microsoft.com/en-us/windows-server/administration/windows-commands/windows-commands/>`_
#. `Various batch files for Windows <https://github.com/Espionage724/Windows/>`_
#. `Collection of Batch scripts, examples <https://github.com/happy05dz/Batch-Script-Collection/>`_
#. `Windows Batch Scripting <https://en.wikibooks.org/wiki/Windows_Batch_Scripting>`_
#. `Batch Script Tutorial <https://www.tutorialspoint.com/batch_script/index.htm>`_



Code
--------------------
::

  git clone git@github.com:eric2003/ModernCmd.git
  or 
  git clone https://github.com/eric2003/ModernCmd.git
  
  
You can edit your Microsoft.VSCode_profile.ps1 file (full path can be found when you type $profile in the editor and remove all code referring to anaconda. To open the profile file, type notepad $profile in the terminal. 

Display all environment variables
::

  SET

Display ProgramData variable
::

  echo %ProgramData%
  
 
Display ProgramFiles(x86) variable
::

  echo %ProgramFiles(x86)%
  
Display ProgramFiles(x86) variable
::

  echo %ProgramFiles(x86)%
  
::

  echo %comspec%
  results:
  C:\Users\eric>echo %comspec%
  C:\WINDOWS\system32\cmd.exe
  
set (environment variable)
::
  
  set testVar=TEST^&1
  set testVar
  
delete environment variable
::

   set testVar=
   
Run Vs2022 bat   
:: 

  C:\Users\eric>"c:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\"VsDevCmd
  **********************************************************************
  ** Visual Studio 2022 Developer Command Prompt v17.7.5
  ** Copyright (c) 2022 Microsoft Corporation
  **********************************************************************
  
  C:\Users\eric>echo %LIB%
  c:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.37.32822\ATLMFC\lib\x86;c:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.37.32822\lib\x86;C:\Program Files (x86)\Windows Kits\NETFXSDK\4.8\lib\um\x86;C:\Program Files (x86)\Windows Kits\10\lib\10.0.19041.0\ucrt\x86;C:\Program Files (x86)\Windows Kits\10\\lib\10.0.19041.0\\um\x86  
  
FOR
-------------------------
::

  @echo off
  for %%i in (1,2,3) do echo %%i
  pause  
  results:
  1
  2
  3  

::

  @echo off
  
  for %%i in (a b c d) do (
    echo %%i
  )
  pause  
  results:
  a
  b
  c
  d  
  
使用for /f命令从文件中读取每一行内容并进行处理：
::

  @echo off
  
  for /f "tokens=*" %%i in (myfile.txt) do (
    echo %%i
  )
  pause
  results:
  1
  2 3 4
  5 a
  
myfile.txt
::

  1
  2 3 4
  5 a  
  
使用for /f命令遍历文件夹中的所有文件：  
::

  @echo off
  
  for /f "tokens=*" %%i in ('dir /b') do (
    echo %%i
  )
  pause
  results:
  d:\work\batch_work\ModernBatchFiles\codes\for\05>testprj.bat
  1.txt
  2.txt
  3.txt
  testprj.bat
  请按任意键继续. . .  
  
vswhere -property installationPath
::

  @echo off
  
  for /f "delims=" %%a in (
  '"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -property installationPath'
  ) do (
      echo %%a
  )
  pause
  results:
  d:\work\batch_work\ModernBatchFiles\codes\for\07>testprj.bat
  C:\Program Files\Microsoft Visual Studio\2022\Community
  请按任意键继续. . .  




