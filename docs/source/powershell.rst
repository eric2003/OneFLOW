PowerShell
==================================

PowerShell 
---------------------------------
#. `PowerShell相关链接 <https://zhuanlan.zhihu.com/p/481907978/>`_
#. `Manage GitHub with PowerShell <https://www.youtube.com/watch?v=w8tFyophdBA/>`_
#. `PowerShell For Beginners Full Course | PowerShell Beginner tutorial Full Course <https://www.youtube.com/watch?v=UVUd9_k9C6A/>`_


You can edit your Microsoft.VSCode_profile.ps1 file (full path can be found when you type $profile in the editor and remove all code referring to anaconda. To open the profile file, type notepad $profile in the terminal. 

`How to Check the PowerShell Version in Windows 10 <https://www.howtogeek.com/731885/how-to-check-the-powershell-version-in-windows-10/>`_

How to Check PowerShell version 1
--------------------------------------
::

  $PSVersionTable.PSVersion
  result:
  PS C:\Users\eric> $PSVersionTable.PSVersion
  
  Major  Minor  Build  Revision
  -----  -----  -----  --------
  5      1      22621  2506
  
How to Check PowerShell version 2
--------------------------------------
::

  $host.version
  result:
  PS C:\Users\eric> $host.version
  
  Major  Minor  Build  Revision
  -----  -----  -----  --------
  5      1      22621  2506
  
How to Check PowerShell version 3
--------------------------------------
::

  $PSVersionTable
  result:
  PS C:\Users\eric> $PSVersionTable
  
  Name                           Value
  ----                           -----
  PSVersion                      5.1.22621.2506
  PSEdition                      Desktop
  PSCompatibleVersions           {1.0, 2.0, 3.0, 4.0...}
  BuildVersion                   10.0.22621.2506
  CLRVersion                     4.0.30319.42000
  WSManStackVersion              3.0
  PSRemotingProtocolVersion      2.3
  SerializationVersion           1.1.0.1  
  
How to Check PowerShell version 4
--------------------------------------
::

  get-host
  result:
  PS C:\Users\eric> get-host
  
  Name             : ConsoleHost
  Version          : 5.1.22621.2506
  InstanceId       : 50fd2681-bc46-4bbd-b763-2a879fd4f6bb
  UI               : System.Management.Automation.Internal.Host.InternalHostUserInterface
  CurrentCulture   : zh-CN
  CurrentUICulture : zh-CN
  PrivateData      : Microsoft.PowerShell.ConsoleHost+ConsoleColorProxy
  DebuggerEnabled  : True
  IsRunspacePushed : False
  Runspace         : System.Management.Automation.Runspaces.LocalRunspace  
  
How to Check PowerShell version 5
--------------------------------------
::

  host
  result:
  PS C:\Users\eric> host
  
  Name             : ConsoleHost
  Version          : 5.1.22621.2506
  InstanceId       : 50fd2681-bc46-4bbd-b763-2a879fd4f6bb
  UI               : System.Management.Automation.Internal.Host.InternalHostUserInterface
  CurrentCulture   : zh-CN
  CurrentUICulture : zh-CN
  PrivateData      : Microsoft.PowerShell.ConsoleHost+ConsoleColorProxy
  DebuggerEnabled  : True
  IsRunspacePushed : False
  Runspace         : System.Management.Automation.Runspaces.LocalRunspace  
  
Search for the latest version of PowerShell
--------------------------------------------
::

  winget search Microsoft.PowerShell
  result:
  PS C:\Users\eric> winget search Microsoft.PowerShell
  名称               ID                           版本      源
  -----------------------------------------------------------------
  PowerShell         Microsoft.PowerShell         7.4.0.0   winget
  PowerShell Preview Microsoft.PowerShell.Preview 7.4.0.101 winget

Install PowerShell or PowerShell Preview using the id parameter
-------------------------------------------------------------------  
::

  winget install --id Microsoft.Powershell --source winget
  PS C:\Users\eric> winget install --id Microsoft.Powershell --source winget
  已找到 PowerShell [Microsoft.PowerShell] 版本 7.4.0.0
  此应用程序由其所有者授权给你。
  Microsoft 对第三方程序包概不负责，也不向第三方程序包授予任何许可证。
  正在下载 https://github.com/PowerShell/PowerShell/releases/download/v7.4.0/PowerShell-7.4.0-win-x64.msi
    ██████████████████████████████   103 MB /  103 MB
  已成功验证安装程序哈希
  正在启动程序包安装...
  已成功安装  

check 64bit Environment
--------------------------------------
::

  [Environment]::Is64BitProcess

Display all environment variables
--------------------------------------
::

  gci env:
  ls env:
  
Some variables
::

  ProgramData                    C:\ProgramData
  ProgramFiles                   C:\Program Files
  ProgramFiles(x86)              C:\Program Files (x86)
  ProgramW6432                   C:\Program Files  

Display ProgramFiles variable
--------------------------------------
::
  
  $env:ProgramFiles  
  
Display ProgramFiles(x86) variable
--------------------------------------
::

  ${Env:ProgramFiles(x86)}
  or
  [Environment]::GetEnvironmentVariable("ProgramFiles(x86)")

show path
--------------------------------------
::

    $env:path

Remove
--------------------------------------
::

  Remove-Item * -Recurse -Force
  
or

::

  rm * -Recurse -Force  


Start testprj.sln from powershell
--------------------------------------
::

  & ./testprj.sln
  
Add Environment variables
--------------------------------------
::

  $env:path = "c:\dev\mingw64\bin\;"+$env:path

  
Find program
--------------------------------
::

  Get-Command gcc
  Results
  PS C:\Users\eric>  Get-Command gcc
  
  CommandType     Name           Version    Source
  -----------     ----           -------    ------
  Application     gcc.exe        0.0.0.0    c:\dev\mingw64\bin\gcc.exe
  
Find powershell
-------------------------------- 
::
 
  PS C:\Users\eric> Get-Command powershell
  
  CommandType     Name           Version    Source
  -----------     ----           -------    ------
  Application     powershell.exe 10.0.22... C:\Windows\System32\WindowsPowerShell\...
  
Find PowerShell 7 (x64)
-------------------------------- 
::
 
  Get-ChildItem -Path c:\ -Recurse pwsh.exe
  results:
  PS C:\Users\eric> Get-ChildItem -Path c:\ -Recurse pwsh.exe
    目录: C:\Program Files\PowerShell\7
  Mode                 LastWriteTime         Length Name
  ----                 -------------         ------ ----
  -a----        2023/11/11      1:59         280000 pwsh.exe

Write-Output
-----------------------------------
::

  $array=@(1,2,3,5,6,7,8);
  Write-Output "array.Length=$($array.Length)"
  typical results:
  array.Length=7
  
Write-Output
-----------------------------------
::

  $pets = @{Cat = 'Frisky'; Dog = 'Spot'; Fish = 'Nimo'; Hamster = 'Whiskers'}
  Write-Output "`$pets=$pets"
  $pets
  Write-Output "`$pets.Cat=$($pets.Cat)"
  Write-Output "`$pets.Dog=$($pets.Dog)"
  Write-Output "`$pets.Fish=$($pets.Fish)"
  Write-Output "`$pets.Hamster=$($pets.Hamster)"
  $pets=System.Collections.Hashtable
  
  Name                           Value
  ----                           -----
  Fish                           Nimo
  Cat                            Frisky
  Hamster                        Whiskers
  Dog                            Spot
  $pets.Cat=Frisky
  $pets.Dog=Spot
  $pets.Fish=Nimo
  $pets.Hamster=Whiskers  


Get-Command -Name
--------------------------------
::
  
  Get-Command -Name "sort.exe"  
  CommandType     Name        Version    Source
  -----------     ----        -------    ------
  Application     sort.exe    10.0.22... C:\Windows\system32\sort.exe  
  
  
msiexec.exe
--------------------------------
::

   msiexec.exe /i D:\work\batch_work\ModernBatchFiles\codes\download\03\msmpisdk.msi
   
How to install MSI file to the custom directory using PowerShell?
----------------------------------------------------------------------   
#. `How to install MSI file to the custom directory using PowerShell? <https://www.youtube.com/watch?v=dnZcZCHdx0k/>`_

To install the MSI file to the custom directory using PowerShell, we can use the TARGETDIR, INSTALLDIR, INSTALLPATH, etc arguments for the custom path, depending upon the MSI file that it supports.
::

  msiexec /i "C:\temp\7z1900-x64.msi" INSTALLDIR="D:\ProgramFiles\7zip" /quiet
  
The above command can be run into PowerShell and cmd both but you can’t control the process that to wait until the installation finishes. To control the above command, we can use the Start-Process cmdlet in PowerShell.
::

  Start-Process -FilePath "C:\windows\system32\msiexec.exe" -ArgumentList "/i
  C:\temp\7z1900-x64.msi INSTALLDIR='D:\ProgramFiles\7zip' /quiet" -Wait