Linux
==================================

#. `Linux Commands <https://www.geeksforgeeks.org/linux-commands/>`_

Ubuntu
----------------
::

  sudo apt-get update
  
::

  sudo apt-get upgrade
  
Linux Commands
---------------------
#. `Linux Commands <https://www.geeksforgeeks.org/linux-commands/>`_

`touch command in Linux with Examples <https://www.geeksforgeeks.org/touch-command-in-linux-with-examples/>`_

The touch command is a standard command used in the UNIX/Linux operating system which is used to create, change and modify the timestamps of a file. Basically, there are two different commands to create a file in the Linux system which are as follows:

cat command: 
`````````````````
It is used to create the file with content.

touch command: 
`````````````````
It is used to create a file without any content. The file created using the touch command is empty. This command can be used when the user doesn’t have data to store at the time of file creation.

::

  sudo vim /etc/apt/sources.list


uninstall cmake
::

  sudo apt-get remove cmake

install cmake
::

  cd /home/eric/Downloads
  sudo sh cmake-3.28.0-linux-x86_64.sh --prefix=/usr/local

check ubuntu ``version`` method 1
::

  lsb_release -a
  Typical results:
  eric@eric-virtual-machine:~$ lsb_release -a
  No LSB modules are available.
  Distributor ID:	Ubuntu
  Description:	Ubuntu 22.04.3 LTS
  Release:	22.04
  Codename:	jammy
  
check ubuntu ``version`` method 2
::
  
  cat /etc/lsb-release
  Typical results:
  eric@eric-virtual-machine:~$ cat /etc/lsb-release
  DISTRIB_ID=Ubuntu
  DISTRIB_RELEASE=22.04
  DISTRIB_CODENAME=jammy
  DISTRIB_DESCRIPTION="Ubuntu 22.04.3 LTS"
  
check ubuntu ``version`` method 3
::

  $  uname -a
  Typical results:
  eric@eric-virtual-machine:~$  uname -a
  Linux eric-virtual-machine 6.2.0-37-generic #38~22.04.1-Ubuntu SMP PREEMPT_DYNAMIC Thu Nov  2 18:01:13 UTC 2 x86_64 x86_64 x86_64 GNU/Linux


Vim
`````````````````
#. `Linux vi/vim <https://www.runoob.com/linux/linux-vim.html>`_
