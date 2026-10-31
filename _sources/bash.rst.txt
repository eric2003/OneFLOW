Bash
==================================

Bash
---------------------------------
#. `Bash Scripting on Linux <https://www.youtube.com/watch?v=2733cRPudvI&list=PLT98CRl2KxKGj-VKtApD8-zCqSaN2mD4w>`_
#. `Beginner's Guide to the Bash Terminal <https://www.youtube.com/watch?v=oxuRxtrO2Ag>`_
#. `3小时Bash脚本开发教程 <https://www.bilibili.com/video/BV1T64y1D7Eq/>`_
#. `【伯乐大典】最实用的Bash脚本知识 <https://www.bilibili.com/video/BV1AT411Y7bq/>`_
#. `Bash Scripting Tutorial <https://ryanstutorials.net/bash-scripting-tutorial/>`_
#. `【油管】Bash Shell Scripting Tutorial | Shell Scripting Tutorial <https://www.bilibili.com/video/BV1Uh411t7Bj/>`_
#. `bash少为人知的技巧 <https://www.bilibili.com/video/BV1NB4y1u74B/>`_
#. `Shell Scripting Tutorial for Beginners <https://www.youtube.com/watch?v=cQepf9fY6cE&list=PLS1QulWo1RIYmaxcEqw5JhK3b-6rgdWO_>`_
#. `A Complete Guide To The Bash Environment Variables <https://www.shell-tips.com/bash/environment-variables/>`_


::

  git clone git@github.com:eric2003/ModernBash.git

window git bash

cd d drive
::

  cd /d
  cd /d/work/bash_work/ModernBash/codes/simple/01/

grep
::

  $ ls -l /usr/bin | grep bash
  -rwxr-xr-x 1 root   root     1396520  1月  7  2022 bash
  -rwxr-xr-x 1 root   root        6818  1月  7  2022 bashbug
  -rwxr-xr-x 1 root   root        4414 11月 11  2021 dh_bash-completion
  lrwxrwxrwx 1 root   root           4  8月 30  2022 rbash -> bash

echo example 1
::
  
  $ echo Hello World! > hello.txt  
  $ cat hello.txt 
  Hello World!

echo example 2
::
  
  $ echo Hello World! > hello.txt
  $ echo Hello World! > hello.txt
  $ echo Hello World 1! >> hello.txt
  $ echo Hello World 2! >> hello.txt 
  $ cat hello.txt 
  Hello World!
  Hello World 1!
  Hello World 2!
  $ wc -w hello.txt 
  8 hello.txt
  
list
::

  $ MY_LIST=(one two three four five)
  $ echo $MY_LIST
  $ echo ${MY_LIST}
  $ echo ${MY_LIST[@]}
  $ echo ${MY_LIST[0]}
  $ echo ${MY_LIST[1]}
  one
  one
  one two three four five
  one
  two
  
cat /etc/shells
::

  $ cat /etc/shells
  # /etc/shells: valid login shells
  /bin/sh
  /bin/bash
  /usr/bin/bash
  /bin/rbash
  /usr/bin/rbash
  /usr/bin/sh
  /bin/dash
  /usr/bin/dash
  
Set Command
-----------------
#. `The set Command in Linux <https://www.baeldung.com/linux/set-command/>`_
#. `Linux set Command Syntax <https://www.baeldung.com/linux/set-command/>`_

set Command Options
---------------------

The -e Option
`````````````````
When a query returns a non-zero status, the -e flag stops the script. It also detects errors in the currently executing script.
::

  #!/bin/bash 
  set -e 
  mkdir newfolder 
  cat filenotindirectory 
  echo 'The file is not in our directory!'
  ---------------------
  $ bash testprj.sh 
  cat: filenotindirectory: No such file or directory
  

The -C Option
`````````````````  
Using the -C flag ensures that we cannot overwrite an existing file with the same name:
::

  #!/bin/bash 
  touch myfile
  echo 'An existing file' > myfile
  set -C
  echo 'Editing an existing file' > myfile
  typical results:
  $ bash testprj.sh 
  testprj.sh: line 5: myfile: cannot overwrite existing file
  
The -f Option  
`````````````````  
As we know, we can easily search for files using wildcard characters such as \?, \*, or []. 
This method is similar to regex, where we try to find similar texts using patterns.
The Bash shell uses the wildcards we specify to generate patterns and match them with filenames.
This feature is called globbing.
Let’s try using globbing by searching for files with the .txt extension:
::

  #!/bin/bash 
  touch files.txt
  ls *.txt
  set -f
  ls *.txt
  typical results:
  $ bash testprj.sh 
  files.txt  README.txt
  ls: cannot access '*.txt': No such file or directory
  
The -x Option   
`````````````````  
We use the -x parameter when debugging our scripts to determine the output of individual commands.

To illustrate, let’s create a Bash script that displays a countdown from 3 to 0:
::

  #!/bin/bash
  set -x
  n=3
  while [ $n -gt 0 ]; do
      n=$[ $n-1 ]
      echo $n
      sleep 1
  done
  typical results:
  $ bash testprj.sh 
  + n=3
  + '[' 3 -gt 0 ']'
  + n=2
  + echo 2
  2
  + sleep 1
  + '[' 2 -gt 0 ']'
  + n=1
  + echo 1
  1
  + sleep 1
  + '[' 1 -gt 0 ']'
  + n=0
  + echo 0
  0
  + sleep 1
  + '[' 0 -gt 0 ']'
  
The -a Option
`````````````````  
We can export variables or functions with this flag, making them reusable in subshells or scripts. First, let’s define variables in our terminal:
::

  1.
  $ set -a 
  $ name='May' 
  $ age=22
  
  2.
  #!/bin/bash 
  echo $name $age
  
  3.
  $ bash testprj.sh 
  may 22
  
The -u Option
`````````````````  
We use this flag to ensure that Bash does not overlook the non-existent variables in our script. We can see that in normal circumstances, Bash ignores the unassigned variables and runs our script without any errors:
::

  1.
  #!/bin/bash
  set -u
  x='Bash scripting is fun'
  echo $x $y
  
  2.
  $ bash testprj.sh
  testprj.sh: line 4: y: unbound variable
  
The +[argument] Option  
````````````````````````  
Running the set command with the +[argument] option unsets the option’s functionality. In essence, it nullifies the effect of the -[argument] option.

We will see some examples of shell scripts run with the set  +[argument] commands and their outputs.

Firstly, we consider set +e:
::

  1.
  #!/bin/bash 
  set +e 
  mkdir newfolder 
  cat filenotindirectory 
  echo 'The file is not in our directory!'
  
  2.
  $ bash testprj.sh
  cat: filenotindirectory: No such file or directory
  The file is not in our directory!
  
Other Special Variables
------------------------------
::

  $0 - The name of the Bash script.
  $1 - $9 - The first 9 arguments to the Bash script. (As mentioned above.)
  $# - How many arguments were passed to the Bash script.
  $@ - All the arguments supplied to the Bash script.
  $? - The exit status of the most recently run process.
  $$ - The process ID of the current script.
  $USER - The username of the user running the script.
  $HOSTNAME - The hostname of the machine the script is running on.
  $SECONDS - The number of seconds since the script was started.
  $RANDOM - Returns a different random number each time is it referred to.
  $LINENO - Returns the current line number in the Bash script.  
  
tr command
-----------------------------------------
#. `tr command in Unix/Linux with examples <https://www.geeksforgeeks.org/tr-command-in-unix-linux-with-examples/>`_
#. `Linux tr命令 <https://www.runoob.com/linux/linux-comm-tr.html>`_ 



cat command
-----------------------------------------
#. `How to Use the Linux cat Command With Examples <https://www.hostinger.com/tutorials/linux-cat-command-tutorial-and-examples/>`_  

Using the cat Command to Create a File
````````````````````````````````````````````
Using the cat command you can quickly create a file and put text into it. To do that, use the > redirect operator to redirect the text in the file.
::

  cat > filename.txt
  
The file is created, and you can begin populating it with text. To add multiple lines of text just press Enter at the end of each line.  Once you’re done, hit CTRL+D to exit the file.
 

  
  
