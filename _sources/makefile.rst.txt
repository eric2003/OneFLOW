Makefile
==================================

#. `Makefile从入门到精通系列链接整理 <https://zhuanlan.zhihu.com/p/398677004/>`_
#. `Step-by-step Makefile tutorial <https://gist.github.com/francois-rozet/c8efb19f66fed666263641d4e40f8863/>`_
#. `makefile 常用函数notdir、wildcard、patsubst <https://cloud.tencent.com/developer/article/2055032/>`_
#. `Writing Makefiles for Modern Fortran <https://aoterodelaroza.github.io/devnotes/modern-fortran-makefiles/>`_


The Wildcard Function 
-------------------------

#. `The Function wildcard <https://www.gnu.org/software/make/manual/html_node/Wildcard-Function.html>`_
#. `Using Makefile Wildcards <https://earthly.dev/blog/using-makefile-wildcards/>`_
#. `跟我一起写Makefile 1.0 文档 <https://seisman.github.io/how-to-write-makefile/introduction.html>`_
#. `Makefile Tutorial By Example <https://makefiletutorial.com/>`_
#. `Makefile的相关资源（不定期更新） <https://zhuanlan.zhihu.com/p/396719316/>`_
#. `makefile详解 <https://www.cnblogs.com/paul-617/p/15501875.html>`_


One use of the wildcard function is to get a list of all the C source files in a directory, like this:
::

  $(wildcard *.c)
  
The Patsubst Function
-----------------------------

#. `Functions for String Substitution and Analysis <https://www.gnu.org/software/make/manual/html_node/Text-Functions.html>`_
#. `patsubst, patsubsti NMAKE functions <https://learn.microsoft.com/en-us/cpp/build/reference/nmake-function-patsubst?view=msvc-170/>`_

::

  $(patsubst pattern,replacement,text)
  
We can change the list of C source files into a list of object files by replacing the ‘.c’ suffix with ‘.o’ in the result, like this:
::

  $(patsubst %.c,%.o,$(wildcard *.c))
  
Thus, a makefile to compile all C source files in the directory and then link them together could be written as follows:
::

  objects := $(patsubst %.c,%.o,$(wildcard *.c))
  
  foo : $(objects)
          cc -o foo $(objects)  

Phony Targets
--------------------
#. `4.6 Phony Targets <https://www.gnu.org/software/make/manual/html_node/Phony-Targets.html>`_
#. `What is the purpose of .PHONY in a Makefile? <https://stackoverflow.com/questions/2145590/what-is-the-purpose-of-phony-in-a-makefile/>`_

Let's assume you have install target, which is a very common in makefiles. If you do not use .PHONY, and a file named install exists in the same directory as the Makefile, then make install will do nothing. This is because Make interprets the rule to mean "execute such-and-such recipe to create the file named install". Since the file is already there, and its dependencies didn't change, nothing will be done.

However if you make the install target PHONY, it will tell the make tool that the target is fictional, and that make should not expect it to create the actual file. Hence it will not check whether the install file exists, meaning: a) its behavior will not be altered if the file does exist and b) extra stat() will not be called.

Generally all targets in your Makefile which do not produce an output file with the same name as the target name should be PHONY. This typically includes all, install, clean, distclean, and so on.

::

  .PHONY: clean
  clean:
          rm *.o temp

Automatic Variables
------------------------
#. `10.5.3 Automatic Variables <https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html>`_

自动变量，最常用包括：

$@用于表示一个规则中的目标。当我们的一个规则中有多个目标时，$@所指的是其中任何造成命令被运行的目标。
$^则表示的是规则中的所有先择条件。
$<表示的是规则中的第一个先决条件。

$@  表示目標文件

$^  表示所有的依賴文件

$<  表示第一個依賴文件

$?  表示比目標還要新的依賴文件列表

Makefile有三个非常有用的变量。分别是$@，$^，$<代表的意义分别是：

$@--目标文件，$^--所有的依赖文件，$<--第一个依赖文件。

Pattern Rule Examples
-----------------------
Here are some examples of pattern rules actually predefined in make. First, the rule that compiles ‘.c’ files into ‘.o’ files:
::

  %.o : %.c
          $(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

defines a rule that can make any file x.o from x.c. The recipe uses the automatic variables ‘$@’ and ‘$<’ to substitute the names of the target file and the source file in each case where the rule applies (see Automatic Variables).		  

print $@
------------------
::

  echo \$$@

::

  @echo \$$@ = $@
  @echo \$$< = $<
  @echo \$$^ = $^
  
Fortran Makefile
------------------
Example 1
::

  FC = gfortran
  FFLAGS = -O3 -Wall
  SRC = $(wildcard *.f90)
  OBJ = $(SRC:.f90=.o)
  $(info SRC=$(SRC))
  $(info OBJ=$(OBJ))
  
  %.o : %.f90
  	$(FC) $(FFLAGS) -o $@ -c $<
  
  testprj: $(OBJ)
  	$(FC) $(FFLAGS) -o $@ $^
  	@echo "echo1= \$$@ = $@ "
  	@echo "echo2= \$$^ = $^ "
  	@echo "echo3= \$$< = $< "
  	$(info $$@ = $@)
  	$(info $$^ = $^)
  	$(info $$< = $<)
  
  
  .PHONY: clean
  clean :
  	rm *.o testprj

windows make
::

  mingw32-make.exe
