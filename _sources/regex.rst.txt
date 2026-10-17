Regular Expressions
==================================

#. `Regular Expressions (RegEx) Tutorial #1 - What is RegEx? <https://www.youtube.com/watch?v=r6I-Ahc0HB4&list=PL4cUxeGkcC9g6m_6Sld9Q4jzqdqHd2HiD/>`_
#. `C++ Regex Tutorial: Regular Expressions In C++ With Examples <https://www.softwaretestinghelp.com/regex-in-cpp/>`_
#. `常用正则表达式示例-CodeSheep <https://r2coding.com/#/README?id=%E6%AD%A3%E5%88%99%E8%A1%A8%E8%BE%BE%E5%BC%8F/>`_
#. `Regular expressions (C++) <https://learn.microsoft.com/en-us/cpp/standard-library/regular-expressions-cpp?view=msvc-170/>`_
#. `C++与正则表达式 <https://paul.pub/cpp-regex/>`_
#. `C++ Tutorial 19 : C++ Regular Expressions <https://www.youtube.com/watch?v=9K4N6MO_R1Y>`_
#. `C++ Tutorial 21 : C++ Regex 3 <https://www.youtube.com/watch?v=tsR2rNTX3D0/>`_
#. `6.2. re — Regular expression operations <https://python.readthedocs.io/en/latest/library/re.html>`_
#. `第22课【C++ 正则表达式】正则表达式基本单元，C++正则替换与捕获，正则替换，enum_class <https://www.bilibili.com/video/BV1Rp4y1j7D9/>`_
#. `C++ Regex 101 <https://www.fluentcpp.com/2020/02/28/c-regex-101-simple-code-for-simple-cases-with-regexes/>`_
#. `【韩顺平讲Java】Java 正则表达式专题 -正则 正则表达式 元字符 限定符 Pattern Matcher 分组 捕获 反向引用等 <https://www.bilibili.com/video/BV1Eq4y1E79W/>`_

Python
------------------------------------------
#. `《孙兴华讲正则表达式》基于Python【已完结】 <https://www.bilibili.com/video/BV1kp4y1C7c8/>`_


Online Regular Expressions
-------------------------------
#. `RegEx101.com <https://RegEx101.com/>`_

::

  git clone git@github.com:eric2003/ModernRegex.git
  
.gitignore
::

  **/build/*  
  
HDF5 CMake Info
::

  -- CONFIG_DATE = 2023-12-18
  -- HDF5_USE_FOLDERS = ON
  -- HDF_CONFIG_DIR = D:/work/hdf5_work/hdf5-1.14.3/config
  -- HDF_RESOURCES_DIR = D:/work/hdf5_work/hdf5-1.14.3/config/cmake
  -- HDF5_SOURCE_DIR = D:/work/hdf5_work/hdf5-1.14.3
  -- HDF5_SRC_DIR = D:/work/hdf5_work/hdf5-1.14.3/src
  -- HDF5_TEST_SRC_DIR = D:/work/hdf5_work/hdf5-1.14.3/test
  -- HDF5_TEST_PAR_DIR = D:/work/hdf5_work/hdf5-1.14.3/testpar
  -- HDF5_TEST_API_SRC_DIR = D:/work/hdf5_work/hdf5-1.14.3/test/API
  -- HDF5_TEST_API_PAR_SRC_DIR = D:/work/hdf5_work/hdf5-1.14.3/testpar/API
  -- HDF_RESOURCES_DIR = D:/work/hdf5_work/hdf5-1.14.3/config/cmake
  -- HDF5_PACKAGE_NAME = HDF5
  
.. list-table::
   :header-rows: 1

   * - Expression
     - Meaning
   * - \`c\`
     - matches any literal character \`c\`
   * - \`\\\\d\`
     - matches any decimal digit
   * - \`\\\\D\`
     - matches any character that's not a decimal digit
   * - \`\\\\f\`
     - matches \`\\f\`
   * - \`\\\\n\`
     - matches \`\\n\`
   * - \`\\\\r\`
     - matches \`\\r\`
   * - \`\\\\s\`
     - matches any ASCII whitespace, including \`\\n\`
   * - \`\\\\S\`
     - matches any character that's not a whitespace
   * - \`\\\\t\`
     - matches \`\\t\`
   * - \`\\\\v\`
     - matches \`\\v\`
   * - \`\\\\w\`
     - matches any letter, \`_\`, or decimal digit
   * - \`\\\\W\`
     - matches any character that \`\\\\w\` doesn't match
   * - \`\\\\c\`
     - matches any literal character \`c\`, which must be a punctuation
   * - \`.\`
     - matches any single character except \`\\n\`
   * - \`A?\`
     - matches 0 or 1 occurrences of \`A\`
   * - \`A*\`
     - matches 0 or many occurrences of \`A\`
   * - \`A+\`
     - matches 1 or many occurrences of \`A\`
   * - \`^\`
     - matches the beginning of a string (not that of each line)
   * - \`$\`
     - matches the end of a string (not that of each line)
   * - \`xy\`
     - matches \`x\` followed by \`y\`

基础正则表达式速查表

字符
------------

.. list-table::
   :header-rows: 1

   * - 表达式
     - 描述
   * - [abc]
     - 字符集。匹配集合中所含的任一字符。
   * - [^abc]
     - 否定字符集。匹配任何不在集合中的字符
   * - [a-z]
     - 字符范围。匹配指定范围内的任意字符。
   * - .
     - 匹配除换行符以外的任何单个字符。
   * - \w
     - 匹配任何字母数字，包括下划线（等价于[A-Za-z0-9\_]）。
   * - \W
     - 匹配任何非字母数字（等价于[^A-Za-z0-9\_]）。
   * - \d
     - 数字。匹配任何数字。
   * - \D
     - 非数字。匹配任何非数字字符。
   * - \s
     - 空白。匹配任何空白字符，包括空格、制表符等。
   * - \S
     - 非空白。匹配任何非空白字符。

分组和引用
------------

.. list-table::
   :header-rows: 1

   * - 表达式
     - 描述
   * - (expression)
     - 分组。匹配括号里的整个表达式。
   * - (?:expression)
     - 非捕获分组。匹配括号里的整个字符串但不获取匹配结果，拿不到分组引用。
   * - \\num
     - 对前面所匹配分组的引用。比如(\\d)\\1可以匹配两个相同的数字，(ABC)(xyz)\\1\\2则可以匹配ABCxyzABCxyz。

锚点/边界
------------

.. list-table::
   :header-rows: 1

   * - 表达式
     - 描述
   * - ^
     - 匹配字符串或行开头。
   * - $
     - 匹配字符串或行结尾。
   * - \\b
     - 匹配单词边界。比如World\\b可以匹配HelloWorld末尾的World，不能匹配HelloWorldHello中的World。
   * - \\B
     - 匹配非单词边界。比如World\\B可以匹配HelloWorldSheep中的World，不能匹配HelloWorld中的World。


数量表示
------------

.. list-table::
   :header-rows: 1

   * - 表达式
     - 描述
   * - ?
     - 匹配前面的表达式0个或1个。即表示可选项。
   * - \+
     - 匹配前面的表达式至少1个。
   * - \*
     - 匹配前面的表达式0个或多个。
   * - \|
     - 或运算符。并集，可以匹配符号前后的表达式。
   * - {m}
     - 匹配前面的表达式m个。
   * - {m,}
     - 匹配前面的表达式最少m个。
   * - {m,n}
     - 匹配前面的表达式最少m个，最多n个。

预查断言
------------

.. list-table::
   :header-rows: 1

   * - 表达式
     - 描述
   * - (?=)
     - 正向预查。比如Hello(?=World)能匹配HelloWorld中的Hello，但不能匹配HelloKitty中的Hello。
   * - (?!)
     - 正向否定预查。比如Hello(?!World)不能匹配HelloWorld中的Hello，但能匹配HelloKitty中的Hello。
   * - (?<=)
     - 反向预查。比如(?<=Hello)World能匹配HelloWorld中的World，但不能匹配HiWorld中的World。
   * - (?<!)
     - 反向否定预查。比如(?<!Hello)World不能匹配HelloWorld中的World，但能匹配HiWorld中的World。

特殊标志
------------

.. list-table::
   :header-rows: 1

   * - 表达式
     - 描述
   * - /.../i
     - 忽略大小写。
   * - /.../g
     - 全局匹配。
   * - /.../m
     - 多行修饰符。用于多行匹配。

特殊字符
------------
#. `特殊字符 <https://paul.pub/cpp-regex/>`_

.. list-table::
   :header-rows: 1

   * - 表达式
     - 描述
   * - .
     - 匹配任意字符
   * - [
     - 字符类的开始
   * - ]
     - 字符类的结束
   * - {
     - 量词重复数开始
   * - }
     - 量词重复数结束
   * - (
     - 分组开始
   * - )
     - 分组结束
   * - \\
     - 转义字符
   * - \\ \\
     - 转义字符自身
   * - \*
     - 量词，0个或者多个
   * - \+
     - 量词，1个或者多个
   * - ?
     - 量词，0个或者1个 
   * - \|
     - 或
   * - ^
     - 行开始；否定
   * - $
     - 行结束
   * - \\n
     - 换行
   * - \\t
     - Tab符
   * - \\xhh
     - hh表示两位十六进展表示的Unicode字符
   * - \\xhhhh
     - hhhh表示四位十六进制表示的Unicode字符

字符类
------------
.. list-table::
   :header-rows: 1

   * - 字符类
     - 简写
     - 说明
   * - [[:alnum:]]
     - 
     - 字母和数字
   * - [_[:alnum:]]
     - \、w
     - 字母，数字以及下划线
   * - [^_[:alnum:]]
     - \\W
     - 非字母，数字以及下划线
   * - [[:digit:]]
     - \\d
     - 数字
   * - [^[:digit:]]
     - \\D
     - 非数字
   * - [[:space:]]
     - \\s
     - 空白字符
   * - [^[:space:]]
     - \\S
     - 非空白字符
   * - [^[:space:]]
     - \\S
     - 非空白字符
   * - [[:lower:]]
     - 
     - 小写字母
   * - [[:upper:]]
     - 
     - 大写字母
   * - [[:alpha:]]
     - 
     - 任意字母
   * - [[:blank:]]
     - 
     - 非换行符的空白字符
   * - [[:cntrl:]]
     - 
     - 控制字符
   * - [[:graph:]]
     - 
     - 图形字符
   * - [[:print:]]
     - 
     - 可打印字符
   * - [[:punct:]]
     - 
     - 标点字符
   * - [[:xdigit:]]
     - 
     - 十六进制的数字字符

重复
----------
.. list-table::
   :header-rows: 1

   * - 表达式
     - 描述
   * - {n}
     - 重复n次
   * - {n,}
     - 重复n或更多次
   * - {n,m}
     - 重复[n ~ m]次
   * - \*
     - 重复0次或多次，等同于{0,}
   * - \+
     - 重复1次或多次，等同于{1,}
   * - ?
     - 重复0次或1次，等同于{0,1}

Example
----------
.. list-table::
   :header-rows: 1

   * - 表达式
     - 描述
   * - .*
     - 匹配任意字符,重复0次或多次
   * - {n,}
     - 重复n或更多次
   * - {n,m}
     - 重复[n ~ m]次
   * - \*
     - 重复0次或多次，等同于{0,}
   * - \+
     - 重复1次或多次，等同于{1,}
   * - ?
     - 重复0次或1次，等同于{0,1}	 
