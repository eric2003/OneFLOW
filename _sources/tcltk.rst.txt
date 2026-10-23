Tcl/tk
==================================

Tcl/tk
---------------------------------
#. `Tcl Developer Xchange <https://www.tcl.tk/>`_
#. `Tcl/Tk 系列链接整理 <https://zhuanlan.zhihu.com/p/512365807/>`_
#. `Windows11+CGNS4.3.0+HDF51.12.2+源码编译CGNSTOOLS读取HDF5格式文件的问题 <https://zhuanlan.zhihu.com/p/517337277/>`_

original cgnsview.tcl
::

  proc file_stats {} {
    global ProgData
    set ProgData(file,size) "[file size $ProgData(file,name)] bytes"
    set ProgData(file,type) [CGNSfile $ProgData(file)]
    set ProgData(file,vers) [CGIOversion]

modified cgnsview.tcl
::

  proc file_stats {} {
    global ProgData
    set ProgData(file,size) "[file size $ProgData(file,name)] bytes"
    if {[catch {CGIOpen $ProgData(file)} msg]} {
      return $msg
    }
    set ProgData(file,type) [CGNSfile $ProgData(file)]
    set ProgData(file,vers) [CGIOversion]