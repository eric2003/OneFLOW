# tetgen.tcl - Tetgen import/export

set Import(tetfile) ""

proc tetgen_import {w name exe} {
  global ProgData Font Import

  set cmd [get_executable $exe 1]
  if {$cmd == ""} return
  set Import(cgnsfile) $ProgData(file,name)

  toplevel $w
  wm title $w "Tetgen File Import"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Import(done) 0}

  import_input $w tetfile Tetgen {.poly .smesh .node}
  import_output $w

  if {[import_buttons $w tetgen_import_check]} {
    lappend cmd $Import(tetfile) $Import(cgnsfile)
    import_run "Tetgen Import" $cmd $Import(cgnsfile)
  }
}

proc tetgen_import_check {w} {
  global Import
  if {[string trim $Import(tetfile)] == "" ||
      [string trim $Import(cgnsfile)] == ""} {
    errormsg "must specify a Tetgen and a CGNS file" $w
    return
  }
  if {![file exists $Import(tetfile)]} {
    errormsg "Tetgen input file doesn't exist" $w
    return
  }
  set Import(done) 1
}

