# vgrid.tcl - VGRID import/export

array set Import {
  vgridfile ""
  vgrid,combine 0
}

proc vgrid_import {w name exe} {
  global ProgData Font Import

  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  toplevel $w
  wm title $w "VGRID File Import"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Import(done) 0}

  import_input $w vgridfile VGRID {.cogsg .bc .mapbc}
  import_output $w 0

  FrameCreate $w.opts -text "Options" -font $Font(bold)
  pack $w.opts -side top -padx 5 -pady 2 -fill x
  set opts [FrameGet $w.opts]

  set f [frame $opts.f]
  pack $f -side left

  checkbutton $f.cb -text "Combine Patches where possible" \
    -variable Import(vgrid,combine) -onvalue 1 -offvalue 0
  pack $f.cb -side top -anchor w

  if {[import_buttons $w vgrid_import_check]} {
    if {$Import(vgrid,combine)} {
      lappend cmd -c
    }
    lappend cmd $Import(vgridfile) $Import(cgnsfile)
    import_run "VGRID Import" $cmd $Import(cgnsfile)
  }
}

proc vgrid_import_check {w} {
  global Import
  if {[string trim $Import(vgridfile)] == "" ||
      [string trim $Import(cgnsfile)] == ""} {
    errormsg "must specify a VGRID and a CGNS file" $w
    return
  }
  if {![file exists $Import(vgridfile)]} {
    errormsg "VGRID input file doesn't exist" $w
    return
  }
  set Import(done) 1
}

