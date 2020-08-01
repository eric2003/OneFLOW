# aflr3.tcl - AFLR3 file import/export

array set Import {
  aflr3file ""
  aflr3,fmt s
  aflr3,endian b
  aflr3,prec 8
}

array set Export {
  aflr3file ""
  aflr3,fmt s
  aflr3,endian b
  aflr3,prec 8
  aflr3,sym ""
}

proc aflr3_import {w name exe} {
  global ProgData Import Font

  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  toplevel $w
  wm title $w "AFLR3 File Import"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Import(done) 0}

  FrameCreate $w.input -text "AFLR3 Input" -font $Font(bold)
  pack $w.input -side top -padx 5 -pady 2 -fill x
  set f [FrameGet $w.input]

  label $f.lab -text Filename -width 8 -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Import(aflr3file) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "aflr3_browse $w"
  pack $f.but -side right -fill y

  import_output $w 0 8

  FrameCreate $w.fmt -text "File Format" -font $Font(bold)
  pack $w.fmt -side top -padx 5 -pady 2 -fill x
  set fmt [FrameGet $w.fmt]

  set f [frame $fmt.fmt]
  pack $f -side top -fill x

  radiobutton $f.b -text binary \
    -variable Import(aflr3,fmt) -value s
  radiobutton $f.u -text unformatted \
    -variable Import(aflr3,fmt) -value u
  radiobutton $f.f -text formatted \
    -variable Import(aflr3,fmt) -value f
  pack $f.b $f.u $f.f -side left -expand 1

  frame $fmt.sep -height 2 -relief sunken -bd 1
  pack $fmt.sep -side top -fill x -pady 3

  set f [frame $fmt.end]
  pack $f -side top -fill x

  radiobutton $f.b -text "big-endian" \
    -variable Import(aflr3,endian) -value b
  radiobutton $f.l -text "little-endian" \
    -variable Import(aflr3,endian) -value l
  checkbutton $f.d -text "real*8" \
    -variable Import(aflr3,prec) -onvalue 8 -offvalue 4
  pack $f.b $f.l $f.d -side left -expand 1

  if {[import_buttons $w aflr3_import_check]} {
    if {$Import(aflr3,fmt) != "s"} {
      lappend cmd -$Import(aflr3,fmt)
    }
    if {$Import(aflr3,endian) != "b"} {
      lappend cmd -$Import(aflr3,endian)
    }
    if {$Import(aflr3,prec) != "8"} {
      lappend cmd -$Import(aflr3,prec)
    }
    lappend cmd $Import(aflr3file) $Import(cgnsfile)
    import_run "AFLR3 Import" $cmd $Import(cgnsfile)
  }
}

proc aflr3_import_check {w} {
  global Import
  if {[string trim $Import(aflr3file)] == "" ||
      [string trim $Import(cgnsfile)] == ""} {
    errormsg "must specify a AFLR3 input and a CGNS output file" $w
    return
  }
  if {![file exists $Import(aflr3file)]} {
    errormsg "AFLR3 input file doesn't exist" $w
    return
  }
  set Import(done) 1
}

proc aflr3_browse {w} {
  global Import tcl_platform
  set fname [FileOpen "AFLR3 File" $Import(aflr3file) $w \
    {{{AFLR3 Files} {.ugrid}} {{All Files} {*}}}]
  if {$fname != ""} {
    if {$tcl_platform(platform) == "windows"} {
      set Import(aflr3file) [join [split $fname /] \\]
    } else {
      set Import(aflr3file) $fname
    }
    set Import(cgnsfile) [file rootname $Import(aflr3file)].cgns
    if {[string first ".r8" $fname] > 0} {
      set Import(aflr3,fmt) u
      set Import(aflr3,prec) 8
      set Import(aflr3,endian) b
    } elseif {[string first ".r4" $fname] > 0} {
      set Import(aflr3,fmt) u
      set Import(aflr3,prec) 4
      set Import(aflr3,endian) b
    } elseif {[string first ".b8" $fname] > 0} {
      set Import(aflr3,fmt) s
      set Import(aflr3,prec) 8
      set Import(aflr3,endian) b
    } elseif {[string first ".b4" $fname] > 0} {
      set Import(aflr3,fmt) s
      set Import(aflr3,prec) 4
      set Import(aflr3,endian) b
    } elseif {[string first ".lr8" $fname] > 0} {
      set Import(aflr3,fmt) u
      set Import(aflr3,prec) 8
      set Import(aflr3,endian) l
    } elseif {[string first ".lr4" $fname] > 0} {
      set Import(aflr3,fmt) u
      set Import(aflr3,prec) 4
      set Import(aflr3,endian) l
    } elseif {[string first ".lb8" $fname] > 0} {
      set Import(aflr3,fmt) s
      set Import(aflr3,prec) 8
      set Import(aflr3,endian) l
    } elseif {[string first ".lb4" $fname] > 0} {
      set Import(aflr3,fmt) s
      set Import(aflr3,prec) 4
      set Import(aflr3,endian) l
    } else {
      set Import(aflr3,fmt) f
      set Import(aflr3,prec) 8
      set Import(aflr3,endian) b
    }
  }
}

proc aflr3_export {w name exe} {
  global ProgData Font Export
  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  toplevel $w
  wm title $w "AFLR3 Export"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Export(done) 0}

  export_input $w 1 1
  export_output $w aflr3file AFLR3 {.ugrid}

  FrameCreate $w.fmt -text "File Format" -font $Font(bold)
  pack $w.fmt -side top -padx 5 -pady 2 -fill x
  set fmt [FrameGet $w.fmt]

  set f [frame $fmt.fmt]
  pack $f -side top -fill x

  radiobutton $f.b -text binary -command aflr3_extension \
    -variable Export(aflr3,fmt) -value s
  radiobutton $f.f -text formatted -command aflr3_extension \
    -variable Export(aflr3,fmt) -value f
  pack $f.b $f.f -side left -expand 1

  frame $fmt.sep -height 2 -relief sunken -bd 1
  pack $fmt.sep -side top -fill x -pady 3

  set f [frame $fmt.end]
  pack $f -side top -fill x

  radiobutton $f.b -text "big-endian" -command aflr3_extension \
    -variable Export(aflr3,endian) -value b
  radiobutton $f.l -text "little-endian" -command aflr3_extension \
    -variable Export(aflr3,endian) -value l
  checkbutton $f.d -text "real*8" -command aflr3_extension \
    -variable Export(aflr3,prec) -onvalue 8 -offvalue 4
  pack $f.b $f.l $f.d -side left -expand 1

  FrameCreate $w.sym -text "Symmetry" -font $Font(bold)
  pack $w.sym -side top -padx 5 -pady 2 -fill x
  set sym [FrameGet $w.sym]

  radiobutton $sym.none -text none -variable Export(aflr3,sym) -value ""
  radiobutton $sym.x -text X -variable Export(aflr3,sym) -value x
  radiobutton $sym.y -text Y -variable Export(aflr3,sym) -value y
  radiobutton $sym.z -text Z -variable Export(aflr3,sym) -value z
  pack $sym.none $sym.x $sym.y $sym.z -side left -expand 1

  set Export(cgnsfile) $ProgData(file,name)
  set Export(aflr3file) [file rootname $ProgData(file,name)].ugrid
  aflr3_extension

  if {[export_buttons $w aflr3_export_check]} {
    if {$Export(basenum) != ""} {
      lappend cmd -B$Export(basenum)
    }
    if {$Export(zonenum) != ""} {
      lappend cmd -Z$Export(zonenum)
    }
    if {$Export(aflr3,fmt) != "s"} {
      lappend cmd -$Export(aflr3,fmt)
    }
    if {$Export(aflr3,sym) != ""} {
      lappend cmd -$Export(aflr3,sym)
    }
    lappend cmd -$Export(aflr3,endian) -$Export(aflr3,prec)
    lappend cmd $Export(cgnsfile) $Export(aflr3file)
    update
    run_command "AFLR3 Export" $cmd
  }
}

proc aflr3_extension {} {
  global Export
  if {$Export(aflr3file) == ""} return
  set basename [file rootname $Export(aflr3file)]
  set ext [file extension $basename]
  if {$ext == ".b4" || $ext == ".b8" || $ext == ".lb4" || $ext == ".lb8"} {
    set basename [file rootname $basename]
  }
  if {$Export(aflr3,fmt) == "s"} {
    append basename "."
    if {$Export(aflr3,endian) == "l"} {append basename l}
    append basename "b$Export(aflr3,prec)"
  }
  set Export(aflr3file) "$basename.ugrid"
}

proc aflr3_export_check {w} {
  global Export
  if {[string trim $Export(cgnsfile)] == "" ||
      [string trim $Export(aflr3file)] == ""} {
    errormsg "must specify a CGNS and an AFLR3 file" $w
    return
  }
  if {![file exists $Export(cgnsfile)]} {
    errormsg "CGNS input file doesn't exist" $w
    return
  }
  set Export(done) 1
}

