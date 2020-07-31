# fast.tcl - FAST file import/export

array set Import {
  fastfile ""
  fast,fmt s
  fast,endian b
  fast,prec 8
}

array set Export {
  fastfile ""
  fast,fmt s
  fast,endian b
  fast,prec 8
  fast,sym ""
}

proc fast_import {w name exe} {
  global ProgData Import Font

  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  toplevel $w
  wm title $w "FAST File Import"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Import(done) 0}

  FrameCreate $w.input -text "FAST Input" -font $Font(bold)
  pack $w.input -side top -padx 5 -pady 2 -fill x
  set f [FrameGet $w.input]

  label $f.lab -text Filename -width 8 -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Import(fastfile) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "fast_browse $w"
  pack $f.but -side right -fill y

  import_output $w 0 8

  FrameCreate $w.fmt -text "File Format" -font $Font(bold)
  pack $w.fmt -side top -padx 5 -pady 2 -fill x
  set fmt [FrameGet $w.fmt]

  set f [frame $fmt.fmt]
  pack $f -side top -fill x

  radiobutton $f.b -text binary \
    -variable Import(fast,fmt) -value s
  radiobutton $f.u -text unformatted \
    -variable Import(fast,fmt) -value u
  radiobutton $f.f -text formatted \
    -variable Import(fast,fmt) -value f
  pack $f.b $f.u $f.f -side left -expand 1

  frame $fmt.sep -height 2 -relief sunken -bd 1
  pack $fmt.sep -side top -fill x -pady 3

  set f [frame $fmt.end]
  pack $f -side top -fill x

  radiobutton $f.b -text "big-endian" \
    -variable Import(fast,endian) -value b
  radiobutton $f.l -text "little-endian" \
    -variable Import(fast,endian) -value l
  checkbutton $f.d -text "real*8" \
    -variable Import(fast,prec) -onvalue 8 -offvalue 4
  pack $f.b $f.l $f.d -side left -expand 1

  if {[import_buttons $w fast_import_check]} {
    if {$Import(fast,fmt) != "s"} {
      lappend cmd -$Import(fast,fmt)
    }
    if {$Import(fast,endian) != "b"} {
      lappend cmd -$Import(fast,endian)
    }
    if {$Import(fast,prec) != "8"} {
      lappend cmd -$Import(fast,prec)
    }
    lappend cmd $Import(fastfile) $Import(cgnsfile)
    import_run "FAST Import" $cmd $Import(cgnsfile)
  }
}

proc fast_import_check {w} {
  global Import
  if {[string trim $Import(fastfile)] == "" ||
      [string trim $Import(cgnsfile)] == ""} {
    errormsg "must specify a FAST input and a CGNS output file" $w
    return
  }
  if {![file exists $Import(fastfile)]} {
    errormsg "FAST input file doesn't exist" $w
    return
  }
  set Import(done) 1
}

proc fast_browse {w} {
  global Import tcl_platform
  set fname [FileOpen "FAST File" $Import(fastfile) $w \
    {{{FAST Files} {.fgrid}} {{All Files} {*}}}]
  if {$fname != ""} {
    if {$tcl_platform(platform) == "windows"} {
      set Import(fastfile) [join [split $fname /] \\]
    } else {
      set Import(fastfile) $fname
    }
    set Import(cgnsfile) [file rootname $Import(fastfile)].cgns
    if {[string first ".r8" $fname] > 0} {
      set Import(fast,fmt) u
      set Import(fast,prec) 8
      set Import(fast,endian) b
    } elseif {[string first ".r4" $fname] > 0} {
      set Import(fast,fmt) u
      set Import(fast,prec) 4
      set Import(fast,endian) b
    } elseif {[string first ".b8" $fname] > 0} {
      set Import(fast,fmt) s
      set Import(fast,prec) 8
      set Import(fast,endian) b
    } elseif {[string first ".b4" $fname] > 0} {
      set Import(fast,fmt) s
      set Import(fast,prec) 4
      set Import(fast,endian) b
    } elseif {[string first ".lr8" $fname] > 0} {
      set Import(fast,fmt) u
      set Import(fast,prec) 8
      set Import(fast,endian) l
    } elseif {[string first ".lr4" $fname] > 0} {
      set Import(fast,fmt) u
      set Import(fast,prec) 4
      set Import(fast,endian) l
    } elseif {[string first ".lb8" $fname] > 0} {
      set Import(fast,fmt) s
      set Import(fast,prec) 8
      set Import(fast,endian) l
    } elseif {[string first ".lb4" $fname] > 0} {
      set Import(fast,fmt) s
      set Import(fast,prec) 4
      set Import(fast,endian) l
    } else {
      set Import(fast,fmt) f
      set Import(fast,prec) 8
      set Import(fast,endian) b
    }
  }
}

proc fast_export {w name exe} {
  global ProgData Font Export
  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  toplevel $w
  wm title $w "FAST Export"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Export(done) 0}

  export_input $w 1 1
  export_output $w fastfile FAST {.fgrid}

  FrameCreate $w.fmt -text "File Format" -font $Font(bold)
  pack $w.fmt -side top -padx 5 -pady 2 -fill x
  set fmt [FrameGet $w.fmt]

  set f [frame $fmt.fmt]
  pack $f -side top -fill x

  radiobutton $f.b -text binary -command fast_extension \
    -variable Export(fast,fmt) -value s
  radiobutton $f.f -text formatted -command fast_extension \
    -variable Export(fast,fmt) -value f
  pack $f.b $f.f -side left -expand 1

  frame $fmt.sep -height 2 -relief sunken -bd 1
  pack $fmt.sep -side top -fill x -pady 3

  set f [frame $fmt.end]
  pack $f -side top -fill x

  radiobutton $f.b -text "big-endian" -command fast_extension \
    -variable Export(fast,endian) -value b
  radiobutton $f.l -text "little-endian" -command fast_extension \
    -variable Export(fast,endian) -value l
  checkbutton $f.d -text "real*8" -command fast_extension \
    -variable Export(fast,prec) -onvalue 8 -offvalue 4
  pack $f.b $f.l $f.d -side left -expand 1

  FrameCreate $w.sym -text "Symmetry" -font $Font(bold)
  pack $w.sym -side top -padx 5 -pady 2 -fill x
  set sym [FrameGet $w.sym]

  radiobutton $sym.none -text none -variable Export(fast,sym) -value ""
  radiobutton $sym.x -text X -variable Export(fast,sym) -value x
  radiobutton $sym.y -text Y -variable Export(fast,sym) -value y
  radiobutton $sym.z -text Z -variable Export(fast,sym) -value z
  pack $sym.none $sym.x $sym.y $sym.z -side left -expand 1

  set Export(cgnsfile) $ProgData(file,name)
  set Export(fastfile) [file rootname $ProgData(file,name)].ugrid
  fast_extension

  if {[export_buttons $w fast_export_check]} {
    if {$Export(basenum) != ""} {
      lappend cmd -B$Export(basenum)
    }
    if {$Export(zonenum) != ""} {
      lappend cmd -Z$Export(zonenum)
    }
    if {$Export(fast,fmt) != "s"} {
      lappend cmd -$Export(fast,fmt)
    }
    if {$Export(fast,sym) != ""} {
      lappend cmd -$Export(fast,sym)
    }
    lappend cmd -$Export(fast,endian) -$Export(fast,prec)
    lappend cmd $Export(cgnsfile) $Export(fastfile)
    update
    run_command "FAST Export" $cmd
  }
}

proc fast_extension {} {
  global Export
  if {$Export(fastfile) == ""} return
  set basename [file rootname $Export(fastfile)]
  set ext [file extension $basename]
  if {$ext == ".b4" || $ext == ".b8" || $ext == ".lb4" || $ext == ".lb8"} {
    set basename [file rootname $basename]
  }
  if {$Export(fast,fmt) == "s"} {
    append basename "."
    if {$Export(fast,endian) == "l"} {append basename l}
    append basename "b$Export(fast,prec)"
  }
  set Export(fastfile) "$basename.fgrid"
}

proc fast_export_check {w} {
  global Export
  if {[string trim $Export(cgnsfile)] == "" ||
      [string trim $Export(fastfile)] == ""} {
    errormsg "must specify a CGNS and an FAST file" $w
    return
  }
  if {![file exists $Export(cgnsfile)]} {
    errormsg "CGNS input file doesn't exist" $w
    return
  }
  set Export(done) 1
}

