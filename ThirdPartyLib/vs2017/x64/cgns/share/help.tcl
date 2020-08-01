
array set HelpData {
  browser ""
  htmlopts ""
  cgns ""
  winhtml 0
  chmfile ""
  cgns,http "http://www.grc.nasa.gov/www/cgns/CGNS_docs_beta/index.html"
}

proc get_browser {} {
  global tcl_platform env

  if {$tcl_platform(platform) != "windows"} {
    foreach name {firefox netscape mozilla} {
      set browser [find_file executable $name]
      if {$browser != ""} {return $browser}
    }
    return ""
  }

  # get browser from registry

  if {![catch {registry get HKEY_CLASSES_ROOT\\.htm {}} htmlfile] ||
      ![catch {registry get HKEY_CLASSES_ROOT\\.html {}} htmlfile]} {
    if {![catch {registry get \
      HKEY_CLASSES_ROOT\\$htmlfile\\shell\\open\\command {}} browser]} {
      return [join [split $browser \\] /]
    }
  }

  # use assoc and ftype commands

  if {[info exists env(COMSPEC)]} {
    set comspec $env(COMSPEC)
  } elseif {[string match "*9?" $tcl_platform(os)]} {
    set comspec command
  } else {
    set comspec cmd
  }
  if {[catch {exec $comspec /c assoc .html} assoc]} {return ""}
  set ftype [string trim [lindex [split $assoc =] 1]]
  if {$ftype == ""} {return ""}
  if {[catch {exec $comspec /c ftype $ftype} cmdstr]} {return ""}
  set cmd [string trim [lindex [split $cmdstr =] 1]]
  return [join [split $cmd \\] /]
}

proc help_defaults {} {
  global HelpData HelpNew cmd_dir tcl_platform
  array set HelpNew {
    browser ""
    htmlopts ""
  }
  set browser [get_browser]
  if {$browser != "" && $tcl_platform(platform) == "windows"} {
    if {[catch {file attributes $browser -shortname} name]} {
      if {[catch {file attributes [lindex $browser 0] -shortname} name]} {
        set name $browser
      } else {
        for {set n 1} {$n < [llength $browser]} {incr n} {
          set opt [lindex $browser $n]
          if {$opt != "%1"} {
            lappend HelpNew(htmlopts) $opt
          }
        }
      }
    }
    set browser $name
  }
  set HelpNew(browser) $browser
  set HelpNew(cgns) $HelpData(cgns,http)
}

proc help_init {args} {
  global HelpData HelpNew
  help_defaults
  foreach i {browser htmlopts cgns} {
    set HelpData($i) $HelpNew($i)
  }
  if {[info procs tclreg_get] != ""} {
    foreach i {browser htmlopts cgns} {
      if {![catch {tclreg_get Help $i} val]} {
        set HelpData($i) $val
      }
    }
  }
  if {[info commands WinHtml] != ""} {
    set HelpData(winhtml) 1
  }
  catch help_menu
}

proc help_setup {} {
  global HelpData HelpNew tcl_platform Font

  set w .hlpsetup
  catch {destroy $w}
  toplevel $w
  wm title $w "Help Setup"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW "destroy $w"

  set lw 7
  set ew 50

  foreach i {browser htmlopts cgns} {
    set HelpNew($i) $HelpData($i)
  }

  FrameCreate $w.hb -text "HTML Browser" -font $Font(bold)
  pack $w.hb -side top -padx 5 -pady 5
  set hb [FrameGet $w.hb]

  foreach j {{browser Path} {htmlopts Options}} {
    set i [lindex $j 0]
    set f [frame $hb.$i]
    pack $f -side top -fill x -expand 1
    label $f.lab -text [lindex $j 1] -width $lw -anchor w
    pack $f.lab -side left
    entry $f.ent -width $ew -textvariable HelpNew($i) -highlightthickness 0
    pack $f.ent -side left -fill x -expand 1 -padx 2
    $f.ent xview end
    if {$i != "htmlopts"} {
      button $f.but -text Browse -padx 0 -pady 0 -command "help_browse $i"
      pack $f.but -side right
    }
  }

  FrameCreate $w.hf -text "CGNS Documentation" -font $Font(bold)
  pack $w.hf -side top -padx 5
  set f [FrameGet $w.hf]
  set lab URL
  if {$HelpData(winhtml)} {append lab /CHM}
  label $f.lab -text $lab -width $lw -anchor w
  pack $f.lab -side left
  entry $f.ent -width $ew -textvariable HelpNew(cgns) -highlightthickness 0
  pack $f.ent -side left -fill x -expand 1 -padx 2
  $f.ent xview end
  button $f.but -text Browse -padx 0 -pady 0 -command "help_browse cgns"
  pack $f.but -side right

  set b [frame $w.but]
  pack $b -side top -pady 5
  button $b.accept -text Accept -width 8 -default active \
    -command "help_check $w"
  button $b.default -text Defaults -width 8 -command help_defaults
  button $b.cancel -text Cancel -width 8 -command "destroy $w"
  pack $b.accept $b.default $b.cancel -side left -padx 5

  bind $w <Return> "help_check $w"

  center_window $w .
  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  catch {grab $w}
  tkwait visibility $w
  focus $w
  tkwait window $w
  catch {focus $oldFocus}
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }
  catch help_menu
}

proc help_browse {what} {
  global HelpData HelpNew tcl_platform
  if {$what == "browser"} {
    if {$tcl_platform(platform) == "windows"} {
      set filelist {{{Executable Files} {.exe .com .bat}}}
    } else {
      set filelist {}
    }
    lappend filelist {{All Files} {*}}
    set fname [FileOpen "Select HTML Browser" $HelpNew(browser) . $filelist]
  } else {
    set filelist {{{HTML Files} {.html .htm}}}
    if {$HelpData(winhtml)} {
      lappend filelist {{CHM Files} .chm}
    }
    if [string match "file://*" $HelpNew($what)] {
      set oldname [string range $HelpNew($what) 7 end]
    } else {
      set oldname $HelpNew($what)
    }
    set fname [FileOpen "CGNS Documentation" $oldname . $filelist]
  }
  if {$fname != ""} {
    if {$tcl_platform(platform) == "windows" &&
      ![catch {file attributes $fname -shortname} name]} {
      set fname $name
    }
    if {$what != "browser"} {
      if {[string tolower [file extension $fname]] == ".chm"} {
        set HelpNew($what) $fname
      } else {
        set HelpNew($what) file://$fname
      }
    }
  }
}

proc help_check {w} {
  global HelpData HelpNew
  set browser [find_file executable $HelpNew(browser)]
  if {$browser == ""} {
    errormsg "can't find HTML browser or it's not executable" $w
    return
  }
  foreach i {browser htmlopts cgns} {
    set HelpData($i) $HelpNew($i)
  }
  if {[info procs tclreg_set] != ""} {
    foreach i {browser htmlopts cgns} {
      catch {tclreg_set Help $i $HelpData($i)}
    }
  }
  destroy $w
}

proc help_valid {} {
  global HelpData
  if {$HelpData(cgns) == ""} {
    return 0
  }
  if {$HelpData(winhtml) &&
      [string tolower [file extension $HelpData($what)]] == ".chm"} {
    return 1
  }
  if {$HelpData(browser) != "" &&
      [file executable $HelpData(browser)]} {
    return 1
  }
  return 0
}

proc help_show {{html ""} {tag ""}} {
  global HelpData tcl_platform

  set htmlroot $HelpData(cgns)
  if {$htmlroot == ""} {
    errormsg "CGNS documentation URL not setup"
    return
  }

  set ext [string tolower [file extension $htmlroot]]
  if {$HelpData(winhtml) && $ext == ".chm"} {
    if {$HelpData(chmfile) != $htmlroot} {
      if {[catch {WinHtml file $htmlroot} msg]} {
        errormsg $msg
        return
      }
      set HelpData(chmfile) $htmlroot
    }
    if {$html == ""} {
      if {[catch {WinHtml index} msg]} {
        errormsg $msg
      }
    } else {
      if {[catch {eval WinHtml topic $html $tag} msg]} {
        errormsg $msg
      }
    }
    return
  }

  if {$HelpData(browser) == ""} {
    errormsg "browser not set up"
    return
  }
  set doc $htmlroot
  if {$ext == ".html" || $ext == ".htm"} {
    set n [string last / $htmlroot]
    if {$n > 0} {set doc [string range $htmlroot 0 $n]}
  }
  if {[string index $doc end] != "/"} {append doc /}
  if {$html == ""} {
    append doc "index.html"
  } else {
    append doc "$html"
  }
  if {$tag != ""} {append doc "\#$tag"}
  set cmd "$HelpData(browser) $HelpData(htmlopts) "
  if {$tcl_platform(platform) == "windows"} {
    append cmd "\"$doc\""
  } else {
    append cmd $doc
  }
  if {[catch {eval exec $cmd &} msg]} {
    errormsg $msg
  }
}

