#TODO double check radian mathmatics!
package require pbctools

#Converts degrees to radians
proc DtoR {d} {
    set pi 3.1415926535897931
    set numerator [expr $d*$pi]
    set out [expr $numerator/180.0]
    return $out
}


;# Centers around the protien and centers the system at 0,0,0
proc Center { stuff } {
    pbc wrap -center com -centersel ${stuff} -orthorhombic  -all
    set nframes [molinfo top get numframes]
    set zero [veczero]
    for {set frames 0} {$frames < $nframes} {incr frames} {
        set sel [atomselect top "all" frame $frames]
        set com [measure center $sel]
        $sel delete
        set sell [atomselect top "all" frame $frames]
        set mov [vecsub $zero $com]
        $sell moveby $mov
        $sell delete
    }
}

;# Alignment based off vmd alignment
proc Align { stuff } {
    set nframes [molinfo top get numframes]
    set ref [atomselect top $stuff frame 0]
    for {set frames 1} {$frames < $nframes} {incr frames} {
        set com [atomselect top $stuff frame $frames]
        set TM [measure fit $com $ref]
        $com delete
        set move_sel [atomselect top "all" frame $frames]
        $move_sel move $TM
        $move_sel delete
    }
    $ref delete
}


proc  Coord_Profile {Xi Xf Yi Yf} {
    if {($Xi < $Xf) && ($Yi < $Yf)} {
        set a "x > $Xi and x < $Xf and y > $Yi and y < $Yf"
    } elseif { ($Xi > $Xf) && ($Yi < $Yf)} {
        set a "x < $Xi and x > $Xf and y > $Yi and y < $Yf"
    } elseif {($Xi < $Xf) && ($Yi > $Yf)} {
        set a "x > $Xi and x < $Xf and y < $Yi and y > $Yf"
    } elseif {($Xi > $Xf) && ($Yi > $Yf)} {
        set a "x < $Xi and x > $Xf and y < $Yi and y > $Yf"
    }
    return $a
}

proc Time_Predict {step_list estimated_steps} { ;#linear, is not really accurate
    set avg 0
    foreach v $step_list { set avg [expr $avg + $v] }
    set avg [expr 1.0*$avg/[llength $step_list]]
    set predict [expr 1.0*$avg*$estimated_steps/60.0]
    puts "\n\nEstimated Time to Run: $predict (min)\n\n"
}
