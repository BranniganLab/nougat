proc print_director_files { mol_id outfile fit_mode normalization_mode } {
  #assign_directors $mol_id $fit_mode $normalization_mode;
  set f [open $outfile "w"]
  set n [molinfo $mol_id get numframes]
  
  #assumes only lipids are in the system and is model dependent,  needs to be generalized
  set sel [atomselect $mol_id "name PO4"] 
  set reslist [$sel get resid]
  puts "[$sel num]"
  $sel delete
  
  for {set i 0} {$i < $n} {incr i} {
    foreach res $reslist {
      set resSel [atomselect $mol_id "resid $res and name PO4"]
      $resSel frame $i
      set director [$resSel get {user user2 user3}]
      puts $f "[lindex $director 0] [lindex $director 1] [lindex $director 2]"
      $resSel delete
    }
  }
  close $f 
}


proc assign_directors { mol_id fit_mode normalization_mode } {
  if {$fit_mode == 0} {
    puts "Using least squares fitting."
  } elseif {$fit_mode == 1} {
    puts "Using PO4 + C2 and C4A C4B"
  }
  set n [molinfo $mol_id get numframes]
  #assumes only lipids are in the system and is model dependent,  needs to be generalized
  set sel [atomselect $mol_id "name PO4"] 
  set reslist [$sel get resid]
  $sel delete
  foreach res $reslist {
    set resSel [atomselect $mol_id "resid $res"]
    set end1_sel [atomselect $mol_id "resid $res and name PO4 C2"]
    set end2_sel [atomselect $mol_id "resid $res and name C4A C4B"]
    for {set i 0} {$i < $n} {incr i} {
      $resSel frame $i
      $end1_sel frame $i
      $end2_sel frame $i
      #      puts "frame $i"
      if {$fit_mode == 0} {
        set director [fit_director_to_sel $resSel $normalization_mode]
      } elseif {$fit_mode == 1} {
        set end1_com [measure center $end1_sel]
        set end2_com [measure center $end2_sel]
        set dist_vector [vecsub $end1_com $end2_com]
        if {$normalization_mode == 0} {
          set normalization [expr 1.0/[lindex $dist_vector 2]]
          set director [vecscale $normalization $dist_vector]
        } else {
          set director [vecnorm $dist_vector]
        }
      }
      #puts $director
      $resSel set user [lindex $director 0]
      $resSel set user2 [lindex $director 1]
      $resSel set user3 [lindex $director 2]
    }
    $resSel delete
    $end1_sel delete
    $end2_sel delete
  }
}


proc fit_director_to_sel { sel normoption} {
  set fitdir  [list [lsq [$sel get x]] [lsq [$sel get y]] [lsq [$sel get z]]]
  if {$normoption == 0} {
    set normalization [expr 1.0/[lindex $fitdir 2]]
    set scaled_fit_dir [vecscale $normalization $fitdir]
  } else {
    set scaled_fit_dir [vecnorm $fitdir]
  }
  
  return $scaled_fit_dir
}

proc fit_vec_to_sel { sel } {
  set fitvec [ vecnorm [list [lsq [$sel get x]] [lsq [$sel get y]] [lsq [$sel get z]]]]
  return $fitvec
}

#  Fit the points x to x = ai + b, i=0...N-1, and return the value of a 
# a = 12/( (N(N^2 - 1)) ) sum[ (i-(N-1)/2) * xi]
# reference: Bevington
proc lsq { x } {
  set N [llength $x]
  set xtot 0
  set d [expr {0.5*($N-1)}]
  
  set i 0.0
  foreach elem $x {
    set xtot [expr {$xtot + ($i - $d) * $elem}]
    set i [expr {$i + 1.0}]
  }
  #puts $xtot
  return $xtot
}

#  least squares procedure from above, optimized for use with vecexpr
proc lsq_vecexpr { tail_length list_of_tail_coords } {
  set d [expr {0.5*($tail_length-1)}]
  set i_list []
  for {set i 0} {$i<$tail_length} {incr i} {
    lappend i_list $i
  }
  vecexpr $i_list $d sub >multiplier
  set vector [vecexpr [vecexpr $list_of_tail_coords $multiplier mult] sum]
  return $vector
}