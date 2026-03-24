proc translation {dx dy} {(
    ((1, 0, $dx),
	(0, 1, $dy),
	(0, 0,   1  ))
)}
proc reflect_X {} {(
    ((1,  0, 0),
	(0, -1, 0),
	(0,  0, 0))
)}
proc reflect_Y {} {(
    ((-1, 0, 0),
	( 0, 1, 0),
	( 0, 0, 0))
)}

proc reflect_O {} {(
    ((-1,  0, 0),
	( 0, -1, 0),
	( 0,  0, 1))
)}

proc shear {sx sy} {(
    ((  1 , $sx, 0 ),
	($sy,   1 , 0 ),
	(  0 ,   0 , 0 ))
)}

proc rotation {angle} {(
    angle=double($angle)/180*acos(-1);
    ((cos($angle), -sin($angle), 0),
	( sin($angle), cos($angle), 0),
	(         0      ,          0      , 1))
)}

puts [rotation 45]
puts [translation 10 20]
puts [reflect_X]
puts [shear 3 4]

foreach ns {matrix tensor vector} {
    namespace eval $ns {
	namespace ensemble create
	namespace export *
    }
}

proc tensor::dim {T {depth 0}} {
    expr {i = 0; k = j = {}}
    lmap e $T {
	if {[llength $e] > 1 } \
	    then {(   j = [dim $e [( $depth + 1 )]]   )} \
	    else {
		if {! ([string is double $e] || [string is entier $e])}\
		    then { error -message "argument <$i> de <$T> non numerique" }\
		    else {(  j = {} )}
	}
	if {$i > 1 && $j !=$k} {error -message "Tensor inhomogene, argument <$i> de < $T >"}
	expr {i=$i+1; k = $j}
    }
    return [list $i {*}$j]
}

proc vector::transpose {V} {
    set Res {}
    foreach e $V {
	lappend $Res [($e, 0, 0)]
    }
    return $Res
}
proc vector::is {V} {
    
}
proc matrix::map M {
    set DIM [tensor dim $M]
    set Alpha [list a b c d e f g h i j k l m n o p]
    if {[llength $DIM] == 1} {
	# simple vecteur
	return [dict create dim [(1, $DIM)] map [lrange $Alpha 0 $DIM-1]]
    }

    set MAP {}
    for {( i = j = 0 )} {$i < [lindex $DIM 0]} {(
	 i = $i+1;
	 j = $i * [lindex $DIM 1]
	 )} {    lappend MAP [lrange $Alpha $j [( $j+[lindex $DIM 1]-1)] ]   }
    
    return [dict create dim $DIM map $MAP]
}

proc matrix::transpose {M} {
    # matrix map retourne un dictionnaire dont les cles sont dim et map
    lmap {var val} [matrix map $M] {($var = $val)}
   
    if {1 in $dim} { return [vector transpose $M] }\
	elseif {[llength $dim] > 2} { error -message "$M n'est pas une matrice"}\
	else {
	    lmap r $map R $M {     lmap var $r val $R {( $var = $val )}   }
	    return [switch $dim {
    	        {2 2} {(   ($a, $c), ($b, $d)   )}
    	        {3 2} {(   ($a, $c, $e) , ($b, $d, $f)   )}
    	        {2 3} {(   ($a, $d) , ($b, $e) , ($c, $f) )}
    	        {3 3} {(   ($a, $d, $g) , ($b, $e, $h) , ($c, $f, $i)  )}
    	        {3 4} {(   ($a, $e, $i) , ($b, $f, $j) , ($c, $g, $k) , ($d, $h, $l)   )}
    	        {4 3} {(   ($a, $d, $g, $j) , ($b, $e, $h, $k) , ($c, $f, $i, $l)  )}
    	        {4 4} {(   ($a, $e, $i, $m) , ($b, $f, $j, $n) , ($c, $g, $k, $o) , ($d, $h, $l, $p)  )}
	    }]
	}}

puts [matrix transpose {{1 2 3 4} {5 6 7 8} {9 10 11 12} {13 14 15 16}}]
puts [matrix map {{1 2 3 4} {5 6 7 8} {9 10 11 12} {13 14 15 16}}]

proc matrix::product {M1 M2} {
    # definition de dim1 et map1
    foreach {var val} [matrix map $M1] {("${var}1" = $val)}
    # definition de dim2 et map2
    foreach {var val} [matrix map $M2] {("${var}2" = $val)}

    
    if {[llength $dim1] > 2 || [llength $dim2] > 2} {error -message "non matrices"}
    if {[lindex $dim1 1] != [lindex $dim2 0]} {
	if {1 in $dim1} {
	    #M1 est un vecteur
	} elseif {1 in $dim2} {
	    # M2 est un vecteur
	} else {
	    error -message "matrice incompatibles"
	}
    }
 
    foreach r $map1 R $M1 {     foreach var $r val $R {( "${var}1" = $val )}   }
    foreach r $map2 R $M2 {     foreach var $r val $R {( "${var}2" = $val )}   }
    set DIMs [list $dim1 $dim2]

    return [switch $DIMs {
	{{1 3} {3 3}} {(
	    (($a1*$a2 + $b1*$d2 + $c1*$g2),
		($a1*$b2 + $b1*$e2 + $c1*$h2),
		($a1*$c2 + $b1*$f2 + $c1*$i2))
	       )}
	{{3 3} {1 3}} {(
	    (($a1*$a2 + $b1*$b2 + $c1*$c2),
		($d1*$a2 + $e1*$b2 + $f1*$c2),
		($g1*$a2 + $h1*$b2 + $i1*$c2))
	       )}
	{{2 2} {2 2}} {(
	    (( $a1 * $a2 + $b1 * $c2, $a1 * $b2 + $b1 * $d2),
		( $c1 * $a2 + $d1 * $c2, $c1 * $b2 + $d1 * $d2))
		 )} 
	{{2 3} {3 2}} {(
	    (($a1*$a2 + $b1*$c2 + $c1*$e2 , $a1*$b2 + $b1*$d2 + $c1*$f2),
		($d1*$a2 + $e1*$c2 + $f1*$e2, $d1*$b2 + $e1*$d2 + $f1*$f2))
		)}
	{{3 2} {2 3}} {(
	    (($a1*$a2 + $b1*$d2 , $a1*$b2 + $b1*$e2 , $a1*$c2 + $b1*$f2),
		      ($c1*$a2 + $d1*$d2 , $c1*$b2 + $d1*$e2 , $c1*$c2 + $d1*$f2),
		($e1*$a2 + $f1*$d2 , $e1*$b2 + $f1*$e2 , $e1*$c2 + $f1*$f2))
		)}
	{{3 3} {3 3}} {(
	    (($a1*$a2 + $b1*$d2 + $c1*$g2, $a1*$b2 + $b1*$e2 + $c1*$h2, $a1*$c2 + $b1*$f2 + $c1*$i2),
		($d1*$a2 + $e1*$d2 + $f1*$g2, $d1*$b2 + $e1*$e2 + $f1*$h2, $d1*$c2 + $e1*$f2 + $f1*$i2),
		($g1*$a2 + $h1*$d2 + $i1*$g2, $g1*$b2 + $h1*$e2 + $i1*$h2, $g1*$c2 + $h1*$f2 + $i1*$i2))
		)}
	defaut {($DIMs)}
    }]
}

puts MP:[matrix product {{1 2 3} {4 5 6} {7 8 9}} {{10 11 12} {13 14 15} {16 17 18}}]
puts MP:[matrix product {1 2 3} {{10 11 12} {13 14 15} {16 17 18}}]

set Point [(100*sqrt(2)/2, 100*sqrt(2)/2, 1)]

puts AffineRotation:[matrix product [rotation +45] $Point ]
puts Affine:[matrix product [translation 0 -100] [matrix product [rotation +45] $Point ]]
set AffineTransform [matrix product [translation 0 -100] [rotation +45]]
puts Affine2:[matrix product $AffineTransform $Point ]

proc matrix::determinant {M} {
    foreach {var val} [matrix map $M] {($var = $val)}
    foreach r $map R $M {     foreach var $r val $R {( $var = $val )}   }
    
    return [switch $dim {
	{3 3} {($a*$e*$i + $b*$f*$g + $c*$d*$h - $c*$e*$g - $b*$d*$i - $a*$f*$h)}
    }]
}

puts determinant:[matrix determinant {{1 2 5} {4 5 6} {7 8 9}}]

proc matrix::comatrix {M} {
    foreach {var val} [matrix map $M] {($var = $val)}
    foreach r $map R $M {     foreach var $r val $R {( $var = double($val) )}   }

    return [switch $dim {
	{3 3} {(
	    (($e*$i - $f*$h, $f*$g - $d*$i, $d*$h - $e*$g),
		($c*$h - $b*$i, $a*$i - $c*$g, $b*$g - $a*$h),
		($b*$f - $c*$e, $c*$d - $a*$f, $a*$e - $b*$d))
	)}
    }]
}

puts comatrice:[matrix comatrix {{1 2 5} {4 5 6} {7 8 9}}]
puts comatrice\ transposee:[matrix transpose [matrix comatrix {{1 2 5} {4 5 6} {7 8 9}}]]

proc matrix::inverse {M} {
    set TCOM  [transpose [comatrix $M]]
    set det [determinant $M]
    if {$det == 0} {error "matrix non inversible"}
    return [lmap row1 $TCOM {
	lmap e1 $row1 {(double($e1/$det))}
    }]
}

puts inverse:[matrix inverse {{1 2 5} {4 5 6} {7 8 9}}]

proc tcl::mathfunc::if {args} { tailcall ::if {*}$args }
proc tcl::mathfunc::for {args} { tailcall ::for {*}$args }
proc tcl::mathfunc::while {args} { tailcall ::for {*}$args }
proc tcl::mathfunc::variable {args} { tailcall ::variable {*}$args }

expr {
      n2=3;
      n=2;
      ynm1 = 4;
      if ($n2==$n+1, {(result=$ynm1)})
  }
puts result:$result

namespace eval ::math {
    namespace eval special {
	    variable pi 3.14
    }}

proc ::math::special::J1/2 {x} {(
     variable("pi");
    #
    # This Bessel function can be expressed in terms of elementary
    # functions. Therefore use the explicit formula
    #
     if( $x != 0.0, {( sqrt(2.0/$pi/$x)*sin($x) )}, {( 0.0 )})
)}

puts J1/2:[::math::special::J1/2 3]

proc ::math::special::Jn {n x} {(
    variable("pi");

    ##nagelfar ignore
    if(![string is integer -strict $n], {
	return -code error "Order argument must be integer"
    });

    #
    # Integrate over the interval [0,pi] using a small
    # enough step - 40 points should do a good job
    # with |x| < 20, n < 20 (an accuracy of 1.0e-8
    # is reported by Davis and Rabinowitz)
    #
     number = 40;
     step = $pi/double($number);
     result = 0.0;
 
     for ({(i=0)}, {$i <= $number}, { incr i }, {(
         t = double($i)*$step;
         f  = cos($x * sin($t) - $n * $t);
         if($i == 0 || $i == $number, {(
               f = $f/2.0
         )});
        result =$result+$f
      )});

    $result*$step/$pi
)}

puts jn:[::math::special::Jn 3 5]
    
proc integralSH { begin end nosteps func } {(
    delta = ($end-$begin)/double($nosteps);
    hdelta = $delta/2.0;
    result = 0.0;
    xval = $begin;
    func_end = [uplevel 1 [($func, $xval)]];
    for({(i=1)},{$i <= $nosteps },{ incr i },{(
	func_begin = $func_end;
	func_middle = [uplevel 1 [($func, $xval+$hdelta)]];
	func_end = [uplevel 1 [($func, $xval+$delta)]];
	result = $result+$func_begin+4.0*$func_middle+$func_end;
	xval = $begin+double($i)*$delta
	)});

    $result*$delta/6.0
    )}

proc integralSH5 { begin end nosteps func } {
    expr {
	  delta = ($end-$begin)/double($nosteps);
	  hdelta = $delta/2.0;
	  result = 0.0;
	  xval = $begin;
	  func_end = [uplevel 1 [list $func $xval]]
      }
    for {(i=1)} {$i <= $nosteps } {(i=$i+1)} {
	expr {
	      func_begin = $func_end;
	      func_middle = [uplevel 1 [list $func [expr {$xval+$hdelta}]]];
	      func_end = [uplevel 1 [list $func [expr {$xval+$delta}]]];
	      result = $result+$func_begin+4.0*$func_middle+$func_end;
	      xval = $begin+double($i)*$delta
	  }
    }

    return [($result*$delta/6.0)]
}

proc integralSH6 { begin end nosteps func } {
    expr {
	  delta = ($end-$begin)/double($nosteps);
	  hdelta = $delta/2.0;
	  result = 0.0;
	  xval = $begin;
	  func_end = [uplevel 1 [list $func $xval]]
      }
    for {set i 1)} {$i <= $nosteps } {incr i} {
	expr {
	      func_begin = $func_end;
	      func_middle = [uplevel 1 [list $func [expr {$xval+$hdelta}]]];
	      func_end = [uplevel 1 [list $func [expr {$xval+$delta}]]];
	      result = $result+$func_begin+4.0*$func_middle+$func_end;
	      xval = $begin+double($i)*$delta
	  }
    }

    return [($result*$delta/6.0)]
}

proc integralSH3 { begin end nosteps func } {
    expr {
	  delta = ($end-$begin)/double($nosteps);
	  hdelta = $delta/2.0;
	  xval = $begin;
	  func_end = [uplevel 1 [($func, $xval)]]
      }  
    
    for {( result = 0.0; i=1)} { $i <= $nosteps } {(i=$i+1)} {(
      	func_begin  = $func_end;
	func_middle= [uplevel 1 [($func, $xval+$hdelta)]];
	func_end    = [uplevel 1 [($func, $xval+$delta)]];
	result =$result+$func_begin+4.0*$func_middle+$func_end;
	xval=$begin+double($i)*$delta
    )}
    return [($result*$delta/6.0)]
}

proc integralSH2 { begin end nosteps func } {

   set delta    [( ($end-$begin)/double($nosteps) )]
   set hdelta   [( $delta/2.0 )]
   set result   0.0
   set xval     $begin
   set func_end [uplevel 1 [($func, $xval)]]
   for { set i 1 } { $i <= $nosteps } { incr i } {
       set func_begin  $func_end
       set func_middle [uplevel 1 [($func, $xval+$hdelta)]]
       set func_end    [uplevel 1 [($func, $xval+$delta)]]
       set result      [($result+$func_begin+4.0*$func_middle+$func_end)]
       set xval        [($begin+double($i)*$delta)]
   }

   return [($result*$delta/6.0)]
}

proc integralSH1 { begin end nosteps func } {

   set delta    [( ($end-$begin)/double($nosteps) )]
   set hdelta   [( $delta/2.0 )]
   set result   0.0
   set xval     $begin
   set func_end [uplevel 1 [list $func $xval]]
   for { set i 1 } { $i <= $nosteps } { incr i } {
       set func_begin  $func_end
       set func_middle [uplevel 1 [list $func [($xval+$hdelta)]]]
       set func_end    [uplevel 1 [list $func [($xval+$delta)]]]
       set result      [($result+$func_begin+4.0*$func_middle+$func_end)]
       set xval        [($begin+double($i)*$delta)]
   }

   return [($result*$delta/6.0)]
}

proc integral { begin end nosteps func } {

   set delta    [expr {($end-$begin)/double($nosteps)}]
   set hdelta   [expr {$delta/2.0}]
   set result   0.0
   set xval     $begin
   set func_end [uplevel 1 [list $func $xval]]
   for { set i 1 } { $i <= $nosteps } { incr i } {
      set func_begin  $func_end
       set func_middle [uplevel 1 [list $func [expr {$xval+$hdelta}]]]
       set func_end    [uplevel 1 [list $func [expr {$xval+$delta}]]]
      set result      [expr  {$result+$func_begin+4.0*$func_middle+$func_end}]

      set xval        [expr {$begin+double($i)*$delta}]
   }

   return [expr {$result*$delta/6.0}]
}

proc f {x} {
    return $x
}
puts Orig:[timerate {integral -2 2 10 f }]
puts length:[(LO=[string length [info body integral]] )]=100%
puts SH1:[timerate {integralSH1 -2 2 10 f }]
puts length:[(L = [string length [info body integralSH1]],  "=", $L/double($LO)*100-100)]%
puts SH2:[timerate {integralSH2 -2 2 10 f }]
puts length:[(L = [string length [info body integralSH2]],  "=", $L/double($LO)*100-100)]%
puts SH3:[timerate {integralSH3 -2 2 10 f }]
puts length:[(L = [string length [info body integralSH3]],  "=", $L/double($LO)*100-100)]%
puts SH:[timerate {integralSH -2 2 10 f }]
puts length:[(L = [string length [info body integralSH]],  "=", $L/double($LO)*100-100)]%
puts SH5:[timerate {integralSH5 -2 2 10 f }]
puts length:[(L = [string length [info body integralSH5]],  "=", $L/double($LO)*100-100)]%
puts SH6:[timerate {integralSH5 -2 2 10 f }]
puts length:[(L = [string length [info body integralSH6]],  "=", $L/double($LO)*100-100)]%

# puts Orig:[time {integral -2 2 10 f } 25 ]
# puts SH4:[time {integralSH4 -2 2 10 f } 25]
# puts SH3:[time {integralSH3 -2 2 10 f } 25]
# puts SH2:[time {integralSH2 -2 2 10 f } 25]
# puts SH:[time {integralSH -2 2 10 f } 25]


proc integralExprSH { begin end nosteps expression } {
   
    expr {
	  hdelta = (delta = ($end-$begin)/double($nosteps)) /2.0;
	  result = 0.0;
	  x = $begin;
	  func_end = [expr $expression]
    }
    # FRINK: nocheck

    for { set i 1 } { $i <= $nosteps } { incr i } {(
	      func_begin = $func_end;
	      x=$x+$hdelta;
	      func_middle=[expr $expression];
	      x=$x+$hdelta;
	      func_end=[expr $expression];
	      result = $result+$func_begin+4.0*$func_middle+$func_end;
	      x = $begin+double($i)*$delta
    )}

    return [($result*$delta/6.0)]
}

proc integralExpr { begin end nosteps expression } {

   set delta    [expr {($end-$begin)/double($nosteps)}]
   set hdelta   [expr {$delta/2.0}]
   set result   0.0
   set x        $begin
   # FRINK: nocheck
   set func_end [expr $expression]
   for { set i 1 } { $i <= $nosteps } { incr i } {
      set func_begin  $func_end
      set x           [expr {$x+$hdelta}]
       # FRINK: nocheck
      set func_middle [expr $expression]
      set x           [expr {$x+$hdelta}]
       # FRINK: nocheck
      set func_end    [expr $expression]
      set result      [expr {$result+$func_begin+4.0*$func_middle+$func_end}]

      set x           [expr {$begin+double($i)*$delta}]
   }

   return [expr {$result*$delta/6.0}]
}


puts Orig:[timerate [list integralExpr  0 4 100 {$x**2}]]
puts SH:[timerate [list integralExprSH 0 4 100 {$x**2}]]

puts [(x=1, x=$x+.5, x=$x+1)]
