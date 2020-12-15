mol new FeOcc.pdb waitfor all
set all [atomselect top all]
set frame 0
set in [open FeOcc.pdb r]
set beta {}
while { [gets $in line] != -1 } {
	switch -- [string range $line 0 3] {
    	END {
    		$all frame $frame
    		$all set user $beta
    		set beta {}
    		incr frame
    	}
    	HETA -
    	ATOM {
    		lappend beta [expr [string range $line 60 65]]
    	}
	}
}
