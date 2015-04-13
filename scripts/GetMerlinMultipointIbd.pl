#!/usr/bin/perl -w

@mkrinfo = ();

%markers = ();
%sex = ();

while(<>) {
    chomp;
    if (/^(\d+)\s+(\d+)\s+(\d+)\s+([\d]+)\.\d+\s+([-\.\d]+)\s+([-\.\d]+)\s+([-\.\d]+)/) {
#	$ped = $1;
	$id1 = $2;
	$id2 = $3;
	$mkr = $4;
	$ibd = $7 + .5*$6;
	
	$markers{$mkr} = 1;

	push (@{$mkrinfo[$mkr]}, $id1 . "  ". $id2 . "  " . $ibd);
    }
}


sub ReadMerlinPedFile {
    open(PED, "merlin.ped");
    
    while(<PED>) {
	@array = split(/\s+/, $_);
	$sex{$array[1]} = $array[4];
    }

    close(PED);
}


ReadMerlinPedFile();

foreach $marker (keys(%markers)) {
    open(MM, ">mibd-mm.X.$marker");
    open(MF, ">mibd-mf.X.$marker");
    open(FF, ">mibd-ff.X.$marker");

    foreach $_ (@{$mkrinfo[$marker]}) {	
	@array = split(/\s+/, $_);
	next if ($array[2] == 0.0);
	if ($sex{$array[0]}==1 && $sex{$array[1]}==1) {
	    print MM $_ . "\n";
	}
	elsif ($sex{$array[0]}==2 && $sex{$array[1]}==2) {
	    print FF $_ . "\n";
	}
	else {
	    # Remember to scale to ibd by two
	    $array[2] *= 2;
	    print MF join("  ",@array) . "\n";	    
	}
    }

#    print(join("\n", @{$mkrinfo[$marker]}));

    close(FF);
    close(MF);
    close(MM);
}


