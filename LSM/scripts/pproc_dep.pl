#!/usr/bin/perl

while(<>) {

    $MOD = "mod"; # could be different for other compilers

    if ( /\s*(\w+)\.$MOD:\s+(\w+)\.o\s*$/ ) {
	print "$1\@$2.smod: $2.o\n";
	print "$1.$MOD: $1\@$2.smod\n";
    } else {
	print ;
    }


}
