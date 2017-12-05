#!/bin/perl

use strict;
use warnings;

my $csfastq = shift;
die unless defined($csfastq);
my $csfasta = $csfastq; $csfasta =~ s/csfastq$/csfasta/;
die unless !($csfastq eq $csfasta);
my $qual = $csfastq; $qual =~ s/.csfastq$/_QV.qual/;
die unless !($csfastq eq $qual);

open(FHcsfastq, "$csfastq") || die;
open(FHcsfasta, ">$csfasta") || die;
open(FHqual, ">$qual") || die;
my $state = 0;
my ($n, $r, $q) = ("", "", "");
while(defined(my $line = <FHcsfastq>)) {
    chomp($line);
    if(0 == $state) {
        &print_out(\*FHcsfasta, \*FHqual, $n, $r, $q);
        $n = $line;
        $n =~ s/^\@/>/;
    }
    elsif(1 == $state) {
        $r = $line;
    }
    elsif(3 == $state) {
        $q = $line;
        # convert back from SANGER phred
        my $tmp_q = "";
        for(my $i=0;$i<length($q);$i++) {
            my $Q = ord(substr($q, $i, 1)) - 33;
            die unless (0 <= $Q);
            if(0 < $i) {
                $tmp_q .= " ";
            }
            $tmp_q .= "$Q";
        }
        $q = $tmp_q;
    }
    $state = ($state+1)%4;
}
&print_out(\*FHcsfasta, \*FHqual, $n, $r, $q);
close(FHcsfasta);
close(FHcsfastq);
close(FHqual);

sub print_out {
    my ($FHcsfasta, $FHqual, $n, $r, $q) = @_;

    if(0 < length($n)) {
        print $FHcsfasta "$n\n$r\n";
        print $FHqual "$n\n$q\n";
    }
}
