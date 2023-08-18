#!/usr/bin/perl -w


require ("funclib.pl");
use strict;
use Getopt::Long;


##########################################################################################
#  parameters & usage
##########################################################################################
my %opts=();
GetOptions(\%opts,"h","in=s","poly=s","ml=i","mp=i","mg=i","mm=i","mtail:i","mper:s","mr:i","odir:s","suf:s","reg:i","deep:s","oraw:s","debug:s","bar:s","review:s");
my $USAGE=<< "_USAGE_";
Require: funclib.pl
Usage: 
  0) Test different parameters -debug T open
  -debug T -oraw T

  1) Recommended parameters for PATseq
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly T -ml 25 -mp 6 -mg 5 -mm 2 -mr 2 -mper 0.75 -mtail 8 -deep T -reg 1 -odir "f:/" -suf "" 
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly T -ml 25 -mp 8 -mg 8 -mm 2 -odir "f:/" -suf "" -reg 1
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly "A|T" -ml 25 -mp 8 -mg 8 -mm 2 -odir "f:/" -suf "AT" -reg 1

do not output the original sequence
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly "A&T" -ml 25 -mp 8 -mg 8 -mm 2 -odir "f:/" -suf "A&T" -reg 1 -oraw F

  debug Output various situations
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly "A&T" -ml 25 -mp 8 -mg 8 -mm 2 -odir "f:/" -suf "A&T" -reg 1 -oraw F -debug T

-h=help
-in=input fa or fq file
-poly=A/T/A&T/A|T
  A|T If there are both A and T，when tail_lenA>=lenT, keep the result of A; a sequence can only contain A or T, not both
  A&T Find A and T at the same time, the two are irrelevant, a sequence may have both A and T at the same time
-ml=min length after trim
-mg=margin from the start (poly=T) or to the end (poly=A) (=5)
-mm=mismatch between TxxTTT (=2)
-mr=minT in reg (=3)
-mp=min length of succesive poly (=8)
-mtail=min length of trimmed tail (=8)
-mper=min percent of A/T in trimmed tail (=0.75)
-deep=T(default)/F If T, when the reg does not match, search in depth, such as this type TTTTTCTTTTCTCTTTTTTTT
-odir=output path (default is the same as input)
-suf=suffix, default: xx.suf.T/A.fq
-reg=1(loose, default)/2
  1='^.{0,mg}?(T{mr,}[^T]{0,mm})*T{mp,}([^T]{0,mm}T{mr,})*'
  2='^.{0,mg}?(T{mp,}([^T]{0,mm}T{mr,})*)'
-oraw=T(default)/F; output raw reads for trimmed A/T file (.A.raw.fq)
-bar=Last bast position of three barcode 
-review=T(default)/If F is T, when reg/deep does not match, search in depth, such as AACCCCTTTTTTTTTTTTTTTTT#seven
_USAGE_


#############################################################################
#  invoke the program                               
#############################################################################
die $USAGE if $opts{'h'}||!$opts{'in'}||!$opts{'poly'}||!$opts{'ml'}||!$opts{'mp'}||!$opts{'mg'}||!$opts{'mm'};

my $poly=$opts{'poly'};
my $in=$opts{'in'};
my $ml=$opts{'ml'};
my $mp=$opts{'mp'};
my $mg=$opts{'mg'};
my $mm=$opts{'mm'};
my $mr=$opts{'mr'};
my $mtail=$opts{'mtail'};
my $mper=$opts{'mper'};
my $deep=$opts{'deep'};
my $debug=$opts{'debug'};
my $suf=$opts{'suf'};
my $odir=$opts{'odir'};
my $reg=$opts{'reg'};
my $oraw=$opts{'oraw'};
my $bar=$opts{'bar'};
my $review=$opts{'review'};
$reg=1 if !$reg;
die "reg=1/2" if $reg!=1 and $reg!=2;


if (!$mtail) {
  $mtail=8;
}
if (!$mper) {
  $mper=0.75;
}
if (!$mr) {
  $mr=3;
}
if (!$deep | $deep eq 'T') {
  $deep=1;
} else {
  $deep=0;
}
if (!$review | $review eq 'T') {
  $review=1;
} else {
  $review=0;
}
if (!$debug | $debug ne 'T') {
  $debug=0;
} else {
  $debug=1;
}
if(!$bar){
 $bar=0;	
}
$oraw=1 if !$oraw or $oraw ne 'F';
$oraw=0 if $oraw eq 'F';

print "poly=$poly\ninfile=$in\noutput dir=$odir\nmin length after trim(ml)=$ml\nmin continue tail length(mp)=$mp\nmax margin to the end(mg)=$mg\nmax mismatches between xTTT(mm)=$mm\n";
print "min reg T(mr)=$mr\nmin tail length(mtail)=$mtail\nmin T/A percent(mper)=$mper\n";
print "deep=$deep\n";
print "DEBUG MODE, output to .test.file too \n" if $debug;

die "error poly(=A T A|T A&T)" if $poly ne 'A' and $poly ne 'T' and $poly ne 'A|T' and $poly ne 'A&T';
my ($findT,$findA)=(0,0);
$findA=1 if ($poly=~ m/A/);
$findT=1 if ($poly=~ m/T/);
#print "$findT\n";

if ($deep & $findA) {
  die "TODO: deepand findA's program has not yet been implemented...";
}

die "only support fastq file" if seqFormat($in) ne 'fq';
open my $fh1, "<$in" or die "cannot read input file ($in)\n";
#open IN, "<$in" or die "cannot read input file ($in)\n";
$suf=".$suf" if $suf;

my ($cntNotailT,$cntShortT,$cntTotal,$cntFinalT,$cntNotailA,$cntShortA,$cntFinalA,$cntMissA,$cntMissT,$cntBadT,$cntBadA)=(0,0,0,0,0,0,0,0,0,0,0);

my $inpath=getDir($in,0); # x/y/
my $infname=getFileName($in,1); #File prefix
if (!$odir) {
  $odir=$inpath;
}

if ($findA) {
  open OA,">$odir$infname$suf.A.fq";
  if($debug){
  	open TESTA,">$odir$infname$suf.testA.fq";
  	}
  if($oraw){
  	open RAWA,">$odir$infname$suf.A.raw.fq" ;
  	}
}


if ($findT) {
  open OT,">$odir$infname$suf.T.fq";#Effective result, #seven with polyT tail removed
  if($debug){
  	open TESTT,">$odir$infname$suf.T.DEBUG.fq";
  	open BEDT,">$odir$infname$suf.badT.fq";
  	open NOT,">$odir$infname$suf.noTail.fq";
  	open SHORTT,">$odir$infname$suf.shortT.fq";
  	}
  if($oraw){
  	open RAWT,">$odir$infname$suf.T.raw.fq";
  	}
}

#print "###############$findT\t$debug\t$oraw\n";

my $regT="^.{0,$mg}?(T{$mr,}[^T]{0,$mm})*T{$mp,}([^T]{0,$mm}T{$mr,})*";
if ($reg==2) {
  $regT="^.{0,$mg}?T{$mp,}([^T]{0,$mm}T{$mr,})*";
}
print "regT=$regT\n";
my $regA="(A{$mr,}[^A]{0,$mm})*A{$mp,}([^A]{0,$mm}A{$mr,})*.{0,$mg}?\$";
if ($reg==2) {
  $regA="A{$mp,}([^A]{0,$mm}A{$mr,})*.{0,$mg}?\$";
}
print "regA=$regA\n";

if ($deep) {
  print "deep reg=^.{0,$mg}?(T{1,}[^T]{0,2})*T{8,}\n";
}
if($review){
   print "review_reg=T{$mp,}([^T]{0,$mm}T{$mr,}){0,}\n";	
}

my ($md1,$md2,$md3)=(0,0,0);


my @fqkeys = qw/id seq pl qual/;
while(!eof $fh1) {
	my (%e1);
	@e1{@fqkeys} = readfq($fh1);#hash stores the id of reads at the same time、seq、pl、qual
	$cntTotal++;
	my ($seqT,$seqA,$qualT,$qualA)=('','','','');
	my ($haveT,$haveA,$shortT,$notailT,$shortA,$notailA,$badTailT,$badTailA)=(0,0,0,0,0,0,0,0);
	if ($findT) {
	    $seqT=$e1{'seq'};
 		if ($seqT=~s/$regT//) {# Replace the matching pattern
		  my $length=length($seqT);
		  #my $Ts=substr($e1{'seq'},0,length($e1{'seq'})-$length); #Intercepted T segment#WXH
		  my $Ts=substr($e1{'seq'},$bar,length($e1{'seq'})-$length-$bar);#It should be called barcode，adapter，Seven 
		  #print "$Ts\n";#Seven test
		  #my $Tcnt = $Ts =~ tr/T/T/; #计算T数
		  my $Tcnt = $Ts=~ tr/T/T/;#计算T数，Seven
		  ##print "Yes: $Ts.$Tcnt\n";
		  if ((length($e1{'seq'})-$length-$bar)<$mtail or $Tcnt/length($Ts)<$mper) { #tail过短 or T%<mper
			$badTailT=1;
			print TESTT $Ts."***BAD***".$seqT."\n"  if $debug;
		  } elsif ($length<$ml) { #The sequence is too short
			$shortT=1;
			print TESTT $Ts."***SHR***".$seqT."\n"  if $debug;
		  } else {
		  	++$md1;
			$haveT=1;
		    $qualT=substr($e1{'qual'},length($e1{'qual'})-$length,$length);#Meet the conditions seven
		  }
		} 

       if (!$haveT & $deep) { #Deep search, judgment fixed regular formula，Find and intercept to the beginning, and then judge mper.. (TTTTTCTTTTCTCTTTTTTTT)
          $seqT=$e1{'seq'};
		      my $deepreg="^.{0,$mg}?(T{1,}[^T]{0,2})*T{$mp,}";
		      if ($seqT=~s/$deepreg//){
			   		 my $length=length($seqT);
			 			 #my $Ts=substr($e1{'seq'},0,length($e1{'seq'})-$length); #Intercepted T fragment WXH
			 			 my $Ts=substr($e1{'seq'},$bar,length($e1{'seq'})-$length-$bar);#It should be culled barcode，adapter，Seven
			 			   #print "$Ts\n";
			  	  my $Tcnt = $Ts =~ tr/T/T/; #计算T数
		      		##print "Deep: $Ts.$Tcnt\n";
			      if ((length($e1{'seq'})-$length-$bar)>=$mtail and $Tcnt/length($Ts)>=$mper and $length>=$ml) { 
				      ++$md2;
				      $haveT=1;
				      $qualT=substr($e1{'qual'},length($e1{'qual'})-$length,$length);
				      print TESTT $Ts."***DEP***".$seqT."\n"  if $debug;
				      $badTailT=0;
				      $shortT=0;
			       }	  
		     }
	     }
	     
	    # Check if conditions are met for processing DNA sequence data
if (!$haveT && $review) {  # If not already processed and review is true
    # Find in-depth degree of depth, as long as it matches specific criteria
    $seqT = $e1{'seq'};
    my $deepmod = "T{$mp,}([^T]{0,$mm}T{$mr,}){0,}";  # Regular expression pattern
    
    # Check if the pattern matches in the sequence
    if ($seqT =~ m/$deepmod/) {
        my $pos_start = $-[0];  # Start position of the match
        my $match_length = length $&;  # Length of the matched portion
        my $pos_end = $pos_start + $match_length;  # End position of the match
        
        # Extract the portion after the match in the sequence
        $seqT = substr($seqT, $pos_end, length($seqT) - $pos_end);
        
        my $length = length($seqT);  # Length of the remaining sequence
        
        # Extract the portion before the match in the 'seq' of the hash 'e1'
        my $Ts = substr($e1{'seq'}, $bar, $pos_end - $bar);  # Extracted sequence
        
        # Count the occurrences of 'T' in the extracted sequence
        my $Tcnt = $Ts =~ tr/T/T/;  # Count of 'T' occurrences
        
        # Check if conditions are met for further processing
        if (length($Ts) >= $mtail && $Tcnt / length($Ts) >= $mper && $length >= $ml) {
            ++$md3;  # Increment counter
            $haveT = 1;  # Set flag indicating 'T' sequence is present
            
            # Extract the portion after the match in the 'qual' of the hash 'e1'
            $qualT = substr($e1{'qual'}, $pos_end, length($e1{'qual'}) - $pos_end);
            
            # Debugging: print the sequences to a file
            if ($debug) {
                print TESTT $Ts . "***DEPSEVEN***" . $seqT . "\n" if $review;
            }
            
            $badTailT = 0;  # Initialize variables
            $shortT = 0;
        }
    }
}

	     

	  if (!$haveT) {
    $notailT = 1;
    if ($debug) {
        # print TESTT "****NOT***$e1{'seq'}\n";#seven
    }
} elsif ($debug) {
    # print TESTT "****YES***$e1{'seq'}\n";
}
# End the process of finding poly(T)

if ($findA) {
    $seqA = $e1{'seq'};
    if ($seqA =~ s/$regA//) {
        my $length = length($seqA);
        print TESTA "$seqA*******" . substr($e1{'seq'}, $length, length($e1{'seq'}) - $length) . "\n" if $debug;
        my $As = substr($e1{'seq'}, $length, length($e1{'seq'}) - $length); # Extracted A segment
        my $Acnt = $As =~ tr/A/A/; # Count the occurrences of A
        if (length($e1{'seq'}) - $length < $mtail or $Acnt / length($As) < $mper) { # If tail is too short
            $badTailA = 1;
        } elsif ($length < $ml) {
            $shortA = 1;
        } else {
            $haveA = 1;
            $qualA = substr($e1{'qual'}, 0, $length);
        }
    } else {
        $notailA = 1;
        print TESTA "*******$e1{'seq'}\n" if $debug;
    }
}
# Process for finding A sequence

# Note: The code snippet provided is incomplete and lacks context, so certain variables and conditions may not be explained accurately.


      # Output 'reads' results
if ($poly eq 'A|T') {
    if ($haveA and $haveT) {
        # If both A and T are present, keep A's result if lenA >= lenT
        if (length($seqT) < length($seqA)) {
            print OT "$e1{'id'}\n$seqT\n$e1{'pl'}\n$qualT\n";
            $cntFinalT++;
            $cntMissA++;
            if ($oraw) {
                print RAWT "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";
            }
        } else {
            print OA "$e1{'id'}\n$seqA\n$e1{'pl'}\n$qualA\n";
            if ($oraw) {
                print RAWA "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";
            }
            $cntFinalA++;
            $cntMissT++;
        }
    } elsif ($shortA or $shortT) {
        # If A or T is too short, don't output anything
        $cntShortA++;
        $cntShortT++;
    } elsif ($badTailT or $badTailA) {
        # If A or T has a bad tail, don't output anything
        $cntBadA++;
        $cntBadT++;
    } elsif ($haveA and !$haveT) {
        # If one side is missing, only output the corresponding side
        print OA "$e1{'id'}\n$seqA\n$e1{'pl'}\n$qualA\n";
        if ($oraw) {
            print RAWA "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";
        }
        $cntFinalA++;
        $cntNotailT++;
    } elsif ($haveT and !$haveA) {
        # If one side is missing, only output the corresponding side
        print OT "$e1{'id'}\n$seqT\n$e1{'pl'}\n$qualT\n";
        if ($oraw) {
            print RAWT "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";
        }
        $cntFinalT++;
        $cntNotailA++;
    } elsif (!$haveT and !$haveA) {
        $cntNotailA++;
        $cntNotailT++;
    }
} else {
    if ($haveA) {
        print OA "$e1{'id'}\n$seqA\n$e1{'pl'}\n$qualA\n";
        if ($oraw) {
            print RAWA "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";
        }
        $cntFinalA++;
    } elsif ($shortA) {
        $cntShortA++;
    } elsif ($badTailA) {
        $cntBadA++;
    } elsif ($notailA) {
        $cntNotailA++;
    }
    
    if ($haveT) {
        print OT "$e1{'id'}\n$seqT\n$e1{'pl'}\n$qualT\n";
        if ($oraw) {
            #print RAWT "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";
            print RAWT "$e1{'seq'}\n";
        }
        $cntFinalT++;
    } elsif ($shortT) {
        $cntShortT++;
        #print SHORTT "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n" if $debug;
        print SHORTT "$e1{'seq'}\n" if $debug;
    } elsif ($badTailT) {
        $cntBadT++;
        #print BEDT "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n" if $debug;
        print BEDT "$e1{'seq'}\n" if $debug;
    } elsif ($notailT) {
        # print NOT "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n" if $debug;
        print NOT "$e1{'seq'}\n" if $debug;
        $cntNotailT++;
    }
}
}#end while 
close($fh1);


if($findA){
   close(OA);
   if($debug){
   	close(TESTA);
   	}	
   	if($oraw){
   	close(RAWA);	
   	}
}


if($findT){
   close(OT);
   if($debug){
   	close(TESTT);
   	close(NOT);
   	close(BEDT);
   	close(SHORTT);
   	}	
   	if($oraw){
   	close(RAWT);	
   	}
}



print "\ntotal\t$cntTotal\n";
if ($findA) {
  print "[polyA]\nfinal\t$cntFinalA\nnotail\t$cntNotailA\ntooshort\t$cntShortA\nbadTail\t$cntBadA\nmissby(A|T)\t$cntMissA\n";
  if ($cntTotal!=$cntFinalA+$cntNotailA+$cntShortA+$cntMissA+$cntBadA) {
	print "cntTotal!=cntFinalA+cntNotailA+cntShortA+cntMissA+cntBadA\n";
  }
}
if ($findT) {
  print "[polyT]\nfinal\t$cntFinalT\nnotail\t$cntNotailT\ntooshort\t$cntShortT\nbadTail\t$cntBadT\nmissby(A|T)\t$cntMissT\n";
  print "region\t$md1\ndeep\t$md2\nreview\t$md3\n";
  if ($cntTotal!=$cntFinalT+$cntNotailT+$cntShortT+$cntMissT+$cntBadT) {
	print "cntTotal!=cntFinalT+cntNotailT+cntShortT+cntMissT+cntBadT\n";
  }
}


###############################################
sub readfq {
        my $fh = pop @_;
        my @entry = ();
        for (qw/id seq pl qual/) {
                my $line = readline($fh); chomp $line; 
                push @entry, $line; warn "line empty" if !$line;
        }
        return @entry;
}

