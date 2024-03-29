#!/usr/bin/perl -w

# Process the SAM file (comparison result of the tail is removed)
# Obtain the PAT, determine direction and coordinates based on provided poly=A or T
# Consider soft clipping (S mark) when calculating coordinates 6
# Allow judging alignment sequence rationality based on S and number of M (when S is too large or M is too small, poorly aligned sequence is considered)

use strict;
use Getopt::Long;
use File::Basename;
# use lib '/media/bmi/8E42CB1742CB0347/APAdb/STAR';  # You can set your library path
require("funclib.pl");

my %opts = ();
GetOptions(\%opts, "h", "sam=s", "poly=s", "m:i", "s:i", "more:s", "ofile:s");
my $USAGE = << "_USAGE_";

Usage:

  Only determine based on flag, without checking M and S positions,
  Coordinates are determined solely based on sequence length.
  MAP_parseSAM2PAT.pl -sam oxt.sam

  Allow further filtering using cigar M and S positions.
  MAP_parseSAM2PAT.pl -sam oxt.sam -m 30 -s 10

  Output columns MAPQ CIGAR CIGAR_M CIGAR_S XS,
  Coordinates are also determined based on M and S positions.
  MAP_parseSAM2PAT.pl -sam xx.sam -more T -ofile xx.sam.PAT

-h = help
-sam = SAM file (mapping result from the trimmed file)
-poly = A/T (SAM is from the ^T or A\$ file)
-ofile = output PA file (chr, strand, coord) (default is xx.PAT)
-m = 0 (default); cigar match (M + D + I) (>= M)
-s = 999 (default); cigar soft clip (S) (<= S)
-more = T/F (default) output more columns (MAPQ CIGAR CIGAR_M CIGAR_S XS)

_USAGE_

die $USAGE if $opts{'h'} || !$opts{'sam'} || !$opts{'poly'};
my $samfile = $opts{'sam'};
my $ofile = $opts{'ofile'};
my $poly = uc($opts{'poly'});
die "poly = A/T" if $poly ne 'A' and $poly ne 'T';
if (!$ofile) {
  $ofile = $samfile;
  if ($ofile =~ /sam$/i) {
    $ofile =~ s/sam$/PAT/i;
  } else {
    $ofile .= '.PAT';
  }
}
my $more = 0;
$more = 1 if $opts{'more'} eq 'T';

my $cigarM = '';
$cigarM = 1 if $opts{'m'} and $opts{'m'} > 0;
my $cigarS = '';
$cigarS = 1 if $opts{'s'} and ($opts{'s'} > 0 or $opts{'s'} eq '0');
my $minM = 0;
my $maxS = 999;
$minM = $opts{'m'} if $cigarM;
$maxS = $opts{'s'} if $cigarS;

die "samfile does not exist!" if !(-e $samfile);
die "sam cannot be the same as ofile" if $samfile eq $ofile;

print "input SAM file (sam) = $samfile\noutput file = $ofile\npoly = $poly\nM >= $minM\nS <= $maxS\nmore output = $more\n";

open(IN, "<$samfile");
my $line;
my $nskip = 0;
while ($line = <IN>) {
  if ($line =~ /^@/) {
    $nskip++;
    next;
  } else {
    last;
  }
}
close(IN);

open(IN, "<$samfile");

while ($nskip) {
  $line = <IN>;
  $nskip--;
}

my ($flagFld, $chrFld, $posFld, $cigarFld, $seqFld, $qFld) = (1, 2, 3, 5, 9, 4);

# MCIC-SOLEXA:2:29:1299:707#0/2  16  NC_003070  1076  60  24M1S  *  0  0  CCCCACACCCCCCCCACCCCCCAAC  #########################  NH:i:1
# MCIC-SOLEXA:2:54:719:1978#0/2  0  NC_003070  1077  60  24M2S  *  0  0  CCCACCCCCCCCCCACCCCCCAAACA  ##########################  NH:i:1

my (@items, $flag, $strand, $chr, $pos, $PATpos);
my ($flag0cnt, $flag1cnt) = (0, 0);
my ($cigarcnt, $noflagcnt) = (0, 0);
$chr = $strand = '';

open(OO, ">$ofile");

if ($more or $cigarM or $cigarS) { # 2015/1/3 If cigar needs to be checked
  while ($line = <IN>) {
    @items = split(/\t/, $line);
    $flag = $items[$flagFld];
    #$chr = $chrs{$items[$chrFld]};
    $chr = $items[$chrFld];  # 2014/12/26 If chr from SAM can be used directly
    $pos = $items[$posFld];
    my $mapq = $items[$qFld];

    my $cigar = $items[$cigarFld];
    my ($leftS, $referLen, $rightS) = parseCigar($cigar);

    if ($leftS + $rightS > $maxS or $referLen < $minM) {
      $cigarcnt++;
      next;
    }

    if ($flag == 0 and $poly eq 'T') {
      $PATpos = $pos - $leftS - 1;
      $strand = '-';
      $flag0cnt++;
    } elsif ($flag == 0 and $poly eq 'A') {
      $strand = '+';
      $PATpos = $pos + $rightS + $referLen;
      $flag0cnt++;
    } elsif ($flag == 16 and $poly eq 'A') {
      $strand = '-';
      $PATpos = $pos - $leftS - 1;
      $flag1cnt++;
    } elsif ($flag == 16 and $poly eq 'T') {
      $strand = '+';
      $PATpos = $pos + $rightS + $referLen;
      $flag1cnt++;
    } else {
      $noflagcnt++;
      next;
    }
    if ($PATpos <= 0) {
      $PATpos = 1;
    }
    if ($more) {
      my $XS = 0;
      $XS = 1 if $line =~ /XS:i:/;
      my $s = $leftS + $rightS;
      print OO "$chr\t$strand\t$PATpos\t$mapq\t$cigar\t$referLen\t$s\t$XS\n";
    } else {
      print OO "$chr\t$strand\t$PATpos\n";
    }
    next;
  }
} else { # Only based on flag
  while ($line = <IN>) {
    @items = split(/\t/, $line);
    $flag = $items[$flagFld];
    #$chr = $chrs{$items[$chrFld]};
    $chr = $items[$chrFld];  # 2014/12/26 If chr from SAM can be used directly
    $pos = $items[$posFld];
    if (($flag == 0 and $poly eq 'A') or ($flag == 16 and $poly eq 'T')) {
      $flag0cnt++;
      $strand = '+';
      $PATpos = $pos + length($items[$seqFld]);
    } elsif (($flag == 16 and $poly eq 'A') or ($flag == 0 and $poly eq 'T')) { # flag == 16
      $flag1cnt++;
      $strand = '-';
      $PATpos = $pos - 1;
    } else { # No flag
      $noflagcnt++;
      next;
    }
    if ($PATpos <= 0) {
      $PATpos = 1;
    }
    print OO "$chr\t$strand\t$PATpos\n";
    next;
  }
}

close(IN);
close(OO);
if ($cigarM or $cigarS) {
  print "M < $minM or S > $maxS\t$cigarcnt\n";
}
print "flag = *\t$noflagcnt\nflag_X\t$flag0cnt\nflag_Y\t$flag1cnt\n";

sub parseCigar {
  my $cigar = shift;
  my $leftS = 0;
  my $rightS = 0;
  my $referLen = 0;
  my @len = $cigar =~ /(^\d+(?=S))/g;
  $leftS = $len[0] if @len;
  @len = $cigar =~ /(\d+(?=S$))/g;
  $rightS = $len[0] if @len;
  my @nums = $cigar =~ /\d+(?=[M D])/g;
  for my $i (0..$#nums) {
    $referLen += $nums[$i];
  }
  return ($leftS, $referLen, $rightS);
}
