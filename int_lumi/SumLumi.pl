#!/usr/bin/perl -w
use Data::Dumper;

# get year
my $year;
if($#ARGV==0) {
  $year = $ARGV[0];
} else {
  print "Usage: $0 [year]\n";
  exit;
}


# run list files (copied to this directory)
my $goodlist = "goodruns_${year}.dat"; # from goodruns.dat
my $pi0list = "analysed_runlist_${year}.list"; # from Asym.C-->runlist.list


# column of lists which is run number
my $goodcol = 0;
my $pi0col = 1;


# build runlists
my @goodruns;
my @pi0runs;
open(GOODFILE,$goodlist) or die("ERROR: ${goodlist} not found\n");
open(PI0FILE,$pi0list) or die("ERROR: ${pi0list} not found\n");
foreach $line (<GOODFILE>) {
  chomp($line);
  my @cols = split " ", $line;
  push(@goodruns,$cols[$goodcol]);
}
foreach $line (<PI0FILE>) {
  chomp($line);
  my @cols = split " ", $line;
  push(@pi0runs,$cols[$pi0col]);
}
close(GOODFILE);
close(PI0FILE);
#print Dumper(\@goodruns);
#print Dumper(\@pi0runs);


# loop through luminosity data files in lumi directory
my @lumifiles = <lumi_${year}/*.txt>;
my $trg;
my @trgname;
my @total_sum;
my @good_sum;
my @pi0_sum;
my @total_sumps;
my @good_sumps;
my @pi0_sumps;
my $idx=0;
my $lumin;
my $luminps;
my $prescale;
my $runnum;
foreach $lumifile (@lumifiles) {
  $trg = $lumifile;
  $trg =~ s/lumi_${year}\/lum_perrun_//;
  $trg =~ s/\.txt//;
  #print "idx=$idx  trg=$trg\n";

  $trgname[$idx] = $trg;
  $total_sum[$idx] = 0;
  $good_sum[$idx] = 0;
  $pi0_sum[$idx] = 0;
  $total_sumps[$idx] = 0;
  $good_sumps[$idx] = 0;
  $pi0_sumps[$idx] = 0;

  open(LUMI,$lumifile) or die("ERROR: ${lumifile} not found\n");
  foreach $line (<LUMI>) {
    chomp($line);
    my @cols = split " ", $line;
    $runnum = $cols[0];
    $lumin = $cols[4];
    $prescale = $cols[5];
    $luminps = $lumin * $prescale;

    $total_sum[$idx] += $lumin;
    $total_sumps[$idx] += $luminps;
    if(grep $_ eq ${runnum}, @goodruns) { 
      $good_sum[$idx] += $lumin;
      $good_sumps[$idx] += $luminps;
    }
    if(grep $_ eq ${runnum}, @pi0runs) {
      $pi0_sum[$idx] += $lumin;
      $pi0_sumps[$idx] += $luminps;
    }
  }
  $idx++;
}
my $n_trgs = $idx-1;


# print output
print("\n\nRUN ${year} LUMINOSITY TOTALS\n\ntriggered luminosity (sum lumi [pb^-1])\n");
print "-" x 100;
print("\n");
print "trigger   lumi_all_runs   lumi_good_runs   lumi_pi0_runs\n";
for my $i (0 .. $n_trgs) {
  print("$trgname[$i]   $total_sum[$i]   $good_sum[$i]   $pi0_sum[$i]\n");
}
print "-" x 100;
print("\n\nde-prescale-ified total luminosity (sum lumi*prescale [pb^-1])\n");
print "-" x 100;
print("\n");
print "trigger   lumi_all_runs   lumi_good_runs   lumi_pi0_runs\n";
for my $i (0 .. $n_trgs) {
  print("$trgname[$i]   $total_sumps[$i]   $good_sumps[$i]   $pi0_sumps[$i]\n");
}
print "-" x 100;
print("\n");
