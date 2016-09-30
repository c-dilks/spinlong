#!/usr/bin/perl -w

my $tight="false";
if($#ARGV>=0) { $tight="true"; }

# creatively get number of classes
my $tmpn="tmpscript.C";
open(TMP,"> $tmpn");
print(TMP "{\n");
print(TMP "gSystem->Load(\"src/RunInfo.so\");\n");
print(TMP "RunInfo * RD = new RunInfo();\n");
print(TMP "if(!(RD->env->success)) return;\n");
print(TMP "EventClass * ev = new EventClass(RD->env,true);\n");
print(TMP "printf(\"%d\\n\",ev->N);\n");
print(TMP "}\n");
close(TMP);
my $nclass = `root -b -q $tmpn | tail -n1`;
chomp($nclass);
system("rm $tmpn");
print "NCLASSES=$nclass\n";

# execute add_diag.C on each class, one at a time
for (my $c=0; $c<$nclass; $c++) {
  print "-"x50 . "\nNOW EXECUTING SCRIPT ON CLASS $c\n" . "-"x50 . "\n";
  system("root -b -q print_diag.C'(${tight},${c})'");
}
