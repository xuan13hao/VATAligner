#!/usr/bin/perl -w
use strict;

my $ref_seq_file = shift;
my $clan_file = shift;

my $perfect = 0;
my $multi = 0;
my $wrong = 0;
my $miss = 0;

my %num_mapped;
my %correct_mapped;
my %checked_hash;

open my $RIN, "<$ref_seq_file" or die "Cannot open file: $!\n";
my @seq_len;
my @ref_pos;
my %id_hash;
my $n = 0;
while(<$RIN>) {
  chomp;
  my $line= $_;
  if($_ = /^>(\d+)\:(\d+)\-(\d+)/) {
    my @pos;
    push @pos, $1; push @pos, $2; push @pos, $3;
    push @ref_pos, \@pos;
    #print "$id\n";
    $line =~ s/^>//g;
    #print "$id\n";
    $id_hash{$line} = $n ++;
  } else  {
    push @seq_len, length($line);
    #print "$line\n";
  }
}
close $RIN;



#foreach(my $x = 0; $x < scalar(@seq_len); ++ $x)   {
#    print "$x   $seq_len[$x]\n";
#}
#die;

#foreach(keys %id_hash)   {
#    print "$_   $id_hash{$_}\n";
#}
#die;

my @content;
open my $CIN, "<$clan_file" or die "Cannot open file: $!\n";
my %score_hash;
my %count_hash;
#my $foo = <$CIN>;
while(<$CIN>) {
  next if(/^\#/);
  my @decom = split /\t/, $_;
  my $id = $id_hash{$decom[0]};
  if(exists $count_hash{$id} && exists $score_hash{$id} && $count_hash{$id} >= 1 && $decom[11] < $score_hash{$id})  {
    #print "$count_hash{$id}\t$score_hash{$id}\n";
    next;
  } else    {
    push @content, $_;
    $score_hash{$id} = $decom[11];
    $count_hash{$id} = 1;
  }
}
close $CIN;
  
for(my $i = 0; $i < scalar(@content); ++ $i)   {
  # record the read is mapped
  
  my @decom = split /\t+/, $content[$i];

  #print "$decom[0]\n";

  my $read = $id_hash{$decom[0]};
  $checked_hash{$read} = 1;
  $num_mapped{$read} ++;
  
  
  #print "read name: $read\n";
  #print "$map_first_count{$read}\n";
  #print "$map_second_count{$read}\n";
  
  # get the first ref info
  my $id1 = $ref_pos[$read][0];
  my $left1 = $ref_pos[$read][1]; 
  my $right1 = $ref_pos[$read][2];
  my $len1 = $right1 - $left1 + 1;
  
  #print "Check ref info:    $id1\t$left1\t$right1\t$id2\t$left2\t$right2\n";
  
  # check if this is an overlap
  my $tid = $decom[1];
  if($id1 eq $tid)  {
    my $l = $decom[8]; my $r = $decom[9];
    #print "Matching first arm:\t$tid\t$l\t$r\n";
    my $l_bound = $left1 > $l ? $left1 : $l;
    my $r_bound = $right1 < $r ? $right1 : $r;
    #print "First arm bounds 1:  $l_bound  $r_bound  $len1\n";
    if($l_bound < $r_bound && ($r_bound - $l_bound + 1) / $len1 >= 0.8)  {
      $correct_mapped{$read} += 1;
    }
  }
} 

for(my $id = 0; $id < scalar(@ref_pos); ++ $id) {
  if(!exists $checked_hash{$id})  {
    ++ $miss;
    #print "both miss: $id\n" if($map_counts <= 2);
  } elsif(exists $correct_mapped{$id}) {
    if ($num_mapped{$id} == 1)   {
      ++ $perfect;
      #print "$id\n";
      #print "perfect: $id\n" if($map_counts <= 2);
    } elsif($num_mapped{$id} > 1)   {
      ++ $multi;
    } 
  } elsif(exists $num_mapped{$id} && !exists $correct_mapped{$id}) {
    ++ $wrong;
    #print "partial wrong: $id\n" if($map_counts <= 2);
  } 
}
close $RIN;

print ">>>Read-level statistics:\n";
print "perfect\t$perfect\n";
print "multi\t$multi\n";
print "wrong\t$wrong\n";
print "miss\t$miss\n";
















