#!/usr/bin/perl -w
use strict;

my $ref_seq_file = shift;
my $clan_file = shift;

my $perfect = 0;
my $partial_multi = 0;
my $both_multi = 0;
my $partial_wrong = 0;
my $both_wrong = 0;
my $partial_miss = 0;
my $both_miss = 0;
my $correct = 0;
my $unique = 0;

my %checked_hash;
my %first_mapped;
my %second_mapped;
my %first_wrong;
my %second_wrong;
my %map_first_count;
my %map_second_count;
my %mapped_first_locs;
my %mapped_second_locs;

open my $RIN, "<$ref_seq_file" or die "Cannot open file: $!\n";
my @seq_len;
my @ref_pos;
my %id_hash;
my $n = 0;
while(<$RIN>) {
  chomp;
  my $line= $_;
  if($_ = /^>(\d+)\:(\d+)\-(\d+)\|\|(\d+)\:(\d+)\-(\d+)/) {
    my @pos;
    push @pos, $1; push @pos, $2; push @pos, $3; push @pos, $4; push @pos, $5; push @pos, $6;
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
  my @decom = split /\t/, $_;
  if(exists $count_hash{$decom[0]} && $count_hash{$decom[0]} >= 2)  {
    next;
  } else    {
    #print "$score_hash{$decom[0]}\t$count_hash{$decom[0]}\n";
    push @content, $_;
    if(exists $count_hash{$decom[0]})   {
        if($decom[11] < $score_hash{$decom[0]})    {
            $score_hash{$decom[0]} = $decom[11];
            ++ $count_hash{$decom[0]};
        }   
    }   else    {
        $score_hash{$decom[0]} = $decom[11];
        $count_hash{$decom[0]} = 1;
    }
  }
}
close $CIN;
  
for(my $i = 0; $i < scalar(@content); ++ $i)   {
  # record the read is mapped
  
  my @decom = split /\t+/, $content[$i];

  my $read = $id_hash{$decom[0]};
  $checked_hash{$read} = 1;
  my $mid = ($decom[6] + $decom[7]) / 2; 
  
  #print "map location:  $mid    $seq_len[$read]\n";
  $map_first_count{$read} ++ if ($mid <= $seq_len[$read] / 2);
  $map_second_count{$read} ++ if ($mid > $seq_len[$read] / 2);
  
  #print "read name: $read\n";
  #print "$map_first_count{$read}\n";
  #print "$map_second_count{$read}\n";
  
  # get the first ref info
  my $id1 = $ref_pos[$read][0];
  my $left1 = $ref_pos[$read][1]; 
  my $right1 = $ref_pos[$read][2];
  my $len1 = $right1 - $left1 + 1;
  
  # get the second ref info
  my $id2 = $ref_pos[$read][3];
  my $left2 = $ref_pos[$read][4];; 
  my $right2 = $ref_pos[$read][5];;
  my $len2 = $right2 - $left2 + 1;
  
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
      $first_mapped{$read} += 1;
    } else  {
      $first_wrong{$read} += 1;
    }
  } elsif($id2 eq $tid) {
    my $l = $decom[8]; my $r = $decom[9];
    #print "Matching second arm: $tid\t$l\t$r\n";
    my $l_bound = $left2 > $l ? $left2 : $l;
    my $r_bound = $right2 < $r ? $right2 : $r;
    #print "Second arm bounds 2:  $l_bound  $r_bound  $len2\n";
    if($l_bound < $r_bound && ($r_bound - $l_bound + 1) / $len2 >= 0.8)  {
      $second_mapped{$read} = +1;
    } else  {
      $second_wrong{$read} = +1
    }
  } else  {
    #print "$mid   $decom[4]\n";
    $first_wrong{$read} += 1 if ($mid <= $seq_len[$read] / 2);
    $second_wrong{$read} += 1 if ($mid > $seq_len[$read] / 2);
    #print "$second_wrong{$read}\n";
  }
} 

for(my $id = 0; $id < scalar(@ref_pos); ++ $id) {
  my $local_correct = 0;
  my $local_unique = 0;
  if(!exists $checked_hash{$id})  {
    ++ $both_miss;
    #print "both miss: $id\n" if($map_counts <= 2);
  } elsif(exists $first_mapped{$id} && exists $second_mapped{$id}) {
    if (!exists $first_wrong{$id} && !exists $second_wrong{$id})   {
      ++ $perfect;
      $correct += 2;
      $local_correct += 2;
      #print "$id\n";
      #print "perfect: $id\n" if($map_counts <= 2);
    } elsif(($second_mapped{$id} && !exists $second_wrong{$id}) || ($first_mapped{$id} && !exists $first_wrong{$id}))   {
      ++ $partial_multi;
      $correct += 1;
      $local_correct += 1;
      #print "partial multi: $id\n" if($map_counts <= 2);
    } else    {
      ++ $both_multi;
      #print "both multi: $id\n" if($map_counts <= 2);
    }
  } elsif((exists $first_mapped{$id} && !exists $first_wrong{$id} && !exists $second_mapped{$id} && exists $second_wrong{$id}) || 
      (!exists $first_mapped{$id} && exists $first_wrong{$id} && exists $second_mapped{$id} && !exists $second_wrong{$id})
  ) {
    ++ $partial_wrong;
    $correct += 1;
    $local_correct += 1;
    #print "partial wrong: $id\n" if($map_counts <= 2);
  } elsif((exists $first_mapped{$id} && !exists $first_wrong{$id} && !exists $second_mapped{$id} && !exists $second_wrong{$id}) || 
      (!exists $first_mapped{$id} && !exists $first_wrong{$id} && exists $second_mapped{$id} && !exists $second_wrong{$id})) {
    ++ $partial_miss;
    $correct += 1;
    $local_correct += 1;
    #print "partial miss: $id\n" if($map_counts <= 2);
  } elsif(!exists $first_mapped{$id} && exists $first_wrong{$id} && !exists $second_mapped{$id} && exists $second_wrong{$id}) {
    ++ $both_wrong;
    #print "both wrong: $id\n" if($map_counts <= 2);
    #print "$id\n";
  } else  {
      # one hit is wrong, the othe hit is multiple
      #print "========\n";
      #print "$id\n";
      #print "first mapped:  $first_mapped{$id}\n" if (exists $first_mapped{$id});
      #print "second mapped: $second_mapped{$id}\n" if (exists $second_mapped{$id});
      #print "first wrong\n" if (exists $first_wrong{$id});
      #print "second wrong\n" if (exists $second_wrong{$id});
      ++ $both_wrong;
      #print "multi miss: $id\n" if($map_counts <= 2);
  }
  ++ $unique if(exists $map_first_count{$id} && $map_first_count{$id} == 1);
  ++ $unique if(exists $map_second_count{$id} && $map_second_count{$id} == 1);
  ++ $local_unique if(exists $map_first_count{$id} && $map_first_count{$id} == 1);
  ++ $local_unique if(exists $map_second_count{$id} && $map_second_count{$id} == 1);
  
  if($local_correct > $local_unique)  {
    #print "$local_correct\t$local_unique\n";
    #print "$map_first_count{$id}\t$map_second_count{$id}\n";
    #print "$id  @{$ref_pos[$id]}\n";
    #die;
  }
}
close $RIN;

print ">>>Read-level statistics:\n";
print "perfect\t$perfect\n";
print "partial_multi\t$partial_multi\n";
print "both_multi\t$both_multi\n";
print "partial_wrong\t$partial_wrong\n";
print "both_wrong\t$both_wrong\n";
print "partial_miss\t$partial_miss\n";
print "both_miss\t$both_miss\n";
print ">>>Strand-level statistics:\n";
print "correct\t$correct\n";
print "unique\t$unique\n";
















