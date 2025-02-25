# action:
#     count the number of unique insert endpoints in each bin
# input:
#     stream of distinct insert endpoints in format "start1\tstrand_flag_1[,strand_flag_2]"
# outputs:
#     unique endpoint counts for each possible bin on chromosome
#     two reads starting at the same position in the same orientation are deduplicated
#     two reads starting at the same position in opposite orientations are each counted

use strict;
use warnings;

# constants
use constant { 
    START1 => 0,
    STRAND_FLAGS => 1
};

# variables
my $BIN_SIZE = $ENV{BIN_SIZE};
my @counts = (0) x ($ENV{MAX_BIN_I0} + 1);

# run the alignment data
while(<STDIN>){ 
    chomp; 
    my @f = split("\t"); 
    foreach my $strandFlag(split(",", $f[STRAND_FLAGS])){ 
        $counts[int(($f[START1] - 1) / $BIN_SIZE)]++;
    }
}

# print the counts
foreach my $binI0(0..$ENV{MAX_BIN_I0}){
    print $counts[$binI0], "\n";
}
