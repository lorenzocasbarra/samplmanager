#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
use Time::HiRes qw( time );
use Time::HiRes qw( gettimeofday );
use Time::HiRes qw( tv_interval );
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

my $input;
my @lines;
my @hosts;
my @guests;
my $start_time = time();

GetOptions (
    'input|i=s' => \$input, ##project list to input
) or die "No projects in input!\n";
open  (OUT,'>',$input."_times");
open (IN,$input);
chomp( @lines = <IN>);
close IN;
foreach my $project (@lines) {
    open (HOST, "$project/host.txt");
    chomp(@hosts = <HOST>);
    close HOST;
    open (GUEST, "$project/guest.txt");
    chomp(@guests = <GUEST>);
    close GUEST;
    say "PROJECT:\t$project";
    say OUT "PROJECT:\t$project";
    foreach my $host (@hosts) {
        say "\tHOST:\t$host";
        opendir(DIR,$project);
        my @list = readdir(DIR);
        closedir(DIR);
        foreach my $file (@list) {
            next if (
                $file !~ /.pdb$/ ||
                grep(/$file/,@hosts) ||
                !(grep(/$file/,@guests))
            );
            say "\tGUEST:\t".$file;
            my $ligand = $file;
            my $receptor = $host;
            $ligand =~ s/\.pdb//;
            $receptor =~ s/\.pdb//;
            my $guest_start = time();
            system("cd $project; docking.bash -l $ligand -r $receptor;");
            my $guest_end = time();
            say OUT "$receptor-$ligand\t".($guest_end - $guest_start);
        }
    }
}
my $end_time = time();
say OUT "END\t".($end_time - $begin_time);
close OUT;
