#!/usr/bin/perl

use strict;

use 5.010;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use File::Basename;
use List::Util qw( min max );



our $verbose;
our $output;
our $debug;
our $evaluate;
our $prepare;
our $map;
our $host;
our $plot;
our $template;
our $consent;
my $input_ref;
my @inputs;
my @selects;
my @flags;
my %references;
my %info;
my %mue;
my %tau;
my %pea;
my %angular;
my %intercept;
my %lin;

GetOptions(
    'input|i=s' => \@inputs,    #input folder or file.txt to parse
    'ref|r=s' => \$input_ref,   #input csv for reference values
    'output|o=s' => \$output,   #output file
    'flags|f=s' => \@flags,     #info to retrieve and write out
    'select|s=s' => \@selects,  #info to filter
    'verbose|v' => \$verbose,   #verbose mode on
    'evaluate|e' => \$evaluate, #evaluate mode one, perform RMSD,KENDALL,PEARSON, LIN
    'prepare|p' => \$prepare,   #enables preparation of a sampl submission, based on all .dlg found in directory input
    'map|m' => \$map,   #sorts an evaluated exp file by lin coefficient
    'template|t=s' => \$template,    #template for plot in .gplt format 
    'plot' => \$plot,   #plot a sorted exp with gnuplot
    'debug|d' => \$debug,
) or die "Usage: $0 --from NAME\n";

if ($prepare) {
    goto PREPARE;
}

if ($map) {
    goto MAP;
}



if ($output) {
    if (-e $output) {`rm $output`;}
}

if ($input_ref) {
    if ($debug){print "References loaded from $input_ref:\n"};
    %references = loadref($input_ref);
    if ($debug){print "\n"};
}

if ($plot) {

    goto PLOT;

}


if (!(grep(/Predictions/, @flags)) && $evaluate) {
    push @flags, "Predictions";
}
if ($flags[0] eq "all") {
    @flags = ("Predictions", "Participant name", "Participant organization", "Name","Software","Category","Ranked");
}


foreach my $input (@inputs) {
    if (-f $input) {
        our ($dirname) = dirname(dirname($input));
        our $file = $inputs[0];
        print("\tWorking on: $file\n") if ($debug || $verbose);
        `dos2unix --quiet $file`;
        open(IN,'<',$file);
        chomp( my @lines = <IN> );
        if (@selects) {
            my $consent = check(\@lines);
        } else {$consent = 1;}
        close IN;
        if ($consent == 1) {
            my %info = parse($file, \@flags);
            if ( (scalar %{$info{"Predictions"}} ) == 0) {
                print "\tWarning: $file does not contain $host\n";
                next;
            }
            if ($debug || $verbose) {hprint(\%info,\%references)}
            if ($evaluate){
                $mue{$file} = mue(\%references,\%info);
                $tau{$file} = tau(\%references,\%info);
                $pea{$file} = pea(\%references,\%info);
                $lin{$file} = lin(\%references,\%info,$pea{$file});
                ($angular{$file} ,$intercept{$file}) = intercept(\%references,\%info,$pea{$file});
            }
            if ($output) {
                prettyprint(\%info,\%references,\%mue,\%tau, \%pea, \%angular, \%intercept, \%lin);
            }
        }
    }

    if (-d $input ) {
        our ($dirname) = dirname($input);
        opendir(DIR, $input);
        foreach our $file (readdir(DIR)) {
            if ($file eq "." ||
                $file eq ".." ||
                $file !~ /.txt/ 
            ) {next}
            if ($verbose || $debug) {say "Working on $file"};
            open(IN,'<',"$input/$file");
            chomp( my @lines = <IN> );
            if (@selects) {
                my $consent = check(\@lines);
            } else {$consent = 1;}
            close IN;
            if ($consent == 1) {
                my %info = parse("$input/$file",\@flags);
                if ( (scalar %{$info{"Predictions"}} ) == 0) {
                    print "\tWarning: $file does not contain $host\n";
                    next;
                }
                if ($debug || $verbose) {hprint(\%info,\%references)}
                if ($evaluate){
                    $mue{$file} = mue(\%references,\%info);
                    $tau{$file} = tau(\%references,\%info);
                    $pea{$file} = pea(\%references,\%info);
                    $lin{$file} = lin(\%references,\%info,$pea{$file});
                    ($angular{$file} ,$intercept{$file}) = intercept(\%references,\%info,$pea{$file});
                }
                if ($output) {
                    prettyprint(\%info,\%references,\%mue,\%tau, \%pea, \%angular, \%intercept, \%lin);
                }
            } else {next}
        }
        closedir(DIR);   
    }
}

PREPARE:
if ($prepare) {
    my %energies;
    my %values;
    my @hosts;
    my @guests;
    my $dg0;
    my $dgs;
    my $dg1;
    my $d;
    my @files;
    foreach my $input (@inputs) {   
        if (-d $input ) {
            
            if ($debug) {say "Current input:$input";}
            opendir(DIR, $input);
            open (IN , "$input/host.txt") or die "Critical: no $input/host.txt\n";
            chomp( my @hosts = <IN> );
            close IN;
            open (IN , "$input/guest.txt") or die "Critical: no $input/guest.txt\n";
            chomp( my @guests = <IN> );
            close IN;
            #<STDIN>;
            open (OUT, '>', "$input/ourdock.txt");
            print OUT "Predictions:\n";
            @files = readdir(DIR);
            my @files_sorted = sort {$b <=> $a || $a cmp $b} @files;
            foreach our $file (@files_sorted) {
                if (-d $file ||
                    $file eq "." ||
                    $file eq ".." ||
                    $file !~ /.dlg/ 
                ) {next}
                say ">$file";
                my $id = $file;
                $id =~ s/.dlg//;
                my ($ligand , $receptor) = split (/_/, $id, 2);
                if ($debug) {say "Guest and Host are: $ligand\t$receptor";}
                if (!(grep(/^$receptor\.pdb/,@hosts)) || !(grep(/^$ligand\.pdb/,@guests))) {next};
                my @res = `grep Binding $input/$file | grep -e ^USER`;
                if (scalar(@res) > 1) {
                    my $i = 0;
                    if ($verbose || $debug) {say ">Values found:";}
                    foreach (@res) {

                        my $tmp = $_;
                        chomp $tmp;
                        $tmp =~ /=.+\s(.+?)\skcal\/mol/;
                        $values{"$ligand\@$receptor\@$i"} = $1;
                        if ($verbose || $debug) {print "\t$receptor-$ligand\t".$values{"$ligand\@$receptor\@$i"}."\n";}
                        $i++;

                    }
                    $dg0 = $values{"$ligand\@$receptor\@0"};
                    if ($verbose || $debug) {say ">Correction:\n\t\$dg0\t\$dg1\t\$dgs += exp(-(\$dg1 -\$dg0)/0.596)\t\$d = 0.6*log(1+\$dgs)\tnew dg";}
                    $dgs = 0;
                    for (my $j = 1; $j <= $i-1; $j++) {
                        $dg1 = $values{"$ligand\@$receptor\@$j"};
                        if ($debug) {print "\nDBG: $dgs + exp( -($dg1 - $dg0)/0.596) = ";}
                        $dgs += exp(-($dg1 -$dg0)/0.596);                         
                        if ($debug) {print "$dgs";}
                        if ($verbose || $debug) {say "\n\t$dg0\t$dg1\t$dgs\t\t";}
                    }
                    $d = 0.6*(log(1 + $dgs));
                    $energies{$ligand."@".$receptor} = $dg0 - $d;
                    if ($verbose || $debug) {say "\n\t$dg0\t\t$dgs\t$d\t".$energies{$ligand."@".$receptor}."\n#\n";}
                    if ($debug) {say "\nDBG: $dg0 - 0.6 * (log( 1 + $dgs )) = ".$energies{$ligand."@".$receptor};}

                    print OUT "$receptor-$ligand".",".$energies{"$ligand"."@"."$receptor"}."\n";
                } else {            
                    my $tmp = $res[0];
                    chomp $tmp;
                    $tmp =~ /=.+\s(.+?)\skcal\/mol/;
                    $energies{"$ligand\@$receptor"} = $1;
                    $values{"$ligand\@$receptor@0"} = $1;
                    if ($verbose || $debug) {print ">Value found:\n\t$receptor-$ligand\t".$energies{"$ligand\@$receptor"}."\n\n#\n";}
                    print OUT "$receptor-$ligand".",".$energies{"$ligand"."@"."$receptor"}."\n";
                }
            }
            closedir(DIR);
            print OUT "#\nName:\nautodock4_default\n#\nParticipant name:\nLorenzo Casbarra\n#\nParticipant organization:\nUniversity of Florence, Italy\n#\nSoftware:\nautodock4_default\n#\nMethod:\nautodock_default\n#\nCategory:\nDocking\n#\nRanked:\nNo\n#";
            close OUT;
            `cp $input/ourdock.txt $input/ANALYSIS`;
            if ($verbose || $debug) {say "My predictions should look like this:\n";}
            foreach my $k (sort {$b <=> $a || $a cmp $b} keys %energies) {
                
                my ($ligand,$receptor) = split /@/, $k;
                if ($verbose || $debug) {say "$receptor-$ligand\t".$energies{$k};}
                if ($debug) {say "Previous $receptor-$ligand\t".$values{"$ligand\@$receptor\@0"};}
            }
        } else {die "Warning: Input is not a directory"}
    }
}

MAP:
if ($map) {
    my %indent;
    my %hash;
    my %ref;
    my $mastercheck;
    my $masterkey;
    my $section;
    my $host;
    my $guest;
    my @hosts;
    my @guests;
    my @infokeys = ("Name","Category","Participant name","Participant organization","Software","Ranked","mue","kendall","pearson","angular","intercept","Lin");
    foreach my $input (@inputs) {
        if (-f $input) {
            print $input."\n";
            
            close IN;
            open (IN, $input);
            chomp( my @lines = <IN> );
            close IN;
            our $dirname = dirname($input);
            our $filename= fileparse($input);
            open (OUT, '>', "$dirname/MAPPED/"."$filename"."_mapped.tsv");
            for (my $i = 0; $i <= $#lines; $i++){
                LABEL:
                next if ($lines[$i] eq "#");
                
                if ($lines[$i] =~ /^>/) {$mastercheck = 0}
                if ($lines[$i] eq ">File") {
                    $section = $lines[$i];
                    $section =~ s/>//;
                    $masterkey = $lines[$i+1];
                    print $masterkey."\n";
                    $i++;
                    next;
                }
                if (
                    $lines[$i] eq ">Ranked" or
                    $lines[$i] eq ">Software" or
                    $lines[$i] eq ">Category" or
                    $lines[$i] eq ">Name" or
                    $lines[$i] eq ">Participant name" or
                    $lines[$i] eq ">Participant organization" or
                    $lines[$i] eq ">Predictions" or
                    $lines[$i] eq ">mue" or
                    $lines[$i] eq ">kendall" or
                    $lines[$i] eq ">pearson" or
                    $lines[$i] eq ">angular" or
                    $lines[$i] eq ">intercept" or
                    $lines[$i] eq ">Lin"
                ) {      
                    $mastercheck = 1;
                    $section = $lines[$i];
                    $section =~ s/>//;
                    $indent{$section} = uc($section);
                    if ($section eq "Predictions") {
                        $i += 2;
                        goto LABEL2;
                    }
                        
                    if ($lines[$i+1] =~ />/) {
                        $i++;
                        goto LABEL;
                    }
                    my $info = $lines[$i+1];
                    $hash{$section."@".$masterkey} = $info;
                                            
                    
                    $i++;
                    next;
                }
                    LABEL2:
                    if ($mastercheck == 1 && $section eq "Predictions") {
                        my %h1;
                        my %h2;
                        next if ($lines[$i] =~ /HG\tSUB\tREF/);
                        until ($lines[$i] =~ /^>/) {
                            my($a, $b, $c) = split(/\t/,$lines[$i]);
                            $ref{$a} = $c;
                            ($host , $guest) = split (/-/,$a);
                            $hash{$a."@".$masterkey} = $b;
                            if (!(grep(/^$host/,@hosts))) {
                                push @hosts, $host;
                            }
                            if (!(grep(/^$guest/,@guests))) {
                                push @guests, $guest;
                            }
                            $i++;
                        }
                        goto LABEL;
                    }                
            }
            my %tmp;
            print OUT "HOST:";
            for (my $i = 0; $i <= $#hosts; $i++) {
                next if ($tmp{$hosts[$i]});
                $tmp{$hosts[$i]} = 1;
                print OUT "\t".$hosts[$i];
            }
            print OUT "\n";
            print OUT "GUESTS:";
            for (my $i = 0; $i <= $#guests; $i++) {
                next if ($tmp{$guests[$i]});
                $tmp{$guests[$i]} = 1;
                print OUT "\t".$guests[$i];
            }
            print OUT "\n";
            
            print OUT "RANK\tFILE";
            foreach (@infokeys) {
                if ($indent{$_}) {
                    print OUT "\t".$indent{$_};
                }
            }
            
            if ($indent{"Predictions"}) {
                foreach my $host (@hosts) {
                    foreach my $guest (@guests) {
                        print OUT "\t".$host."-".$guest;
                    }
                }
            }                
            print OUT "\n";
            my $i=1;
            foreach my $key (sort { $hash{$b} <=> $hash{$a} } keys %hash) {
                if ($key =~ /Lin@/) {
                    my ($r,$file) = split /@/, $key;
                    print OUT "$i\t$file";
                    foreach (@infokeys) {
                        if ($indent{$_}) {         
                            print OUT "\t".$hash{$_."@".$file};
                        }
                    }
                    
                    if ($indent{"Predictions"}) { 
                        foreach my $host (@hosts) {
                            foreach my $guest (@guests) {
                                print OUT "\t".$hash{$host."-".$guest."@".$file};
                            }
                        }
                        print OUT "\n";
                    }
                    $i++;
                }
            }   
            close OUT;
        }
    }
}





PLOT:

if ($plot) {
    my @exps;
    my @splitted;
    my @ourvalues;
    my @bestvalues;
    my $intplace;
    my $angplace;
    my $j=0;
    my %references = loadref($input_ref);
    foreach my $input (@inputs) {
        if (-f $input) {
            open (IN,$input);
            chomp( my @lines = <IN> );
            close IN;
            my @hosts = split /\t/,$lines[0];
            my @guests = split /\t/,$lines[1];
            $input =~ s/\.tsv//;
            my @tmp = split /\t/, $lines[2];
            my $hg_start = (scalar @tmp) - (scalar @hosts - 1) * (scalar @guests - 1);
            for (my $j = 0; $j <= $#tmp; $j++) {
                if ($tmp[$j] eq "INTERCEPT") {$intplace = $j;}
                if ($tmp[$j] eq "ANGULAR") {$angplace = $j;}    
            }
            my @hgs;
            for (my $i = $hg_start; $i <= $#tmp; $i++) {
                push @hgs , $tmp[$i];
            }
            foreach my $hg (@hgs) {
                push @exps, $references{$hg};  #array with experimental values ordered
                print "$hg\t".$references{$hg}."\n";
            }
            
            my @best = split /\t/, $lines[3];   #array with best score
            for (my $i = 4; $i <= $#lines; $i++) {                
                @splitted = split /\t/, $lines[$i];  #array with our values
                if ($splitted[1] eq "ourdock.txt") {
                    last;
                }
            }
            my $ourname = $splitted[1];
            $ourname =~ s/\.txt//;
            my $bestname = $best[1];
            $bestname =~ s/\.txt//;
            my $bestint = $best[$intplace];
            my $bestang = $best[$angplace]; #line will look like this ($bestang * x + $bestint)
            
            my $ourint = $splitted[$intplace];
            my $ourang = $splitted[$angplace];
            for (my $i = $hg_start; $i <= $#splitted; $i++) {
                push @ourvalues, $splitted[$i];
                push @bestvalues, $best[$i];
            }
            my $b_min = min(@bestvalues);
            my $e_min = min(@exps);
            my $o_min = min(@ourvalues);
            my $b_max = max(@bestvalues);
            my $e_max = max(@exps);
            my $o_max = max(@ourvalues);
            my $max = max($b_max,$e_max,$o_max) + 0.5;
            my $min = min($b_min,$e_min,$o_min) - 0.5;
            my $bestout = dirname($input)."/".fileparse($input)."_best";
            my $ourout = dirname($input)."/".fileparse($input)."_our";
            open (OUT, '>', $bestout);
            open (OUT2, '>', $ourout);
            for (my $i = 0; $i <= $#bestvalues; $i++) {                
                print OUT $hgs[$i]."\t".$bestvalues[$i]."\t".$exps[$i]."\n";
                print OUT2 $hgs[$i]."\t".$ourvalues[$i]."\t".$exps[$i]."\n";
            }
            close OUT;
            close OUT2;
            my $outgplt = dirname($input)."/".fileparse($input).".gplt";
            my $title = fileparse($input);
            open(OUT, '>', $outgplt);
            open(IN,"$template");
            my @lines3 = <IN>;
            close IN;
            for (my $i = 0; $i <= $#lines3; $i++) {
                print $lines3[$i]."\n";                
                $lines3[$i] =~ s/NAME/$outgplt/g;
                $lines3[$i] =~ s/BEST/$bestout/g;
                $lines3[$i] =~ s/TITLE/$title/g;
                $lines3[$i] =~ s/OUR/$ourout/g;                    
                $lines3[$i] =~ s/XMI/$min/g;
                $lines3[$i] =~ s/XMA/$max/g;
                $lines3[$i] =~ s/AYX/$bestang/g;
                $lines3[$i] =~ s/BYX/$bestint/g;
                $lines3[$i] =~ s/CYX/$ourang/g;
                $lines3[$i] =~ s/DYX/$ourint/g;
                $lines3[$i] =~ s/Best/$bestname/g;
                $lines3[$i] =~ s/Our/$ourname/g;
                print OUT $lines3[$i];
            }
            close OUT;
            `gnuplot -c $outgplt`;
            
        }
    }





}

##SUBROUTINES

sub hprint {
    my $h = shift @_;
    my %h = %$h;
    my $r = shift @_;
    my %r = %$r;
    our $file;
    
    foreach my $k (sort keys %h) {
        if ($k eq "Software") {
            print "\t$k\t".$h{$k}."\n";
            goto LABEL3;
        }
        if ($k eq "Predictions") {
            say "\t$k\t".%{$h{$k}}."/".%r;
            say "\tHG\tREF\tSUB";
            foreach my $k2 (sort keys %{$h{$k}}) {
                say "\t".$k2."\t".$r{$k2}."\t".%{$h{$k}}{$k2};
            }
        } else {
            print "\t$k:\t";
            foreach my $p (@{$h{$k}}) {
                print $p;
            }       
        }
        LABEL3:
        print "\n";
    }
    print"\n";
}






#Makes hash from a reference values csv 
sub loadref {
    my $csv = shift @_; #experimental_values.csv
    my $l;
    my @array;
    my %hash;
    
    open(IN,$csv) or die "Error:  No $csv\n";
    while ($l = <IN>) {
        chomp $l;
        $l =~ s/NaN/0/g;
        my @array = split ';', $l;
        next if ($l =~ /ID/);
        my $n = scalar(@array);
        $hash{$array[0]} = $array[($n-2)];
        if ($debug){print $array[0]."\t".$array[($n-2)]."\n"};
    }
    close IN;
    return %hash;
}


sub parse {
    my $i = shift @_;
    my $f = shift @_;
    my @f = @$f;
    my %h;
    my %p;
    my @a;
    our $verbose;
    
    foreach my $g (@f) {
        if ($g =~ /Predict/) {
            %p = getprediction($i);
            $h{$g} = \%p;
        } else {
            my @a = getarray($i,$g);
            if ($g =~ /Software/) {
                 my $s = join(";",@a);
                 $h{$g} = $s;   
            } else {            
                $h{$g} = \@a;
            }
        }        
    }
    #foreach my $k (sort keys %{$hash{"Predictions"}}) {print %{$hash{"Predictions"}}{$k}."\n"}
    return %h;
}


sub getprediction {
    my $i = shift @_;
    my $g = "Predictions:";
    my %h;
    my $c = 0;
    our $dirname;
    open (IN , "$dirname/host.txt") or die "host.txt NOT FOUND in $dirname";
    chomp( our @hosts = <IN> );
    close IN;
    open (IN , "$dirname/guest.txt") or die "guest.txt NOT FOUND in $dirname";
    chomp( our @guests = <IN> );
    close IN;
    open(IN,'<',$i);
    chomp( my @ls = <IN> );
    close IN;
    if ($g =~ /Prediction/) {
        foreach my $l (@ls) {
            $l =~ s/\x0d{0,1}\x0a\Z//s;
            $l =~ s/\t//g;
            $l =~ s/\s//g;
            $l =~ s/kcal\/mol//g;
            if ($l =~ $g) {
                $c = 1;
                next;
            }
            if ($l eq "#" && $c == 1) {
                last
            }
            if ( $l !~ /^#$/ && $c == 1) {
                next if ($l =~ /^#/ || $l eq "");
                last if ($l eq "#" && $c == 1);
                $l =~ s/\s//;
                my @a = split(',' , $l);
                my $k = shift(@a);
                my ($x,$y) = split(/-/,$k);
                if ((grep(/^$x\.pdb/,@hosts)) && (grep(/^$y\.pdb/,@guests))) {
                    $h{$k} = $a[0];
                }
            }
        }
        return %h;
    }
}

#Generates an array for infos
sub getarray {
    my $i = shift @_;    #file to open as array
    my $g = shift @_;    #info to capture
    my $c = 0;
    my @a;
    
    open(IN,'<',$i);
    chomp( my @ls = <IN> );
    close IN;
    foreach my $l (@ls) {
        $l =~ s/\x0d{0,1}\x0a\Z//s;
        $l =~ s/\t//g;
        if ($l =~ $g) {
            $c = 1;
            next;
        }
        if (($l eq "#" || $l eq "" || $l eq "\n") && $c == 1) {
            last
        }
        if ($l !~ /^#$/ && $c == 1) {
            next if ($l =~ /^#/ || $l eq "" || $l eq "\n");
            last if ($l eq "#" && $c == 1);
            push @a, $l;
        }
    }
    return @a;
}

#RMSD calculation
sub mue { 
    my $r = shift @_;   #reference values hash
    my %r = %$r;
    my $h = shift @_;   #submission values hash
    my %h = %$h;     
    my $s = 0;
    my $m;
    our $verbose;
    our $debug;    
    foreach my $k ( sort keys %{$h{"Predictions"}} ) {
        if ($debug) {print "$k\t$s + sqrt(($r{$k} - ".%{$h{"Predictions"}}{$k}.")^2) = "}
        $s += sqrt(($r{$k} - %{$h{"Predictions"}}{$k})**2);
        if ($debug){print "$s\n";}
    }
    my $n = scalar( keys %{$h{"Predictions"}});
    $m = $s / $n;
    if ($debug){say "Number of Predictions: $n\nMUE: $s / $n = $m"};
    if ($verbose || $debug) {print"\tMUE:\t$m\n\tMue evaluation is done\n\n"};
    return $m;
}

#Kendall test
sub tau {
    my $r = shift @_;   #reference values hash
    my %r = %$r;
    my $h = shift @_;   #submission values hash
    my %h = %$h;    
    my @ar1;    #reference values array
    my @ar2;    #submission values array
    my $xy = 0;
    my $t;
    our $verbose;
    our $debug;
    if ($debug) {say "REF[i]\tSUB[j]";}
    foreach my $k ( sort keys %{$h{"Predictions"}} ) {
        push @ar1 , $r{$k};
        push @ar2 , $h{"Predictions"}{$k};
        if ($debug) {say $r{$k}."\t".$h{"Predictions"}{$k}};         
    }
    if ($debug) {say "\$xs = ( REF[i]-REF[j] ) / ( sqrt( (REF[i]-REF[j])^2+10^(-15) ) )"};
    if ($debug) {say "\$ys = ( SUB[i]-SUB[j] ) / ( sqrt( (SUB[i]-SUB[j])^2+10^(-15) ) )"};
    if ($debug) {say "\$xy += ( \$xs * \$ys )"};
    if ($debug) {say "REF[i]\tREF[j]\tSUB[i]\tSUB[j]\txs\tys\txy"};
    for (my $i = 0; $i < $#ar2; $i++){
        for (my $j = $i + 1; $j <= $#ar2; $j++) {
            my $xs = ($ar1[$i]-$ar1[$j])/(sqrt(($ar1[$i]-$ar1[$j])**2+10**(-15)));
            if ($debug) {say "DBG: REF[i]= $ar1[$i]\tREF[j]= $ar1[$j]\t\$xs = ($ar1[$i]-$ar1[$j])/(sqrt(($ar1[$i]-$ar1[$j])**2+10**(-15))) = $xs"};
            my $ys = ($ar2[$i]-$ar2[$j])/(sqrt(($ar2[$i]-$ar2[$j])**2+10**(-15)));
            if ($debug) {say "DBG: SUB[i]= $ar2[$i]\tSUB[j]= $ar2[$j]\t\$ys = ($ar2[$i]-$ar2[$j])/(sqrt(($ar2[$i]-$ar2[$j])**2+10**(-15))) = $ys"};
            if ($debug) {print "$xy + "};
            $xy += ($xs*$ys);
            if ($debug) {say "($xs * $ys) = $xy"};
            if ($debug) {say "$ar1[$i]\t$ar1[$j]\t$ar2[$i]\t$ar2[$j]\t$xs\t$ys\t$xy"};
        }        
    }
    my $n = scalar(@ar2);
    if ($debug) {say "\$n = ".$n};
    $t = 2*$xy/($n*($n-1));
    if ($debug) {say "TAU: 2 * $xy / ( $n * ($n - 1))"};
    if ($verbose || $debug) {print"\tTAU:\t$t\n"};
    if ($verbose || $debug) {say "\tTau evaluation is done\n"};
    return $t;
}

#Pearson Correlation
sub pea {       
    my $r = shift @_;   #reference values hash
    my %r = %$r;
    my $h = shift @_;   #submission values hash
    my %h = %$h;     
    my $s1 = 0;
    my $s2 = 0;
    my $num = 0;
    my $den1 = 0;
    my $den2 = 0;
    my $rxy;
    our $verbose;
    our $debug;
    my $n = scalar( keys %{$h{"Predictions"}});    
    
    if ($debug){print "HG\tREF\tSUB\n";}
    foreach my $k (sort keys %{$h{"Predictions"}}) {
        $s1 += %{$h{"Predictions"}}{$k};
        $s2 += ($r{$k});
        if ($debug){say $k."\t".$r{$k}."\t".%{$h{"Predictions"}}{$k}};        
    }
    if ($debug) {
        say "\$n = $n";
        say "\$s1 += SUB";
        say "\$s2 += REF";
        say "\$num += (REF - (\$s2 / \$n)) * (SUB - (\$s1 / \$n))";
        say "\$den2 += (SUB - (\$s1 / \$n))^2";
        say "\$den1 += (REF - (\$s2 / \$n))^2";
        say "\$rxy = \$num / ((sqrt(\$den1)) * (sqrt(\$den2)))";
    }
    if ($debug) {say "\$num\t\$den2\t\$den1"} 
    foreach my $k (sort keys %{$h{"Predictions"}}) {
        $num += ($r{$k} - ($s2 / $n)) * (%{$h{"Predictions"}}{$k} - ($s1 / $n));
        $den2 += (%{$h{"Predictions"}}{$k} - ($s1 / $n))**2;
        $den1 += ($r{$k} - ($s2 / $n))**2;
        if ($debug){print $num."\t".$den2."\t".$den1."\n"};
    }
    $rxy = $num / ((sqrt($den1)) * (sqrt($den2)));
    if ($debug) {say "\$rxy = $num / ((sqrt($den1)) * (sqrt($den2))) = $rxy"}
    if ($verbose || $debug) {say "\tPEARSON:\t$rxy\n\tPearson evaluation is done\n"};
    return $rxy;
}


sub intercept {
    my $r = shift @_;   #reference values hash
    my %r = %$r;
    my $h = shift @_;   #submission values hash
    my %h = %$h;  
    my $rxy = shift @_;
    our $verbose;
    our $debug; 
    my $n = scalar(%{$h{"Predictions"}});
    my $xm = 0;
    my $ym = 0;
    my $xm2 = 0;
    my $ym2 = 0;

    foreach my $k (sort keys %{$h{"Predictions"}}) {
        $xm += $r{$k};
        $ym += %{$h{"Predictions"}}{$k};
        $xm2 += ($r{$k})**2;
        $ym2 += (%{$h{"Predictions"}}{$k})**2;
    }
    $xm = $xm/$n;
    $xm2 = $xm2/$n; 
    $ym = $ym/$n;
    $ym2 = $ym2/$n;

    my $a = $rxy * sqrt($ym2 - $ym**2) / sqrt($xm2 - $xm**2);
    my $b = $ym - $a * $xm;
    if ($verbose || $debug) {print"\tAngular coefficient:\t$a\n\tIntercept:\t$b\n"};
    return ($a,$b);
}

sub lin {

    my $r = shift @_;   #reference values hash
    my %r = %$r;
    my $h = shift @_;   #submission values hash
    my %h = %$h;
    my $rxy = shift @_;
    my $s1 = 0; #sum of predictions
    my $s2 = 0; #sum of experimentals
    my $ss1 = 0;
    my $ss2 = 0;
    our $debug;
    our $verbose;
    my $n = scalar(%{$h{"Predictions"}});


    foreach my $k (sort keys %{$h{"Predictions"}}) {
            $s1 += %{$h{"Predictions"}}{$k};
            $s2 += ($r{$k});
            $ss1 += (%{$h{"Predictions"}}{$k})**2;
            $ss2 += (($r{$k}))**2;
            if ($debug){print "ref: ".$r{$k}."\n"}; #ref value 
            if ($debug){print "sub: ".%{$h{"Predictions"}}{$k}."\n"};        
            if ($debug){print "ss1: ".$ss1."\n"};
            if ($debug){print "ss2: ".$ss2."\n"};
    }

    my $sx = sqrt($ss1/$n - ($s1/$n)**2);
    if ($debug){print "Varianza x: ".$sx."\n"};
    my $sy = sqrt($ss2/$n - ($s2/$n)**2);
    if ($debug){print "Varianza y: ".$sy."\n"};
    if ($debug){print "Media x: ".$s1/$n."\n"};
    if ($debug){print "Media y: ".$s2/$n."\n"};
    if ($debug){print "Pears: ".$rxy."\n"};

    my $rc = (2 * $rxy * $sx * $sy) / ( $sx**2 + $sy**2 + ($s1/$n - $s2/$n)**2 );
    if ($verbose || $debug) {print"\tLIN:\t$rc\n"};
    return $rc;

}

sub prettyprint {

    my $h = shift(@_);
    my %h = %$h;
    my $r = shift(@_);
    my %r = %$r;
    our $file;
    our $evaluate;
    our $verbose;
    our $output;
    if ($evaluate) {
        my $mue = shift(@_);
        my %mue = %$mue;
        my $tau= shift(@_);   
        my %tau = %$tau;
        my $pea= shift(@_);   
        my %pea = %$pea;
        my $angular = shift(@_);
        my %angular = %$angular;
        my $intercept = shift(@_);
        my %intercept = %$intercept;
        my $lin = shift(@_);
        my %lin = %$lin;
    }

    
    open (OUT, '>>', $output);
    print OUT ">File\n$file\n";
    foreach my $k ( sort keys %h ) {
        print OUT ">".$k."\n";
        if ($k =~ /Prediction/) {
            print OUT "HG\tSUB\tREF\n";
            foreach my $k2 ( sort keys %{ $h{$k} } ) {
                print OUT $k2."\t".%{$h{$k}}{$k2}."\t".$r{$k2}."\n";        
            }
            if ($evaluate) {
                print OUT ">mue\n".$mue{$file}."\n";
                print OUT ">kendall\n".$tau{$file}."\n";
                print OUT ">pearson\n".$pea{$file}."\n";
                print OUT ">angular\n".$angular{$file}."\n";
                print OUT ">intercept\n".$intercept{$file}."\n";
                print OUT ">Lin\n".$lin{$file}."\n";
            }
        } else {
          if ($k eq "Software") {
            print OUT $h{$k}."\n";
          } else {
            foreach my $e (@{ $h{$k} }) {
                print OUT $e."\n";
            }  
          }
        }
    }
    print OUT "#\n";
}


#Check if @file contains selections (e.g Category:Alchemical)
sub check {
    my $arr = shift @_;
    my @arr = @$arr;
    my @ok;
    my %hash;
    my %checks;
    my $cons;
    our $debug;
    our $verbose;
    our $file;
    if ($debug || $verbose) {say "\tDoing Checks for:$file";}
    foreach my $s (@selects) {
        if ($debug || $verbose) {say "$s";}
        my ($a,$b) = split /:/,$s;
        $hash{$a} = $b;
        $checks{$a} = 0;
        foreach my $l (@arr) {
            if ($l =~ /$a/) {
                $checks{$a} = 1;
                next;
            }
            if ( $checks{$a} == 1 && $l eq "#") {
                last;
            }
            if ( $l !~ /^#$/ && $checks{$a} == 1 && $l =~ $b) {
                push @ok, $b;                
                last
            }
        }
    }
    if (scalar(@ok) != scalar(@selects)) {
        my $cons = 0;
        if ($verbose || $debug) {say "\tNegative check for $file";}
        return $cons;
    } else {
        my $cons = 1;
        if ($verbose || $debug) {say "\tPositive check for $file";}
        return $cons;
    }    
}

