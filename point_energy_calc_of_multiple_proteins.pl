#!/bin/perl                                                                                                                                                                                                                                   
use warnings;

# These sequences of commands simulate multiple (three) proteins in a vacuum using Gromacs.                                                                                                                                                   
# Make sure that the necessary .mdp file is ported for the run to work.                                                                                                                                                                       

# Perl script: loop over files in a folder                                                                                                                                                                                                    

print "Watch and Learn";

my @files = <*.pdb>;

foreach my $file (@files) {

    my ($base,$ext) = split (/\./,$file);
    print "filename is $file\n";
    print "basename = $base \n";
    print "extension = $ext \n";

system("gmx pdb2gmx -ff charmm36-nov2016 -f $file -o $base.gro -v -ignh -water none");

system("gmx editconf -f $base.gro -o $base.gro -c -d 10 -bt dodecahedron");

system("gmx grompp -f modified_mdp_file_1_for_vacuum.mdp -p topol.top -c $base.gro -o md_vacuum.tpr");

system("gmx mdrun -s md_vacuum.tpr -rerun $base.gro");

system("mv md.log $base.log");

}

