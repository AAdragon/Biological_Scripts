#!/bin/perl
use strict;
use warnings;

# This is a perl script. 
# Created by Adron Ung, a senior Biochemistry Major, at South Dakota State University on December 3rd, 2016
# with help from Dr. Brian Moore at South Dakota State University, as well. 
# If you have questions, Contact author: adron.ung@jacks.sdstate.edu 
# Please cite this script in publications to encourage me to write more scripts


#This script works for Gromacs verson 5.1.3. 
# These sequences of commands simulate multiple (three) proteins in a vacuum using Gromacs.                                                                                                                                                   
# Make sure that the necessary .mdp file is ported for the run to work.                                                                                                                                                                       
# To execute this script with prior assumptions: Source GMXRC, then obtain the files that you would like to analyze via curl. 
# Type something like this directory into a linux terminal to retrieve several proteins. 
# This is a bash script so it's commented out in the perl script:

# for id in 1AKI 1UBQ 2N9M
# do
# curl "http://files.rcsb.org/view/$id.pdb" - o $id.pdb
# done

# If you do that you will get these files: 1UBQ.pdb, 1AKI.pdb, and 2N9M.pdb.
# It does not matter what the names of your PDB files are for this perl script to work. 
# Then remove the Heteroatom HOH lines of water in a text editor and you are able to execute this script. 
# Note: I may have to make some changes to the .mdp file in a future to model the physics more accurately.
# Nevertheless, this is a great demonstration of how to call perl within Gromacs. 

# The force-field that I specified is CHARMM36. This is not integrated into Gromacs, so I suggest 
# obtaining it from the authors via here: http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs
# Then you will have to untar it where all the other default force-fields are. 

# Perl script: loop over files in a folder                                                                                                                                                                                                    

print "Watch and Learn";

my @files = <*.pdb>;

# Splits up the pdb files in a folder to the base names and extensions, 
# so 1UBQ.pdb becomes basename: 1UQB with extensions .pdb, etc. 
# foreach loop starts here: 

foreach my $file (@files) {

    my ($base,$ext) = split (/\./,$file);
    print "filename is $file\n";
    print "basename = $base \n";
    print "extension = $ext \n";

# Iteratively, the .pdb file is fed into gromacs, the .pdb extensions is rewritten as .gro with gromacs commands.
# system calls are used to tell perl to write these commands in the terminal. 

system("gmx pdb2gmx -ff charmm36-nov2016 -f $file -o $base.gro -v -ignh -water none"); # think of how you might use other force-fields. 

system("gmx editconf -f $base.gro -o $base.gro -c -d 10 -bt dodecahedron");

system("gmx grompp -f modified_mdp_file_1_for_vacuum.mdp -p topol.top -c $base.gro -o md_vacuum.tpr");

system("gmx mdrun -s md_vacuum.tpr -rerun $base.gro");

system("mv md.log $base.log"); # by default each protein gets a log file as md.log, so as to not overwrite
# the md.log files, it is neccessary to rename them as 1UBQ.pdb, 2N9M.pdb, 1AKI.pdb, etc ... 

} #foreach loop ends here. 
# If you would like to perform solvations, add more systems calls within the foreach loop
# , and you will have to use a different .mdp file. 



