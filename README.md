# Distributed & Parallel Computation


Final Project as part of course taken Distributed & Parallel Computation

## Description

This is a parallel and distributed implementation of sequence alignment using MPI, OpenMP and CUDA.
The requirment - for given 2 sequences (seq1,seq2) of unkown lengths find the offset (from the start of seq1) and mutant (by changing chars in seq2),
that produces the best alignment score (max/min).
Using MPI in order to communicate between 2 process and sending/reciving data - one of the project requirment was to be able to run the
program on 2 computers (using virtual machines) and reciving the same result.

## Getting Started
* Using Windows:
  * my advice is to use VScode - in order to be able to run it best use WSL2, install Microsoft MPI and CUDA.
* Using UNIX:
  * sudo apt-get update -y
  * sudo apt-get install -y mpi
  * for cuda - https://docs.vmware.com/en/VMware-vSphere-Bitfusion/3.0/Example-Guide/GUID-ABB4A0B1-F26E-422E-85C5-BA9F2454363A.html
### Executing program

* There is a makefile - in the same dirctory as the executable open a terminal
  * make build - to compile and build the project
  * make run - executes the program on the same machine using 2 process
  * make runOn2 - allows to run on 2 VM (there's a file "mf" where you need to add the ip's of the 2 computer using hostname -I)


## Authors

Daniel Aizenband

## Version History

* 0.1
    * Initial Release
