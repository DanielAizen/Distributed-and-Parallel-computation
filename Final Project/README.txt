@author: Daniel Aizenband
@id: 205455520
@due date: 1/8/21

DOCUMANTAION:
-------------

This is a parallel and distributed impllamentation of sequence alignment using MPI, OpenMP and CUDA.

The requirment - for given 2 sequences (seq1,seq2) find the offset and mutant that produces the best alignment score (max/min).

By using MPI I split the workload between each process after reading the data in the main process ("master") and sending it to
the secondary process ("slave") by using MPI_Bcast (I used Bcast because each process needs to recive the entire seq1,seq2 and other
parameters in order to continue the work).

In order to parallize the program I used OpenMP, in the main, I first loop over the offset of each process (ex. seq1 length=46 seq2 length= 31,
the total offset will be 15 (seq1_count - seq2_count +1)-> which then I divide by 2 -> offset =7 (rounded down). Which means that the master will
be in charge offsets 0-7 (7 not included) and the slave will be in charge of offsets 7-15). Then by using * #pragma omp parallel for * I 
iterate over half of seq2 and send it to be mutated - each theard is in charge of mutating only one char from seq2 based on the constraints given in the
assignment.
Given sort order "maximum" - first I check that the chars at the current offset and index of seq2 aren't in the same conservative group, if not change the
							 char from seq2 to be the same char as seq1 (who's at the same offset + index).
Given sort order "minimun" - first  I check that the chars at the current offset and index of seq2 aren't in the same conservative group, if not I check
							 the weights (between semi-conservative and no match) then what char is allowed, and check if it's not in the same conservative group.

By using CUDA i send the other half of the workload (seq1 and seq2) to the GPU, that follows the same logic as I did in the CPU.

The Complexity evaluation is as such: O(n*m) - where n is the offset between seq1 and seq2, and m is the length of seq2.