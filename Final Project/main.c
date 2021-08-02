#include <mpi.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "myProto.h"

int main(int argc, char *argv[])
{
	int size = 0, rank, i, numOfProc;
	MPI_Status status;

	/*
	 * Setting the max size allowed for both sequences
	 * Initializing the parameters that are going to read from the file,
	 * allocating space inside the master
	 */

	int sort_Order = 1;
	double weights[4];
	char *seq1_data, *seq2_data;
	int seq1_count, seq2_count;
	
	int score_offset =0;
	double score;
	char *score_mutantSeq;
	char *mutantSeq2;

	int totOffset = 0, offset =0;
	int startIdx = 0;
	int endIdx = 0;
	double result=0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
	if (numOfProc != 2)
	{
		printf("Run the example with two processes only\n");
		MPI_Abort(MPI_COMM_WORLD, __LINE__);
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		
	if (rank == 0)
	{
		//first allocate MAX_SEQ1 & MAX_SEQ2 to the seqeunces the realloc to the actual size
		seq1_data = (char *)malloc(MAX_SEQ1*sizeof(char));
		if(seq1_data == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		seq2_data = (char *)malloc(MAX_SEQ2*sizeof(char));
		if(seq2_data == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		
		sort_Order = readFromFile("Input1.txt", weights, seq1_data, seq2_data);
		
		seq1_count = strlen(seq1_data);
		if ((realloc(seq1_data,seq1_count * sizeof(char))) == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		
		seq2_count = strlen(seq2_data);
		if ((realloc(seq2_data,(seq2_count) * sizeof(char))) == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		
		totOffset = seq1_count - seq2_count + 1; //calculate the total number of iterations based on the offset
		offset = totOffset / 2;
		endIdx = offset; // index for master when to end
		
		score_mutantSeq = (char *)malloc((seq2_count) * sizeof(char));
		if (score_mutantSeq == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
	}

	MPI_Bcast(&sort_Order, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(&offset, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(&totOffset, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(&seq1_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(&seq2_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(weights, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(rank != 0){
		startIdx = offset;
		endIdx = totOffset;

		if ((seq1_data = (char *)malloc(seq1_count * sizeof(char))) == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		if ((seq2_data = (char *)malloc(seq2_count * sizeof(char))) == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		
		score_mutantSeq = (char *)malloc(seq2_count * sizeof(char));
		if (score_mutantSeq == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);

	}
	MPI_Bcast(seq1_data, seq1_count, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(seq2_data, seq2_count, MPI_CHAR, 0, MPI_COMM_WORLD);
	

	mutantSeq2 = (char *)malloc(seq2_count * sizeof(char));
	if (mutantSeq2 == NULL)
		MPI_Abort(MPI_COMM_WORLD, __LINE__);
	score = CalculateScore(weights, seq1_data,  seq2_data, 0); // calculate score for the first time with original seq2
	for (int idx = startIdx; idx < endIdx; idx++)
	{
		strcpy(mutantSeq2, seq2_data);
	
#pragma omp parallel for
		for (int i = 0; i < seq2_count/2; i++)
		{ // send to omp half of workload
			CreateMutant(seq1_data, mutantSeq2, i, sort_Order, weights, idx);
		}
		
		//sending to CUDA the 2nd half for mutation
		if(seq2_count % 2 == 0 ){
			//				(----2nd half of seq1-----, ----2nd half of seq2-------,offset,max/min, -weights[4]-, --size of half of seq1--,  -size of half of seq2-)
			if(computeOnGPU(seq1_data + (seq2_count/2), mutantSeq2 + (seq2_count/2), idx, sort_Order, weights,   seq1_count - (seq2_count/2), (seq2_count/2)) != 0)
				MPI_Abort(MPI_COMM_WORLD, __LINE__);
		}
		else{
			//			   (----2nd half of seq1------, ----2nd half of seq2--------,offset, max/min, -weights[4]-, ----size of half of seq1 +1 ----, size of half of seq1 +1)
			if(computeOnGPU(seq1_data + (seq2_count/2), mutantSeq2 + (seq2_count/2), idx, sort_Order, weights, (seq1_count - (seq2_count/2) +1 ), ((seq2_count/2) +1))!= 0)
				MPI_Abort(MPI_COMM_WORLD, __LINE__);
		}
		result = CalculateScore(weights, seq1_data, mutantSeq2, idx);
		compareResults(sort_Order,score_mutantSeq,mutantSeq2,&score,&result, &idx, &score_offset);
				
	}
	if (rank == 0)
	{
		int slave_scoreOffset;
		double slave_score;
		char *slave_mutantSeq;
		slave_mutantSeq = (char *)malloc(seq2_count * sizeof(char));
		if (slave_mutantSeq == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);

		MPI_Recv(slave_mutantSeq, seq2_count, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&slave_scoreOffset, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&slave_score, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
		compareResults(sort_Order, score_mutantSeq, slave_mutantSeq, &score,&slave_score, &slave_scoreOffset, &score_offset);
		
		// writeToFile("Output.txt",score_mutantSeq,score_offset,score);
		writeToFile("Output1.txt",score_mutantSeq,score_offset,score);
		// writeToFile("Output2.txt",score_mutantSeq,score_offset,score);
		free(slave_mutantSeq);
	}
	else
	{
		MPI_Send(score_mutantSeq, seq2_count, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&score_offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&score, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	
	free(seq1_data);
	free(seq2_data);
	free(score_mutantSeq);
	free(mutantSeq2);

	MPI_Finalize();

	return 0;
}

