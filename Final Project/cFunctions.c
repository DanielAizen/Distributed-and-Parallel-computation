#include "myProto.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int readFromFile(char const *fileName, double *weights, char *temp_seq1, char *temp_seq2)
{
	FILE *fp = fopen(fileName, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "Failed to open Input.txt\n");
		exit(EXIT_FAILURE);
	}

	fscanf(fp, "%lf %lf %lf %lf", &weights[0], &weights[1], &weights[2], &weights[3]);
	fscanf(fp, "%s", temp_seq1);
	fscanf(fp, "%s", temp_seq2);
	
	char temp[8];
	fscanf(fp, "%s", temp);
	
	if (strcmp(temp, "maximum") == 0)
	{
		return 1;
	}
	fclose(fp);
	return -1;
}

void CreateMutant(char *seq1_data, char *temp_Seq2, int indexSeq2, int sort_order, double *weights, int offset)
{
	
	char tmp;
	if (sort_order == 1)
	{ //max - check only if both chars aren't in the same conservative group, if not I
	  //switch to the char from seq1 in order to add to the score (ASTRIX)
		if (IsConservative(seq1_data[indexSeq2 + offset], temp_Seq2[indexSeq2]) == 0)
		{
			temp_Seq2[indexSeq2] = seq1_data[indexSeq2 + offset];
		}
	}
	else
	{ 
		if (IsConservative(seq1_data[indexSeq2 + offset], temp_Seq2[indexSeq2]) == 0){
			//   w3-semi		   w4-no match
			if (weights[2] > weights[3])
			{
				tmp = IsSemiConservative(seq1_data[indexSeq2 + offset]);
				if (tmp != '$' && IsConservative(temp_Seq2[indexSeq2],tmp) == 0)
					temp_Seq2[indexSeq2] = tmp;
			}
			else{
				temp_Seq2[indexSeq2] = 'Z';
			}
		}
			
	}
}

int IsConservative(char c1, char c2)
{
	const char *conservative_Group[CONSERVATIVE] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
	const char *cmp1, *cmp2;
	for (int i = 0; i < CONSERVATIVE; i++)
	{
		cmp1 = strchr(conservative_Group[i], c1);
		cmp2 = strchr(conservative_Group[i], c2);
		if (cmp1 != NULL && cmp2 != NULL)
			return 1;
	}
	return 0;
}

char IsSemiConservative(char a)
{
	const char *semi_Conservative_Group[SEMI_CONSERVATIVE] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};
	const char *cmp1;
	
	for (int i = 0; i < SEMI_CONSERVATIVE; i++)
	{
		cmp1 = strchr(semi_Conservative_Group[i], a);
		if (cmp1 != NULL){
			for(int j=0; j< strlen(semi_Conservative_Group[i]); j++){
				if(semi_Conservative_Group[i][j] != a && IsConservative(a , semi_Conservative_Group[i][j]) == 0)
					return semi_Conservative_Group[i][j];	
			}
		}
	}
	return '$';
}

double CalculateScore(double *weights, char *seq1, char *temp_Seq2, int offset)
{
	char c1, c2;
	int numOfWeights[4] = {0,0,0,0};
	for (int idx = 0; idx < strlen(temp_Seq2); idx++)
	{
		c1 = seq1[idx + offset];
		c2 = temp_Seq2[idx];
		int index = CheckSign(c1, c2);
		numOfWeights[index]++;
	}
	float res= (numOfWeights[0] * weights[0] - numOfWeights[1] * weights[1] - numOfWeights[2] * weights[2] - numOfWeights[3] * weights[3]);
	return res;
}

int CheckSign(char c1, char c2)
{
	const char *conservative_Group[CONSERVATIVE] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
	const char *semi_Conservative_Group[SEMI_CONSERVATIVE] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};
	const char *cmp1, *cmp2;

	if (c1 == c2) 
		return ASTRIX;

	for (int i = 0; i < CONSERVATIVE; i++)
	{
		cmp1 = strchr(conservative_Group[i], c1);
		cmp2 = strchr(conservative_Group[i], c2);
		if (cmp1 != NULL && cmp2 != NULL)
			return COLON;
	}
	for (int i = 0; i < SEMI_CONSERVATIVE; i++)
	{
		cmp1 = strchr(semi_Conservative_Group[i], c1);
		cmp2 = strchr(semi_Conservative_Group[i], c2);
		if (cmp1 != NULL && cmp2 != NULL)
			return DOT;
	}
	return SPACE;
}

void compareResults(int sort_order, char *score_mutantSeq, char *mutantSeq2, double *best_score, double *result, int *curr_offset, int *score_offset){
	if (sort_order == 1)
		{
			if (*result > *best_score)
			{
				*best_score = *result;
				*score_offset = *curr_offset;
				strcpy(score_mutantSeq, mutantSeq2);
			}
		}
		else
		{
			if (*result < *best_score)
			{
				*best_score = *result;
				*score_offset = *curr_offset;
				strcpy(score_mutantSeq, mutantSeq2);
			}
		}
}

void writeToFile(char const *fileName, char *mutantSeq, int scoreOffset, double score){
	FILE *fp = fopen(fileName, "w");
	if (fp == NULL)
	{
		fprintf(stderr, "Failed to open Output.txt\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fp,"Mutant sequence: %s\n",mutantSeq);
	fprintf(fp,"offset is: %d\tallinment score: %lf",scoreOffset,score);
	fclose(fp);
}