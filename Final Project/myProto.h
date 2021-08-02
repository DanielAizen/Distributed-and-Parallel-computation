#pragma once


#define MAX_SEQ1 10000
#define MAX_SEQ2 5000
#define ASTRIX 0
#define COLON 1
#define DOT 2
#define SPACE 3
#define CONSERVATIVE 9
#define SEMI_CONSERVATIVE 11



/*
 * readFromFile:
 * This fucntion recives the desired file name to be opened, and an array of weights to be read.
 * Additionally it recives 2 char pointer (seq1, seq2) that will hold the
 * Seq1 and Seq2 that are given in the input file.
 * Based on the fourth line the programs will know what score to give back,
 * maximin or minimum. 
 */
int readFromFile(char const *fileName, double* weights, char *temp_seq1, char *temp_seq2);

/*
 * createMutant:
 * This function recives the sequnces and specific offset from the beggining of seq1,
 * a copy array of seq2 to be mutated based on the sort order that is given (max/min). 
 * In the temp_seq2 I specify the index of which char I want to mutate.
 * This is where I parallize the program using OpnenMP, by sending each char of temp_seq2. 
 * Complexity evaluation: O(C)
 */
void CreateMutant(char *seq1_data, char *mutant_seq, int indexSeq2, int sort_order, double *weights, int offset);

/*
 * IsConservative:
 * This function recives chars from both sequnces, using strchr I find if the char from each seq
 * exists in the same conservative group.
 * returns 1 if both chars are in the same group, 0 if not
 * Complexity evaluation: O(C)
 */
int IsConservative(char c1, char c2);

/*
 * IsSemiConservative:
 * This function recives char from seq1 at the index of the offset + the index of seq2.
 * By using strchr I find if this char is in any semi-conservative group, if it is, I check
 * the entire group for a char that is not the same char as i recived and also not in the same
 * conservative group.
 * returns the char the can be substituted if succesful, and $ if not.
 *  Complexity evaluation: O(C)
 */
char IsSemiConservative(char a);

/*
 * CalculateScore:
 * This function recives chars from both sequnces, the weights and the offset from seq1.
 * for each char from seq2 that is under seq1 (meaning running of the length of seq2), I send
 * both char to CheckSign that returns the type of match are those chars.
 * returns the calculated allingment score of each mutant at a specific offset,
 * based on an array of ints that represents the signs. 
 * Complexity evaluation: O(n)
 */
double CalculateScore(double *weights, char *seq1, char *temp_Seq2, int offset);

/*
 * CheckSign:
 * Recives 2 chars and checks what type of score this 2 chars gives.
 * Returns a defined value which later is used in CalculateScore.
 * Complexity evaluation: O(C)
 */
int CheckSign(char c1, char c2);

/*
 * compareResults:
 * The purpose of this function is to compare the results in 2 instances.
 * 1)for each process between it's best results (max/min)
 * 2)between the slave and the master after both finished their workload (each works on different offsets) then choose the highest score
 *  where in case both return the same score the master will be chosen by default.  
 */
void compareResults(int sort_order, char *score_mutantSeq, char *mutantSeq2, double *best_score, double *result, int *curr_offset, int *score_offset);

/*
 * writeToFile:
 * This fucntion recives the desired *mutated* sequence after comparing with the slave process.
 * First line is used to write to the file the best mutant seq.
 * Second line is for the offset of that occurence and the highest/lowest score.
 */
void writeToFile(char const *fileName, char *mutantSeq, int scoreOffset, double score);

/*
 * computeOnGPU:
 * This function recives half of the workload and allocates memory on the GPU.
 * After coping the arrays the computeOnGP invokes CreateMutationGPU which works
 * identically to CreateMutant.
 * Becuase the gpu doesn't know the string.h library I added impllematation to 2
 * function (strlen & strchr), additionally i added IsConservativeGPU and IsSemiConservativeGPU
 * to mimic the same idea I used on the gpu.
 */
int computeOnGPU(char *h_Seq1, char *h_Seq2, int h_Offset, int sort_Order, double *h_Weights, int h_Len_Seq1, int h_Len_Seq2);

//GPU functions:
/*
 * __device__ int IsConservativeGPU(char c1, char c2);

 * __device__ char IsSemiConservativeGPU(char s);

 * __device__ int MyStrlen(const char *str);

 * __device__ int MyStrchr(const char *s, char c);
 */

