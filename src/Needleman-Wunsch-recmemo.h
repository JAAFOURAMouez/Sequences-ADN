/**
 * \file Needleman-Wunsch-recmemo.h
 * \brief recursive implementation with memoization of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences 
 * \version 0.1
 * \date 03/10/2022 
 * \author Jean-Louis Roch (Ensimag, Grenoble-INP - University Grenoble-Alpes) jean-louis.roch@grenoble-inp.fr
 */

#include <stdlib.h> /* for size_t */

/**
 * \brief Calculates the minimum of three values.
 *
 * This macro evaluates and returns the smallest of three values.
 *
 * \param a The first value.
 * \param b The second value.
 * \param c The third value.
 */
#define MIN(a, b, c) ((a < b && a < c) ? a : ((b < c) ? b : c))

/**
 * \brief Returns the maximum of two values.
 * \param a The first value.
 * \param b The second value.
 * \return The larger of the two values.
 */
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/*
 * Costs for operations on canonical bases
 * Three  operations: insertion and sustitution of one base by an another 
 * Note= substitution of an unknown base N by another one (known or unknown) as the same cost than substitution between 2 different known bases
 */
/** \def SUBSTITUTION_COST
 *  \brief Cost of substitution of one canonical base by another
 */
#define SUBSTITUTION_COST	1

/** \def SUBSTITUTION_UNKNOWN_COST
 *  \brief Cost of substitution of an unknown base (N) by another one (canonical or unknown)
 */
#define SUBSTITUTION_UNKNOWN_COST	1  /* Cost for sustitition of an Unknown bas N by another on -known or unkown- */ 

/** \def INSERTION_COST
 *  \brief Cost of insertion of a canonical base 
 */
#define INSERTION_COST		2

/** 
 * \def K
 * \brief Size of the block for processing
 */
#define K 64

#define S 128

/********************************************************************************
 * Recursive implementation of NeedlemanWunsch with memoization
 */
/**
 * \fn long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB);
 * \brief computes the edit distance between A[0 .. lengthA-1] and B[0 .. lengthB-1]
 * \param A  : array of char represneting a genetic sequence A 
 * \param lengthA :  number of elements in A 
 * \param B  : array of char represneting a genetic sequence B
 * \param lengthB :  number of elements in B 
 * \return :  edit distance between A and B }
 *
 * editDistance_RecMemo is a memoized recursive immplementatioin of Needleman-Wunsch algorithm.
 * It allocates the data structure for memoization table and calls the internal recursive function _editDistance_memo
 * that fills in the memoization table.
 * 
 * If lengthA < lengthB, the sequences A and B are swapped.
 *
 */
long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB);

long EditDistance_NW_Iter(char* A, size_t lengthA, char* B, size_t lengthB);

long EditDistance_NW_Aware(char* A, size_t lengthA, char* B, size_t lengthB);

long EditDistance_NW_Oblivious(char* A, size_t lengthA, char* B, size_t lengthB);
