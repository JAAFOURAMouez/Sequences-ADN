/**
 * \file Needleman-Wunsch-recmemo.c
 * \brief recursive implementation with memoization of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences 
 * \version 0.1
 * \date 03/10/2022 
 * \author Jean-Louis Roch (Ensimag, Grenoble-INP - University Grenoble-Alpes) jean-louis.roch@grenoble-inp.fr
 *
 * Documentation: see Needleman-Wunsch-recmemo.h
 * Costs of basic base opertaions (SUBSTITUTION_COST, SUBSTITUTION_UNKNOWN_COST, INSERTION_COST) are
 * defined in Needleman-Wunsch-recmemo.h
 */


#include "Needleman-Wunsch-recmemo.h"
#include <stdio.h>  
#include <stdlib.h> 
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */

/*****************************************************************************/
   
/* Context of the memoization : passed to all recursive calls */
/** \def NOT_YET_COMPUTED
 * \brief default value for memoization of minimal distance (defined as an impossible value for a distance, -1).
 */
#define NOT_YET_COMPUTED -1L 

/** \struct NW_MemoContext
 * \brief data for memoization of recursive Needleman-Wunsch algorithm 
*/
struct NW_MemoContext 
{
    char *X ; /*!< the longest genetic sequences */
    char *Y ; /*!< the shortest genetic sequences */
    size_t M; /*!< length of X */
    size_t N; /*!< length of Y,  N <= M */
    long **memo; /*!< memoization table to store memo[0..M][0..N] (including stopping conditions phi(M,j) and phi(i,N) */
} ;

/*
 *  static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
 * \brief  EditDistance_NW_RecMemo :  Private (static)  recursive function with memoization \
 * direct implementation of Needleman-Wursch extended to manage FASTA sequences (cf TP description)
 * \param c : data passed for recursive calls that includes the memoization array 
 * \param i : starting position of the left sequence :  c->X[ i .. c->M ] 
 * \param j : starting position of the right sequence :  c->Y[ j .. c->N ] 
 */ 
static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
/* compute and returns phi(i,j) using data in c -allocated and initialized by EditDistance_NW_Rec */
{
   if (c->memo[i][j] == NOT_YET_COMPUTED)
   {  
      long res ;
      char Xi = c->X[i] ;
      char Yj = c->Y[j] ;
      if (i == c->M) /* Reach end of X */
      {  if (j == c->N) res = 0;  /* Reach end of Y too */
         else res = (isBase(Yj) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i, j+1) ;
      }
      else if (j == c->N) /* Reach end of Y but not end of X */
      {  res = (isBase(Xi) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i+1, j) ;
      }
      else if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
      {  ManageBaseError( Xi ) ;
         res = EditDistance_NW_RecMemo(c, i+1, j) ;
      }
      else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
      {  ManageBaseError( Yj ) ;
         res = EditDistance_NW_RecMemo(c, i, j+1) ;
      }
      else  
      {  /* Note that stopping conditions (i==M) and (j==N) are already stored in c->memo (cf EditDistance_NW_Rec) */ 
         long min = /* initialization  with cas 1*/
                   ( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                          : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
                   )
                   + EditDistance_NW_RecMemo(c, i+1, j+1) ; 
         { long cas2 = INSERTION_COST + EditDistance_NW_RecMemo(c, i+1, j) ;      
           if (cas2 < min) min = cas2 ;
         }
         { long cas3 = INSERTION_COST + EditDistance_NW_RecMemo(c, i, j+1) ;      
           if (cas3 < min) min = cas3 ; 
         }
         res = min ;
      }
       c->memo[i][j] = res ;
   }
   return c->memo[i][j] ;
}

/* EditDistance_NW_Rec :  is the main function to call, cf .h for specification 
 * It allocates and initailizes data (NW_MemoContext) for memoization and call the 
 * recursivefunction EditDistance_NW_RecMemo 
 * See .h file for documentation
 */
long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB)
{
   _init_base_match() ;
   struct NW_MemoContext ctx;
   if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
   {  ctx.X = A ;
      ctx.M = lengthA ;
      ctx.Y = B ;
      ctx.N = lengthB ;
   }
   else
   {  ctx.X = B ;
      ctx.M = lengthB ;
      ctx.Y = A ;
      ctx.N = lengthA ;
   }
   size_t M = ctx.M ;
   size_t N = ctx.N ;
   {  /* Allocation and initialization of ctx.memo to NOT_YET_COMPUTED*/
      /* Note: memo is of size (N+1)*(M+1) but is stored as (M+1) distinct arrays each with (N+1) continuous elements 
       * It would have been possible to allocate only one big array memezone of (M+1)*(N+1) elements 
       * and then memo as an array of (M+1) pointers, the memo[i] being the address of memzone[i*(N+1)].
       */ 
      ctx.memo = (long **) malloc ( (M+1) * sizeof(long *)) ;
      if (ctx.memo == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo" ); exit(EXIT_FAILURE); }
      for (int i=0; i <= M; ++i) 
      {  ctx.memo[i] = (long*) malloc( (N+1) * sizeof(long));
         if (ctx.memo[i] == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo[i]" ); exit(EXIT_FAILURE); }
         for (int j=0; j<=N; ++j) ctx.memo[i][j] = NOT_YET_COMPUTED ;
      }
   }    
   
   /* Compute phi(0,0) = ctx.memo[0][0] by calling the recursive function EditDistance_NW_RecMemo */
   long res = EditDistance_NW_RecMemo( &ctx, 0, 0 ) ;
    
   { /* Deallocation of ctx.memo */
      for (int i=0; i <= M; ++i) free( ctx.memo[i] ) ;
      free( ctx.memo ) ;
   }
   return res ;
}

long EditDistance_NW_Iter(char *A, size_t lengthA, char *B, size_t lengthB)
{
   _init_base_match();
   size_t sizeM = lengthA;
   size_t sizeN = lengthB;

   /* Initialisation of arrays to store intermediate results */
   long *columnM = malloc((sizeM + 1) * sizeof(long));
   long *columnN = malloc((sizeN + 1) * sizeof(long));
   columnM[sizeM] = 0;
   columnN[sizeN] = 0;

   for (int i = sizeM - 1; i >= 0; i--)
   {
      columnM[i] = (isBase(A[i]) ? INSERTION_COST : 0) + columnM[i + 1];
   }

   for (int j = sizeN - 1; j >= 0; j--)
   {
      columnN[j] = (isBase(B[j]) ? INSERTION_COST : 0) + columnN[j + 1];
   }

   int col_index = sizeN - 1;
   long right_element = columnN[sizeN - 1];
   long diagonal_element = columnN[sizeN];
   long result = 0;
   long mismatch_cost = 0;

   while (col_index >= 0)
   {
      if (!isBase(B[col_index]))
      {
         columnN[col_index + 1] = columnM[0];
         col_index--;
      }
      else
      {
         right_element = columnN[col_index];
         diagonal_element = columnN[col_index + 1];
         int row_index = sizeM - 1;

         while (row_index >= 0)
         {
            if (!isBase(A[row_index])){
               columnM[row_index] = right_element;
               row_index--;
            }else{
               mismatch_cost = (isUnknownBase(A[row_index]) || !isSameBase(A[row_index], B[col_index])) ? 1 : 0;
               result = MIN(2 + columnM[row_index], 2 + right_element, diagonal_element + mismatch_cost);
               diagonal_element = columnM[row_index];
               right_element = result;
               columnM[row_index] = result;
               row_index--;
            }
         }
         columnN[col_index + 1] = diagonal_element;
         col_index--;
      }
   }

   long final_result = columnM[0];
   free(columnM);
   free(columnN);
   return final_result;
}

long EditDistance_NW_Aware(char *A, size_t lengthA, char *B, size_t lengthB){
   _init_base_match();
   size_t m = lengthA;
   size_t n = lengthB;
   long *columnM = malloc((m + 1) * sizeof(long));
   long *columnN = malloc((n + 1) * sizeof(long));
   columnM[m] = 0;
   columnN[n] = 0;

   for (int i = m - 1; i >= 0; i--)
   {
      columnM[i] = (isBase(A[i]) ? INSERTION_COST : 0) + columnM[i + 1];
   }

   for (int j = n - 1; j >= 0; j--)
   {
      columnN[j] = (isBase(B[j]) ? INSERTION_COST : 0) + columnN[j + 1];
   }

   long result = 0;
   long mismatch_cost = 0;
   long diagonal_element = 0;
   long right_element = 0;

   for (int outerJ = m - 1; outerJ >= 0; outerJ -= K)
   {
      for (int outerI = n - 1; outerI >= 0; outerI -= K)
      {
         right_element = columnN[outerI];
         diagonal_element = columnN[outerI + 1];

         int indexB = outerI;

         while (indexB >= 0 && indexB > outerI - K)
         {
            if (isBase(B[indexB]))
            {
               int indexA = outerJ;
               right_element = columnN[indexB];
               diagonal_element = columnN[indexB + 1];

               while (indexA >= 0 && indexA > outerJ - K)
               {
                  if (isBase(A[indexA]))
                  {
                     mismatch_cost = (isUnknownBase(A[indexA]) || !isSameBase(A[indexA], B[indexB])) ? 1 : 0;
                     result = MIN(2 + columnM[indexA], 2 + right_element, diagonal_element + mismatch_cost);
                     diagonal_element = columnM[indexA];
                     right_element = result;
                     columnM[indexA] = result;
                  }
                  indexA--;
               }
               columnN[indexB + 1] = diagonal_element;
            }
            else
            {
               columnN[indexB + 1] = columnM[outerJ < K ? 0 : outerJ - K + 1];
            }
            indexB--;
         }
      }
   }

   long final_result = columnM[0];
   free(columnM);
   free(columnN);
   return final_result;
}

void  Process_Block(int start_blockA, int end_blockA, int start_blockB, int end_blockB, long *columnN, long *columnM, char *A, char *B)
{
    _init_base_match();
    long right_element = 0;
    long diagonal_element = 0;
    int indexB = end_blockB;
    long mismatch_cost = 0;
    long result = 0;
    while (indexB >= start_blockB) {
        if (!isBase(B[indexB])) {
            columnN[indexB + 1] = columnM[start_blockA];
        } else {
            int indexA = end_blockA;
            right_element = columnN[indexB];
            diagonal_element = columnN[indexB + 1];
            while (indexA >= start_blockA) {
                if (!isBase(A[indexA])) {
                    columnM[indexA] = right_element;
                } else {
                    mismatch_cost = (isUnknownBase(A[indexA]) || !isSameBase(A[indexA], B[indexB])) ? 1 : 0;
                    result = MIN(2 + columnM[indexA], 2 + right_element, diagonal_element + mismatch_cost);
                    diagonal_element = columnM[indexA];
                    right_element = result;
                    columnM[indexA] = result;
                }
                --indexA;            
            }
            columnN[indexB + 1] = diagonal_element;
        }
        --indexB;
    }
}

void EditDistance_NW_Oblivious_Rec(int start_blockA, int end_blockA, int start_blockB, int end_blockB, long *columnN, long *columnM, char *A, char *B)
{
   int n =  end_blockB - start_blockB;
   int m =  end_blockA - start_blockA;
   if (n <= S && m <= S){
      Process_Block(start_blockA, end_blockA, start_blockB, end_blockB, columnN, columnM, A,B);
   }else if (n > S){
      int mid = end_blockB - n/2;
      EditDistance_NW_Oblivious_Rec(start_blockA, end_blockA, mid  , end_blockB, columnN, columnM, A, B);
      EditDistance_NW_Oblivious_Rec(start_blockA, end_blockA, start_blockB, mid - 1, columnN, columnM, A, B);
   }else if ( m > S){
      int mid = end_blockA - m/2;
      EditDistance_NW_Oblivious_Rec(mid , end_blockA, start_blockB, end_blockB, columnN, columnM, A, B);
      EditDistance_NW_Oblivious_Rec(start_blockA, mid - 1, start_blockB, end_blockB, columnN, columnM, A, B);
   }
}

long EditDistance_NW_Oblivious(char *A, size_t lengthA, char *B, size_t lengthB)
{
   _init_base_match();
   int n = (int) lengthB;
   int m = (int) lengthA;
   long *columnM = malloc((m + 1) * sizeof(long));
   long *columnN = malloc((n + 1) * sizeof(long));
   columnM[m] = 0;
   columnN[n] = 0;
   for (int i = m - 1; i >= 0; i--){
      columnM[i] = (isBase(A[i]) ? INSERTION_COST : 0) + columnM[i + 1];
   }
   for (int j = n - 1; j >= 0; j--){
      columnN[j] = (isBase(B[j]) ? INSERTION_COST : 0) + columnN[j + 1];
   }
   EditDistance_NW_Oblivious_Rec(0, m - 1, 0, n - 1, columnN, columnM, A, B);
   long result = columnM[0];
   free(columnM);
   free(columnN);
   return result;
}