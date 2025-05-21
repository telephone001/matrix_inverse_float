

#include "matrix_inv.h"

/// @brief finds P_map and PA and initializes U for the next step.
/// @param n size of the matrices
/// @param U the U matrix we are trying to find
/// @param P_map the P map we update
static void matrix_find_p_map_and_swap(int n, float *A, float *PA, int *P_map)
{
        //find the biggest element in the column corresponding to the row we are finding
        for (int col = 0; col < n; col++) {
                float max_val = 0;
                int max_in_col = 0; //row with the max element
            
                for (int row = 0; row < n; row++) {
                        //Usually we have to copy the A matrix into U and swap the rows
                        // but this approach stops that.
                        int already_taken = 0;
                        for (int i = 0; i < col; i++) {
                                if (P_map[i] == row && col != 0) {
                                    already_taken = 1;
                                }
                        }
                        if (already_taken == 1) {
                                continue;
                        }
                        
                        //find the maximum valid pivot point for the current column
                        if (abs(A[n * row + col]) > abs(max_val)) {
                                max_in_col = row;
                                max_val = A[n * row + col];
                        }
                } 
                
                P_map[col] = max_in_col; 
        }
        
        
        //fill PA
        for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                        PA[n * row + col] = A[n * P_map[row] + col];
                }
        }
}

/// @brief fills out L and U matrix given an A matrix where PA = LU
///            A, L and U must be allocated before using this. 
///             Also pivoting must be done before using this.
/// @param n size of the matrix
/// @param A Input matrix pointer (doesn't allocate)
/// @param P_map output Pivot indices. shows what rows get swapped. index of array swaps with value (doesn't allocate)
/// @param l output L pointer (doesn't allocate)
/// @param U output U pointer (doesn't allocate)
static void matrix_lu_decompose(int n, float *A, int *P_map, float *L, float *U)
{
        //copy PA into U
        matrix_find_p_map_and_swap(n, A, U, P_map);
        
        
        //decompose A into inv(P)*LU
        for (int row = 0; row < n; row++) {
                
                //compute a row of U and L
                for (int col = 0; col < row; col++) {
                        //we try to make U[n * row + col] go to 0.
                        
                        float l_element; //the element that L[n*row+col] is
                        l_element = U[n * row + col] / U[n * col + col];

                        //elementary row operation to make it 0
                        for (int i = 0; i < n; i++) {
                                U[n * row + i] -= l_element * U[n * col + i]; 
                        }

                        //assign the L element
                        L[n * row + col] = l_element;
                }
                //assign the rest of L
                L[n * row + row] = 1; //diagonal is 1
                for (int col = row + 1; col < n - 1; col++) {
                        L[n * row + col] = 0;
                }
        }
}

/// @brief inverses a lower triangle matrix L. Modifies the array you give it.
/// @param n size of matrix (matrix is n*n)
/// @param L pointer to matrix
/// @param invL pointer to the output matrix
static void matrix_inverse_L(int n, float *L, float *invL)
{
        //forward substitution
        //find out what multiplied with L gives an identity matrix
        //[1 0 0]   [1 0 0]   [1 0 0]
        //[a 1 0] * [y 1 0] = [0 1 0]
        //[b c 1]   [z u 1]   [0 0 1]

        //each solved column is a system of equations
        //ex. solve xyz
        //x = 1
        //ax + y = 0
        //bx + cy + z = 0
        //
        //x = 1
        //y = -ax
        //z = -bx -cy

        //
        //
        //


        for (int x = 0; x < n; x++) {
                invL[x*n + x] = 1; //make diagonal correct term

                //if you think of x as our y coordinate, we can fill the zeros for each iteration  
                for (int i = x + 1; i < n; i++) {
                        invL[x*n + i] = 0;
                }


                for (int y = x+1; y < n; y++) {
                        //set to 0 cuz we subtract from this
                        invL[y*n + x] = 0;

                        for (int i = 0; i < y; i++) {
                                invL[y*n + x] -= L[y*n + i] * invL[i*n + x];
                        }
                }
        }
}


/// @brief inverses a lower triangle matrix U. Modifies the array you give it.
/// @param n size of matrix (matrix is n*n)
/// @param U pointer to matrix
/// @param invU output
static void matrix_inverse_U(int n, float *U, float *invU)
{
        //back substitution
        //find out what multiplied with U gives an identity matrix
        //[a b c]   [x w q]   [1 0 0]
        //[0 e f] * [y v r] = [0 1 0]
        //[0 0 i]   [z u s]   [0 0 1]

        //each solved column is a system of equations
        //ex. solve qrs
        //s*i = 1
        //r*e + s*f = 0
        //q*a + r*b + s*c = 0
        //
        //s= 1/i
        //r= -s*f / e
        //q= (-s*c - r*b) / a 

        for (int x = 0; x < n; x++) {
                invU[x*n + x] = 1 / U[x*n + x]; //make diagonal correct term

                for (int y = x - 1; y >= 0; y--) {
                        //set to 0 cuz we subtract from this
                        invU[y*n + x] = 0;

                        for (int i = n-1; i > y; i--) {
                                invU[y*n + x] -= U[y*n + i] * invU[i*n + x];
                        }
                        invU[y*n + x] /= U[y*n + y];
                }

                //if you think of x as our y coordinate, we can fill the zeros for each iteration

                for (int i = 0; i < x; i++) {
                        invU[x*n + i] = 0;
                }
        }
}

/// @brief gets inverse of A
///        decompose to PA = LU
///        inv(A) * inv(P) = inv(U) * inv(L)
///        inv(A) = inv(U) * inv(L) * P
/// @param n the side length of matrix
/// @param A is an n*n matrix represented as an array
/// @param invA inverted A output
void matrix_inverse(int n, float *A, float *invA)
{
        float *L = malloc(sizeof(float) * n * n);
        //U also acts like the inverse matrix
        float *U = malloc(sizeof(float) * n * n);
        float *invU = malloc(sizeof(float) * n * n);
        float *invL = malloc(sizeof(float) * n * n);
        int *P_map = malloc(sizeof(int) * n);

        //get L and U and P from A
        matrix_lu_decompose(n, A, P_map, L, U);

        int *inv_P_map = malloc(sizeof(int) * n);
        for (int i = 0; i < n; i++) {
                inv_P_map[P_map[i]] = i;
        }

        //inverse U and L
        matrix_inverse_U(n, U, invU);

        matrix_inverse_L(n, L, invL);


        free(P_map);
        free(L);
        free(U);



        //uses the P_map to swap the columns
        for (int y = 0; y < n; y++) {
                for (int x = 0; x < n; x++) {
                        invA[n*y + x] = 0;
                        for (int i = 0; i < n; i++) {
                                invA[n*y + x] += invU[n*y + i] * invL[n * i + inv_P_map[x]];
                        }
                }
        }
        // get inv()

        free(inv_P_map);
        free(invU);
        free(invL);


}



