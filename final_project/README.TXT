#7 - Cholesky decomposition of a real symmetric positive-definite matrix

A real symmetric positive-definite matrix, A, has been factorized using Cholesky decomposition. 
The matrix L, for which A=LL^T, has been deterimined using a Chilesky-Crout algorithm. In this algorithm the elements of L are calutlated starting with the 
upper-left-most element and continued by calulating elements column by column. 

In my program the functions used are found in the file "funcs.c" and the decomposition is done in "main.c". The Matrix, A, used in the project is defined in the
main file as well, A has been selected to be real, symmentric and positive-definite. 

Matrix A and L are printed along with LL^T to show that the decomposition has been done correctly.



Furthermore i have used the cholesky decompositipon af A to solve Ax=b, where b is a random vector. This is done by implementing an algorithm using both forward and 
back substitution. All functions are again found in "funcs.c". 
Solving Ax=b is done using the following; A=LL^T => LL^Tx=b => Ly=b where y=L^Tx. In my program I first solve Ly=b using forward substitution. 
Then L^Tx=y is solved using back substitution. This is all printet when running make and found in the "out.txt"-file as well.

Solving Ax=b using the decomposition does not in fact yield the correct solution as it returns a Ly=/=b as only the first element is the same as in b,
this must be an error in my program, sadly I could not find it. 

