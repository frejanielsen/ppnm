#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <stdlib.h>
#include <math.h>


void timesJ(gsl_matrix* A, int p, int q, double phi){
	double c=cos(phi), s=sin(phi);
	int n=A->size1;
	double newA_ip=0, newA_iq=0;
	for (int i=0;i<n;i++) {
		newA_ip = c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		newA_iq = s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,newA_ip);
		gsl_matrix_set(A,i,q,newA_iq);
	}
}

void Jtimes(gsl_matrix* A, int p, int q, double phi){
        double c=cos(phi), s=sin(phi);
        int n=A->size2; //only square matricies are relevant for this excersize, so either size could have been used. however picking the sizes as done here generalizes the functions to work for non-square matricies   
        double newA_pi=0, newA_qi=0;
        for (int i=0;i<n;i++) {
                newA_pi = c*gsl_matrix_get(A,p,i)+s*gsl_matrix_get(A,q,i);
                newA_qi = -s*gsl_matrix_get(A,p,i)+c*gsl_matrix_get(A,q,i);
                gsl_matrix_set(A,p,i,newA_pi);
                gsl_matrix_set(A,q,i,newA_qi);
        }
}

void Jt_A_J(gsl_matrix* A, int p, int q, double phi) {
	double c=cos(phi);
	double s=sin(phi);
	double A_ip, A_iq, A_pi, A_qi;
	double A_pq=gsl_matrix_get(A, p, q);
   	double A_pp=gsl_matrix_get(A, p, p);
   	double A_qq=gsl_matrix_get(A, q, q);
   	gsl_matrix_set(A, p, q, 0);
   	gsl_matrix_set(A, p, p, c*c*A_pp - 2*s*c*A_pq + s*s*A_qq);
   	gsl_matrix_set(A, q, q, s*s*A_pp + 2*s*c*A_pq + c*c*A_qq);

 	for (int i=0;i<p;i++){
      		A_ip=gsl_matrix_get(A,i,p);
     		A_iq=gsl_matrix_get(A,i,q);
      		gsl_matrix_set(A,i,p,c*A_ip-s*A_iq);
      		gsl_matrix_set(A,i,q,s*A_ip+c*A_iq);
   	}

   	for (int i=p+1;i<q;i++){
      		A_pi=gsl_matrix_get(A,p,i);
      		A_iq=gsl_matrix_get(A,i,q);
      		gsl_matrix_set(A,p,i,c*A_pi-s*A_iq);
      		gsl_matrix_set(A,i,q,s*A_pi+c*A_iq);
   	}

   	for (int i = q + 1; i < A->size1; ++i){
      		A_pi=gsl_matrix_get(A,p,i);
      		A_qi=gsl_matrix_get(A,q,i);
      		gsl_matrix_set(A,p,i,c*A_pi-s*A_qi);
      		gsl_matrix_set(A,q,i,s*A_pi+c*A_qi);
   	}

}


void jacobi_diag (gsl_matrix* A, gsl_matrix* V){
	int changed=1;
	double A_pq, A_pp, A_qq;
	double phi ,c ,s;
	double newA_pp, newA_qq;
	int n=A->size1; //size of square matrix
	while(changed!=0){
		changed=0;
		for (int p=0;p<n;p++){
			for (int q=p+1; q<n;q++){
				A_pp=gsl_matrix_get(A,p,p);
				A_qq=gsl_matrix_get(A,q,q);
				A_pq=gsl_matrix_get(A,p,q);
				phi=0.5*atan2(2*A_pq,A_qq-A_pp); //calculate tan^-1(2*A_pq/(A_qq-A_pp))
				c=cos(phi);
				s=sin(phi);
				newA_pp=c*c*A_pp-2*s*c*A_pq+s*s*A_qq;
				newA_qq=s*s*A_pp+2*s*c*A_pq+c*c*A_qq;
				if (newA_pp!=A_pp || newA_qq!=A_qq){
					changed=1;
					timesJ(A,p,q,phi);
					Jtimes(A,p,q,-phi);
					timesJ(V,p,q,phi);
				}
			}
		}
	}
}


void jacobi_diag2 (gsl_matrix* A, gsl_matrix* V){
        int changed=1;
        double A_pq, A_pp, A_qq;
        double phi ,c ,s;
        double newA_pp, newA_qq;
        int n=A->size1; //size of square matrix
        while(changed!=0){
                changed=0;
                for (int p=0;p<n-1;p++){
                        for (int q=p+1; q<n;q++){

                                A_pp=gsl_matrix_get(A,p,p);
                                A_qq=gsl_matrix_get(A,q,q);
                                A_pq=gsl_matrix_get(A,p,q);
                                phi=0.5*atan2(2*A_pq,A_qq-A_pp); //calculate tan^-1(2*A_pq/(A_qq-A_pp))
                                c=cos(phi);
                                s=sin(phi);
                                newA_pp=c*c*A_pp-2*s*c*A_pq+s*s*A_qq;
                                newA_qq=s*s*A_pp+2*s*c*A_pq+c*c*A_qq;
                                if (newA_pp!=A_pp || newA_qq!=A_qq){
                                        changed=1;
                                        Jt_A_J(A,p,q,phi);
                                        timesJ(V,p,q,phi);
                                }
                        }
                }
        }
}







void matrixprint(gsl_matrix* A){
int n=A->size1;


for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
		printf("%10.2g",gsl_matrix_get(A,i,j));
	}
	printf("\n");
}

printf("\n\n");
}



