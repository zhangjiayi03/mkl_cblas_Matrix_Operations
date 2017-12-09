#define Max(a,b) (a>b)?a:b
#define Min(a,b) (a<b)?a:b
#include "stdio.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "stdafx.h"
#include "mkl.h"
#include "mkl_cblas.h"
//#define n 3
#define dim_x0 4
#define dim_y0 3
#define dim_x1 3
#define dim_y1 3
#define matrix_order LAPACK_ROW_MAJOR
//#include "mkl_lapacke.h"
int Matrix_eig(double S[dim_x0*dim_y0]);
int Matrix_inv(double S[dim_x0*dim_y0]);
int Matrix_qr(double S[dim_x0*dim_y0]);
int Matrix_mul(double A[dim_x0*dim_y0], double B[dim_x1*dim_y1]);


int main()
{
	FILE *file0,*file1;//如果不是VS中编译，FILE *file0 = fopen("d:\\test.txt","r");
	errno_t err= fopen_s(&file0,"d:\\test.txt","r");
	errno_t err1= fopen_s(&file1,"d:\\test1.txt","r");

	int m,k,re_eig,re_inv;
	int re_qr;
	double S[(dim_x0+1)*dim_y0],S1[(dim_x1+1)*dim_y1];
	printf("源点阵S\n");
	for(m=0;m<dim_x0;m++)
	{
        for(k=0;k<dim_y0;k++)
          {
              fscanf_s(file0,"%lf",&S[k+m*dim_x0]);//不是VS中则 fscanf(file0,"%lf",&S[p+m*n]);
			  printf("%lf\t",S[k+m*dim_x0]);
          }
		printf("\n");
	}
	printf("\n");

	printf("目标点阵S1\n");
	for(m=0;m<dim_x1;m++)
	{
        for(k=0;k<dim_y1;k++)
          {
              fscanf_s(file1,"%lf",&S1[k+m*dim_x1]);//不是VS中则 fscanf(file0,"%lf",&S[p+m*n]);
			  printf("%lf\t",S1[k+m*dim_x1]);
          }
		printf("\n");
	}
	printf("\n");
	//调用特征向量函数
	//re_eig = Matrix_eig(S1);
	//if (re_eig == 0) { printf("特征向量求解成功\n");}


	//re_inv = Matrix_inv(S);
	//if (re_inv == 0) { printf("逆矩阵求解成功\n");}
	//getchar();//important
        
	//调用QR分解
	//re_qr = Matrix_qr(S);
	//if(re_qr == 0){ printf("QR分解成功\n");}
	
	//调用乘法 
	//Matrix_mul(S,S1);*/



	getchar();//important
	return 0;
}


int Matrix_eig(double S[dim_x0*dim_x0])
{
	//int matrix_order = LAPACK_ROW_MAJOR;
    char jobvl = 'N';
    char jobvr = 'V';
	int lda = dim_x0;
    double wr[dim_x0] = {0};
    double wi[dim_x0] = {0};
    double vl[dim_x0*dim_x0];
    int ldvl = dim_x0;
    double vr[dim_x0*dim_x0];
    int ldvr = dim_x0;
	int info;
	clock_t start,end;
	start = clock( );
	info = LAPACKE_dgeev(matrix_order,jobvl,jobvr,dim_x0,S,lda,wr,wi,vl,ldvl,vr,ldvr);	
	end = clock();
	printf("Time: %lf\n\n", (double) (end - start) / CLOCKS_PER_SEC );
	if(info==0){
        int i = 0;
        int j = 0;
		int flag = 0;//区分复特征值的顺序
        for(i=0;i<dim_x0;i++){
            printf("eigenvalue %d:\t",i+1);
            printf("%.6g + %.6gi\t",wr[i],wi[i]);
			printf("\n\n");
            printf("right eigenvector:\t ");
			if(wi[i]==0)
			{
					for(j=0;j<dim_y0;j++){
					printf("%.6g\t",vr[i+j*dim_x0]);
					//p[i] = vr[i*n+j];
				}
			}
			else if(flag==0)
			{
				flag=1;
				for(j=0;j<dim_y0;j++)
				{
					printf("%.6g + %.6gi\t",vr[i+j*dim_x0],vr[(j+1)*dim_x0+i]);
				}
			}
			else if(flag==1)
			{
				flag=0;
				for(j=0;j<dim_y0;j++)
				{
					printf("%.6g + %.6gi\t",vr[(j-1)*dim_x0+i],-vr[j*dim_x0+i]);
				}
			}
            printf("\n\n");
		}
	}
	return 0;
}

int Matrix_inv(double S[dim_x0*dim_y0])
{
	int i,j;
	//int matrix_order = LAPACK_ROW_MAJOR;
	int lda = dim_x0;
	int info;
	int ipiv[Max(1,(Min(dim_x0,dim_y0)))];
	char trans = 'N';
	clock_t start,end;
	start = clock( );
	info = LAPACKE_dgetrf(matrix_order,dim_x0,dim_y0,S,lda,ipiv);
	end = clock();
	printf("Time: %lf\n\n", (double) (end - start) / CLOCKS_PER_SEC );
	if (info == 0)
	{
		double I[dim_x0*dim_y0]={0};
		int ldb = lda;
		for(i = 0;i<dim_x0;i++)
		{
				I[i*(dim_x0+1)] = 1;
		}
		LAPACKE_dgetrs(matrix_order,trans,dim_x0,dim_x0,S,lda,ipiv,I,ldb);

        for(i=0;i<dim_x0;i++){
			for(j=0;j<dim_y0;j++){
					printf("%.6g\t",I[i*dim_x0+j]);	
				}
			printf("\n");
		}
            printf("\n\n");
	}
	return 0;
}

int Matrix_qr(double S[dim_x0*dim_y0])
{
	int i,j;
	//int matrix_order = LAPACK_ROW_MAJOR;
	int lda = dim_x0;
	int info;
	double tau[Max(1,(Min(dim_x0,dim_y0)))] = {0};
	double A[dim_x0*dim_x0] = {0};
	clock_t start,end;
	/*lapack_int LAPACKE_dgeqrf( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, double* tau );*/
	
	info = LAPACKE_dgeqrf(matrix_order, dim_x0, dim_y0, S, lda, tau);

	if (info == 0)
	{
		 for(i=0;i<dim_x0;i++){
			for(j=0;j<dim_y0;j++){
					A[i*dim_x0+j] = S[i*dim_x0+j];
					printf("%.6g\t",S[i*dim_x0+j]);	
				}
			printf("\n");
		}
            printf("\n\n");

		start = clock( );
		/*lapack_int LAPACKE_dorgqr(int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int k, double* a, lapack_int lda,
                           const double* tau );*/
		LAPACKE_dorgqr(matrix_order,dim_x0,dim_x0,dim_y0,A,lda,tau);
		end = clock();
		printf("Time: %lf\n\n", (double) (end - start) / CLOCKS_PER_SEC );
        for(i=0;i<dim_x0;i++){
			for(j=0;j<dim_y0;j++){
					printf("%.6g\t",A[i*dim_x0+j]);	
				}
			printf("\n");
		}
            printf("\n\n");
	}
	return 0;
}

int Matrix_mul(double A[dim_x0*dim_y0], double B[dim_x1*dim_y1])
{
	/*cblas_dgemm(const  CBLAS_LAYOUT Layout, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                 const MKL_INT K, const double alpha, const double *A,
                 const MKL_INT lda, const double *B, const MKL_INT ldb,
                 const double beta, double *C, const MKL_INT ldc);*/
	double C[(dim_x0+1)*dim_y1];
	//const char TransA = 'N';
	//const char TransB = 'N';
	double alpha = 1.0;
	double beta = 0.0;
	int lda = Max(1,dim_x0);
	int ldb = Max(1,dim_y0);
	int ldc = Max(1,dim_x0);
	int i,j;
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_x0, dim_y1, dim_y0, alpha, A, lda, B, ldb, beta, C, ldc);
	//printf("Time: %lf\n\n", (double) (end - start) / CLOCKS_PER_SEC );
    for(i=0;i<dim_x0;i++){
		for(j=0;j<dim_y1;j++){
				printf("%.6g\t",C[i*dim_x0+j]);	
			}
			printf("\n");
		}
            printf("\n\n");
	return 0;
}