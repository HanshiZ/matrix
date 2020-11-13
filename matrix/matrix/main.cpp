#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vector>
#include "matrix.h"
/*
void * Malloc(int memsize);

struct Matrix
{

public:
    Matrix(int r=0, int c=0);
    Matrix(int r, int c, double array[]);
    Matrix(int r, int c, int size, ...);

    void SetData(double * data,int size){ memcpy(matrix,data,sizeof(double)*size);};
    ~Matrix();
private:
    double * matrix;
    int row;
    int col;
};

const int sd = sizeof(double);

struct Matrix mcreate(int row, int col);
short eq(struct Matrix *a, struct Matrix *b);
struct Matrix madd(struct Matrix *a, struct Matrix *b);
struct Matrix scmult(double scale, const struct Matrix m);
struct Matrix Ma(int r, int c, double array[]);
struct Matrix M(int r, int c, int size, ...);
struct Matrix I(int n);
struct Matrix T(struct Matrix *m);
struct Matrix rswap(struct Matrix *m, int r1, int r2);
void aradd(double a[], double add, int beg, int end);
void armult(double a[], double prod, int beg, int end);
void addar(double a[], double b[], int len, int pos, int scale);
struct Matrix mmult(struct Matrix *a, struct Matrix *b);
double * r(struct Matrix *m, double * p, int r);
void *del(struct Matrix *m);
double * c(struct Matrix *m, double *p, int c);
struct Matrix cpy(const struct Matrix const *m);

struct Matrix inv(struct Matrix *m);

void ptm(const struct Matrix m);
*/
int main()
{
    int test[4] = { 2,3,0,4 };
    matrix<int> M = matrix<int>(1, 4, test);
    std::vector<double> x = {2.3,6.8,8,0.34 };
    matrix <double> M2 = matrix<double>(2, 2, M);
    matrix <double> M3 = matrix <double>(2, 2, x);

    auto m4 = M2 | 1;
    auto m5 = (M2 | 0);
    std:: cout << m4 - m5;
    return 0;

}
/*
void *Calloc(int num, int size)
{
    void * temp = calloc(num,size);
    if (temp)
        return temp;
    else
        printf("***Warning:Memory allocation failed!***");
    return 0;
}

void * Malloc(int memsize)
{
    void * temp = malloc(memsize);
    if (temp) return temp;
    printf("***Warning:Memory allocation failed!***");
    return 0;
}

Matrix ::Matrix(int r, int c)
{
  row =r; col =c;
  matrix = (double*) Calloc(row*col, sizeof(double));

}

bool eq(struct Matrix *a, struct Matrix *b)
{   if (a->row != b->row || b->col != a->col) {return 0;}
    for (int i =0; i < a->row*a->col; i++)
    {if (a->matrix[i]!=b->matrix[i]) {return 0;}}
    return 1;
}

struct Matrix madd(struct Matrix *a, struct Matrix *b)
{  struct Matrix c= mcreate(a->row, a->col);
    if (a->row != b->row || a->col != b->col)
    { printf ("\n***Invalid operation***\n dimensions != \n");
       c.row= c.col = 0;
      return (c);
    }
  double temp[(a->row)*(a->col)];
  for (int i =0; i < a->row*a->col; i++)
    {  temp[i] = a->matrix[i]+b->matrix[i];}

   memcpy (c.matrix,temp, sd*a->row*a->col);
   return c;
}

struct Matrix scmult(double scale, const struct Matrix m)
{   double a[m.row*m.col];
    struct Matrix otput = mcreate(m.row, m.col);

    for (int i = 0; i < m.row*m.col; i++)
            { a[i]= scale * m.matrix[i];}
   memcpy(otput.matrix, a, sd*m.row*m.col);

   return otput;
}

Matrix::Matrix(int r, int c, double array[])
{
    Matrix::Matrix(r,c);

    SetData(array, r*c);
}

Matrix::Matrix(int r, int c, int size, ...)
{   matrix = (double*) Malloc(size*sizeof(double));
    double entry;
    va_list alist;


    va_start(alist, size);
    for (int i= 0; i< size; i++)
    {  entry = va_arg(alist, double);
       matrix[i] = entry;
     }
    va_end(alist);
    row =r;
    col =c;
    return  otput;
}
struct Matrix I(int n)
{   int count =0;
    double * a = (double *) Malloc(sd* n*n);
    struct Matrix otput =mcreate(n,n);
    for (int i = 0; i < n*n ; i++)
        {  if (count) { a[i] = 0; }
           else {a[i]= 1;}
           count = count+1;
           if (!(count % (n+1))) {count =0;}
        }

    memcpy(otput.matrix, a, sd*n*n);
    free(a); a=0;
    return otput;
}


void ptm(const struct Matrix m)
{
    for (int i=0; i < m.row*m.col; i++)
    {  if (i % m.col ==0) {printf("\n");}
        if (!i ) {printf("[ %.3f,", m.matrix[i]);}
        else if ((i+1) % m.col){printf("  %.3f,", m.matrix[i]);}
        else if (i == m.col*m.row-1) printf("  %.3f]", m.matrix[i]);
        else {printf("  %.3f,", m.matrix[i]);}

    }
}


struct Matrix rswap(struct Matrix *m, int r1, int r2)
{   r1 = r1-1;
    r2=r2-1;
    double * temp = (double *) Malloc(sd*m->col);
    memcpy(temp, &(m->matrix[r1*m->col]), m->col*sd);
    memmove(&(m->matrix[r1*m->col]), &(m->matrix[r2*m->col]), m->col * sd);
    memcpy(&(m->matrix[r2*m->col]),temp, m->col*sd);
    free(temp); temp = 0;
    return *m;
}

struct Matrix cpy(const struct Matrix const *m)
{
    struct Matrix c = mcreate(m->row, m->col);
    memcpy(c.matrix, m->matrix,sd*(m->row)*m->col);
    return c;
}

double * r(struct Matrix *m, double *p, int r)
{
    memcpy(p, &(m->matrix[(r-1)*m->col]), sd*m->col);
    return p;
}

double * c(struct Matrix *m, double *p, int c)
{
    int count =0;
    for (int i= c-1; i < m->row*m ->col; i=i+m ->col)
        {  p[count]=m->matrix[i];

            count = count +1;
        }
    return p;
}

void aradd(double a[], double add, int beg, int end)
{
    for(int i=beg; i < end+1; i++) {a[i]= a[i]+add;}
}

void armult(double a[], double prod, int beg, int end)
{
    for(int i=beg; i < end+1; i++) {a[i]= a[i]*prod;}
}

struct Matrix mmult(struct Matrix *a, struct Matrix *b)
{   if (a->col != b->row)
        {
            printf("\n***Invalid operation***\n col != row \n");
            return *a;
        }
    struct Matrix otput = mcreate(a->row,b->col);
    for (int i =0; i < a->row; i++)
        {
            for (int j =0; j < b->col; j++)
                {  double s =0;
                   for (int k =0; k < a->col; k++)
                    { s = s + (a->matrix[a->col*i+k]*b->matrix[j + k * b->col]);}
                    otput.matrix[b->col*i+j] =s;
                }

        }

    return otput;
}

struct Matrix inv(struct Matrix *m)
{  if (m->row != m->col)
        {  printf("Invalid operation");
           return mcreate(m->row, m->col);
        }
   struct Matrix i = I(m->row);
   for (int j = 0; j < m->col-1 ; j++)
        {  for(int k=j; k < m->row-1; k++)
              {  	int add =1;
        if (m->matrix[j+j*m->col])
                 {
                     if(m->matrix[(j+j*m->col)+(add)*m->col] && (j+j*m->col)+(add)*m->col < m->col*m->row-1) /*(k off)*/
                     /*              {  const double r = (m->matrix[j+j*m->col])/(m->matrix[(j+j*m->col)+(add)*m->col]);
                                      printf("%.2f",r);
                                      armult(&(m->matrix[(j+j*m->col)+(add)*m->col]), -r,
                                      0,m->col-j);
                                      armult(&(i.matrix[(add+j)*m->col]),-r,
                                      0,m->col-j);
                                      addar(&m->matrix[(j+j*m->col)+(add)*m->col], &m->matrix[(j+j*m->col)],m->row-j,0,1);
                                      addar(&i.matrix[j*m->col + (add)*m->col], &i.matrix[j*m->col] ,j+1,0,1);
                                   }
                                   else {continue;} add = add+1;
                                }
                                else {rswap(m, k+1,k+2); rswap(&i, k+1,k+2);}
                               }
                            continue;
                       }
                  int loc = 1;
                  for (int l = m->row * m->col-1; l >= 0 ; l = l- (m->col+1))
                       { if (m->matrix[l])
                           {
                               armult(&i.matrix[l - i.col + loc-1] , 1/m->matrix[l], 0, m->col );
                               for (int o=l; o > i.col; o = o - i.col)
                                   {addar(&i.matrix[o-2*i.col+ loc], &i.matrix[o-i.col+loc], i.col,0,-(m->matrix[o-i.col]));}
                           }
                         else {printf ("Not invertible"); return *m;}
                         loc =loc+1;
                       }

                  return i;
               }

               void addar(double a[], double b[], int len, int pos, int scale)
               {
                   for (int i =0; i< len; i++) {a[i]=b[i]*pow(-1, pos)*scale+a[i];}
               }

               struct Matrix T(struct Matrix *m)
               {  struct Matrix t = mcreate(m->col, m->row);
                  for (int i =1; i<=m->col; i++)
                   {   double temp[m->col];
                       memcpy(&(t.matrix[(i-1)*m->row] ), c(m, temp ,i),sd*m->row); }
                  return t;
               }

               Matrix :: ~Matrix()
               {
                    free(matrix);
               }
               */