#include "matrix.h"
#include "matrix.h"
#include <cstdlib>
#include <cctype>
#include <iostream>
using std::ostream;
using std::istream;


#define TEM template<class numeric>

void* Calloc(int num, int size)
{
    void* temp = calloc(num, size);
    if (temp)
        return temp;
    else
        printf("***Warning:Memory allocation failed!***");
    return 0;
}

void* Malloc(int memsize)
{
    void* temp = malloc(memsize);
    if (temp) return temp;
    printf("***Warning:Memory allocation failed!***");
    return 0;
}

void* Realloc(void* old, int size)
{
    void* temp = realloc(old, size);
    if (temp)
        return temp;
    else
        printf("***Warning:Memory allocation failed! Old array returned!***");
    return old;
}


template <class numeric>
matrix<numeric>::matrix(int r, int c, ...)
{
    row = r;
    col = c;
    int size = r * c;
    data = (numeric*)Calloc(size, sizeof(numeric));
    numeric entry;
    va_list vl;
    va_start(vl, size);
    for (int i = 0; i < r * c; i++)
    {
        entry = va_arg(vl, numeric);
        *(data + i) = entry;
    }
    va_end(vl);

}


TEM
matrix<numeric> matrix<numeric>::T(matrix<numeric>& m)
{
    matrix<numeric> t = matrix<numeric>(m.col, m.row);
    /*for (int i = 1; i <= m->col; i++)
    {
        double temp[m->col];
        memcpy((t.matrix[(i - 1) * m->row]), c(m, temp, i), sd * m->row);
    }*/
    return t;
}

TEM 
const matrix<numeric>& matrix<numeric>::operator=(const matrix& old)
{
    if (&old != this && row == old.row && col == old.col)
    {
        memcpy(data,old.data,row*col*sizeof(numeric));
    }
    else if (&old == this)
    {
        std::cout << "Attempted self assignment.";
    }
    else
    {
        std::cout << "Wrong dimensions.";
    }
    return *this;
}
/*
TEM
matrix<numeric> matrix<numeric> ::operator+ (const matrix& add)
{
    matrix ret;
    if (add.col == col && add.row == row)
    {
        ret.col = col;  
        ret.row = row;
        numeric* dat =(numeric*) Malloc(sizeof(numeric)*row*col);

        for (matrix<numeric>::iterator i = add.begin(); i != add.end(); ++i)
            *(dat++) = *(data++) + *i;
        ret.data = dat;
    }
    else
        std::cout << "\nWrong dimensions\n";
    return ret;
}*/
/*TEM
matrix<numeric> matrix<numeric>:: operator- (const matrix& sub)
{
    matrix ret(row, col, nullptr);
    if (sub.col == col && sub.row == row)
    {
        for (int i = 0; i < row*col; ++i)
            *(ret.data + i) = *(data + i) - *(sub.data + i);
    }
    else
        std::cout << "\nWrong dimensions\n";
    return ret;

}*/
/**TEM
matrix<numeric> matrix<numeric>::operator[](int row)
{
    if (IsVector())
    {
        return
    }
    matrix ret = matrix(1, col);
    ret.SetData((row-1)*col,col);
    return ret;
}

TEM
matrix<numeric> matrix<numeric>::operator()(int col)
{
    matrix ret = matrix(row, 1);
    
    
    for (int i = c - 1; i < m->row * m->col; i = i + m->col)
    {
        p[count] = m->matrix[i];

        count = count + 1;
    }
    return p;
}
*/