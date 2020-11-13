#ifndef MATRIX
#define MATRIX
/**
 Header file for a generic matrix class
 implemeted with arrays
*/
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <ostream>

void* Malloc(int memsize);

void* Calloc(int num, int size);

void* Realloc(void* old, int size);

template <class numeric>
class matrix;

//template<class numeric>
//std::istream operator>> (std::istream& i, matrix<numeric>& M);

template <class numeric>
class matrix
{
public:
    friend std::ostream& operator<< (std::ostream& o, const matrix<numeric>& M)
    {
        o << "[\n";
        for (int count = 0; count < M.row * M.col; ++count)
        {
            o << M.data[count];
            if ((count + 1) % M.col)
                o << ", ";
            else
                o << "\n";
        }
        o << "]";
        return o;
    }
    
    const matrix& operator=(const matrix & old);

    // friend istream operator>> (istream &i, const matrix<numeric>& M);
    matrix(void)
    {
        row = 0;
        col = 0;
        data = nullptr;
    }
    matrix(int rows, int columns, numeric* Array)
    {
        row = rows;
        col = columns;
        data = (numeric*)Calloc(row * col, sizeof(numeric));
        if (Array) SetData(Array, row * col);
    };
    template <class scalar>
    matrix(scalar scale)
    {
        row = 1;
        col = 1;
        data = (numeric*)Calloc(row * col, sizeof(numeric));
        data[0] = scale;
    }

    template<class iterable>
    matrix(int rows, int columns, iterable& datalist)
    {
        row = rows;
        col = columns;
        data = (numeric*) Calloc(row * col, sizeof(numeric));
        int location = 0;
        for (auto i = datalist.begin(); i != datalist.end() && location < row * col; ++i)
        {
            *(data + location) = *i;
            ++location;
        }

    }
    matrix(int r, int c, ...);

    
    matrix(const matrix& copy)
    {
        row = copy.row;
        col = copy.col;
        data = (numeric*)Malloc(col * row * sizeof(numeric));

        std::memcpy(data, copy.data, row * col * sizeof(numeric));
    }


    matrix operator[](int row);
    matrix operator()(int col);

    matrix operator+ (matrix& add)
    {
        matrix ret(row, col, nullptr);
        if (add.col == col && add.row == row)
        {
            for (int i = 0; i < row * col; ++i)
                ret.data[i] = *(data + i) + *(add.data + i);
        }
        else
            std::cout << "\nWrong dimensions\n";
        return ret;
    }

    matrix operator- (matrix& sub)
    {
        matrix ret(row, col, nullptr);
        if (sub.col == col && sub.row == row)
        {
            for (int i = 0; i < row * col; ++i)
                *(ret.data + i) = *(data + i) - *(sub.data + i);
        }
        else
            std::cout << "\nWrong dimensions\n";
        return ret;

    }

    matrix operator- (int row)
    {
        return *row_iterator(this, row);
    }
    matrix operator|(int col)
    {
        return *col_iterator(this, col);
    }

    // allow 1 by 1 to be treated as scalar
    bool IsScalar(void) { return (row == 1 && col == 1); }
    
    bool IsVector(void) { return (row == 1 || col == 1); }

    matrix T(matrix& m);

    void print(void) { std::cout << *this; };
    
    int size(void) { return row * col; }
    
    

    void SetData(numeric* Data, int size) { std::memcpy(data, Data, sizeof(numeric) * size); };
   
    class col_iterator
    {
    public:
        col_iterator(matrix* Matrix, int col = 0)
        {
            m = Matrix;
            column = col;
        }
        const matrix<numeric> operator*(void)
        {
            matrix ret(m->rows(), 1, nullptr);
            
            for (int i = 0; i < m->row; ++i)
                *(ret.data+i) = *(m->data + (column + i * m -> row));
            
            return ret;
        }

        const col_iterator& operator++(void)
        {
            ++column;
            return *this;
        }
        const col_iterator& operator--(void)
        {
            --column;
            return *this;
        }
        bool operator ==(const col_iterator& i)
        {
            return i.m == m && i.column == column;
        }
        bool operator !=(const col_iterator& i)
        {
            return i.m != m || i.column != column;
        }


    private:
        matrix* m;
        int column;
    };

    class row_iterator
    {
    public:
        row_iterator(matrix* Matrix, int r = 0)
        {
            m = Matrix;
            row = r;
        }
        const matrix<numeric> operator*(void)
        {
            return matrix(1,m->col, m->data + row*m->col);
        }

        const row_iterator& operator++(void)
        {
            ++row;
            return *this;
        }
        const row_iterator& operator--(void)
        {
            --row;
            return *this;
        }
        bool operator ==(const row_iterator& i)
        {
            return i.m == m && row == i.row;
        }
        bool operator !=(const row_iterator& i)
        {
            return i.m != m || i.row != row;
        }

    private:
        matrix* m;
        int row;
    };

    class iterator
    {
    public:
        iterator(matrix* Matrix, int location = 0)
        {
            p = Matrix->data+location;
        }

        const numeric operator*(void)
        {
            return *p;
        }

        const iterator& operator++(void)
        {
            ++p;
            return *this;
        }
        const iterator& operator--(void)
        {
            --p;
            return *this;
        }
        bool operator ==(const iterator& i)
        {
            return i.p == p;
        }
        bool operator !=(const iterator& i)
        {
            return i.p != p;
        }

    private:
        numeric* p;
    };

    iterator begin(void) { return iterator(this, 0); };
    iterator end(void) { return iterator(this, row*col); };

    row_iterator row_begin(void) { return row_iterator(this , 0); }
    row_iterator row_end(void) { return row_iterator(this, row); }

    col_iterator col_begin(void) { return col_iterator(this, 0); }
    col_iterator col_end(void) { return col_iterator(this, col); }
    
    int rows(void) { return row; }
    int cols(void) { return col; }

    ~matrix() { std::cout << "Destructor: \n" << (*this); free(data); };

private:
    int row;
    int col;
    numeric* data;
};


#endif