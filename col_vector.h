#ifndef MX2_COL_VECTOR_H
#define MX2_COL_VECTOR_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

template <typename T>
struct col_vector {
    //column vector of values of fixed size

    int n;
    T *values;

    //constructor, destructor, copyconstructor
    col_vector();
    col_vector(int size);
    col_vector(const col_vector &source);
    col_vector(const vector<T> &source);
    ~col_vector();

    //operators
    //subscripting
    T &operator[](int i);
    //binary with constant
    void operator=(const T &k);
    void operator+=(const T &k);
    void operator-=(const T &k);
    void operator*=(const T &k);
    void operator/=(const T &k);
    col_vector operator+(const T &k);
    col_vector operator-(const T &k);
    col_vector operator*(const T &k);
    col_vector operator/(const T &k);
    //binary with col_vector
    void operator=(const col_vector &source);
    void operator+=(const col_vector &source);
    void operator-=(const col_vector &source);
    void operator*=(const col_vector &source);
    void operator/=(const col_vector &source);
    col_vector operator+(const col_vector &source);
    col_vector operator-(const col_vector &source);
    col_vector operator*(const col_vector &source);
    col_vector operator/(const col_vector &source);
    //stream
    friend ostream& operator<<(ostream &output, const col_vector &source) {
        for (int i = 0; i < source.n; ++i) output << source.values[i] << " ";
        return output;
    };

    //methods
    double sum();
    double asum();
    double normSq();

};

#include "col_vector.tpp"

#endif //MX2_COL_VECTOR_H
