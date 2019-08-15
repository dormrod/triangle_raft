#include "col_vector.h"

//##### CONSTRUCTORS, DESTRUCTORS #####
template <typename T>
col_vector<T>::col_vector() {
    //default of size 1
    n=1;
    values=new T[n]();
}

template <typename T>
col_vector<T>::col_vector(int size) {
    //use given size
    n=size;
    values=new T[n]();
}

template <typename T>
col_vector<T>::col_vector(const col_vector &source) {
    //make deep copy of values
    n=source.n;
    values=new T[n]();
    for(int i=0; i<n; ++i) values[i]=source.values[i];
}

template <typename T>
col_vector<T>::col_vector(const vector<T> &source) {
    //make copy from STL vector
    n=source.size();
    values=new T[n]();
    for(int i=0; i<n; ++i) values[i]=source[i];
}

template <typename T>
col_vector<T>::~col_vector() {
    //clear allocated memory
    delete[] values;
}

//##### SUBSCRIPTING OPERATORS
template <typename T>
T& col_vector<T>::operator[](int i) {
    return values[i];
}

//##### BINARY OPERATORS WITH CONSTANT #####
template <typename T>
void col_vector<T>::operator=(const T &k) {
    for(int i=0; i<n; ++i) values[i]=k;
}

template <typename T>
void col_vector<T>::operator+=(const T &k) {
    for(int i=0; i<n; ++i) this->values[i]+=k;
}

template <typename T>
void col_vector<T>::operator-=(const T &k) {
    for(int i=0; i<n; ++i) this->values[i]-=k;
}

template <typename T>
void col_vector<T>::operator*=(const T &k) {
    for(int i=0; i<n; ++i) this->values[i]*=k;
}

template <typename T>
void col_vector<T>::operator/=(const T &k) {
    for(int i=0; i<n; ++i) this->values[i]/=k;
}

template <typename T>
col_vector<T> col_vector<T>::operator+(const T &k) {
    col_vector v(this->n);
    for(int i=0; i<n; ++i) v[i]=this->values[i]+k;
    return v;
}

template <typename T>
col_vector<T> col_vector<T>::operator-(const T &k) {
    col_vector v(this->n);
    for(int i=0; i<n; ++i) v[i]=this->values[i]-k;
    return v;
}

template <typename T>
col_vector<T> col_vector<T>::operator*(const T &k) {
    col_vector v(this->n);
    for(int i=0; i<n; ++i) v[i]=this->values[i]*k;
    return v;
}

template <typename T>
col_vector<T> col_vector<T>::operator/(const T &k) {
    col_vector v(this->n);
    for(int i=0; i<n; ++i) v[i]=this->values[i]/k;
    return v;
}

//##### BINARY OPERATORS WITH COL VECTOR #####
template <typename T>
void col_vector<T>::operator=(const col_vector &source) {
    if (this == &source) return;
    if(n!=source.n){
        delete[] values;
        n=source.n;
        values=new T[n]();
    }
    for(int i=0; i<n; ++i) this->values[i]=source.values[i];
    return;
}

template <typename T>
void col_vector<T>::operator+=(const col_vector &source) {
    for(int i=0; i<n; ++i) this->values[i]+=source.values[i];
}

template <typename T>
void col_vector<T>::operator-=(const col_vector &source) {
    for(int i=0; i<n; ++i) this->values[i]-=source.values[i];
}

template <typename T>
void col_vector<T>::operator*=(const col_vector &source) {
    for(int i=0; i<n; ++i) this->values[i]*=source.values[i];
}

template <typename T>
void col_vector<T>::operator/=(const col_vector &source) {
    for(int i=0; i<n; ++i) this->values[i]/=source.values[i];
}

template <typename T>
col_vector<T> col_vector<T>::operator+(const col_vector &source) {
    col_vector v(this->n);
    for(int i=0; i<n; ++i) v[i]=this->values[i]+source.values[i];
    return v;
}

template <typename T>
col_vector<T> col_vector<T>::operator-(const col_vector &source) {
    col_vector v(this->n);
    for(int i=0; i<n; ++i) v[i]=this->values[i]-source.values[i];
    return v;
}

template <typename T>
col_vector<T> col_vector<T>::operator*(const col_vector &source) {
    col_vector v(this->n);
    for(int i=0; i<n; ++i) v[i]=this->values[i]*source.values[i];
    return v;
}

template <typename T>
col_vector<T> col_vector<T>::operator/(const col_vector &source) {
    col_vector v(this->n);
    for(int i=0; i<n; ++i) v[i]=this->values[i]/source.values[i];
    return v;
}

//##### METHODS #####
template <typename T>
double col_vector<T>::sum() {
    //sum of all values
    double sum=0.0;
    for(int i=0; i<n; ++i) sum+=values[i];
    return sum;
}

template <typename T>
double col_vector<T>::asum() {
    //sum of all absolute values
    double sum=0.0;
    for(int i=0; i<n; ++i) sum+=fabs(values[i]);
    return sum;
}

template <typename T>
double col_vector<T>::normSq(){
    //norm square
    double nSq=0.0;
    for(int i=0; i<n; ++i) nSq+=values[i]*values[i];
    return nSq;
}
