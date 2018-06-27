//Easier writing to console and reading/writing to files

#ifndef MX2_EASYIO_H
#define MX2_EASYIO_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

using namespace std;

//#### writing to console ####
template <typename T>
void consoleValue(T value){
    //write single value to console
    cout<<value<<endl;
}
template <typename T>
void consoleVector(vector<T> values){
    //write vector of values to console
    for(int i=0; i<values.size(); ++i) cout<<values[i]<<" ";
    cout<<endl;
}
template <typename T>
void consoleMatrix(vector< vector<T> > values){
    //write matrix of values to console
    for(int i=0; i<values.size(); ++i){
        for(int j=0; j<values[i].size(); ++j) cout<<values[i][j]<<" ";
        cout<<endl;
    }
}
template <typename T>
void consoleArray(T values, int n){
    //write array of values to console
    for(int i=0; i<n; ++i) cout<<values[i]<<" ";
    cout<<endl;
}

//#### writing to file ####
void writeFileDashedLine(ofstream &file, int n=70);
void writeFileIndent(ofstream &file, int w=4);

template <typename T>
void writeFileValue(ofstream &file, T value){
    //write single value to file
    file<<value<<endl;
}
template <typename S, typename T>
void writeFileValue(ofstream &file, S value0, T value1, int width=20){
    //write two values to file
    file<<setw(width)<<left<<value0<<setw(width)<<left<<value1<<endl;
}
template <typename S, typename T, typename U>
void writeFileValue(ofstream &file, S value0, T value1, U value2, int width=20){
    //write three values to file
    file<<setw(width)<<left<<value0<<setw(width)<<left<<value1<<setw(width)<<left<<value2<<endl;
}
template <typename S, typename T, typename U, typename V>
void writeFileValue(ofstream &file, S value0, T value1, U value2, V value3, int width=20){
    //write four values to file
    file<<setw(width)<<left<<value0<<setw(width)<<left<<value1<<setw(width)<<left<<value2<<setw(width)<<left<<value3<<endl;
}
template <typename T>
void writeFileVector(ofstream &file, vector<T> values, int width=10){
    //write vector of values to file
    for(int i=0; i<values.size(); ++i){
        file<<setw(width)<<left<<values[i];
    }
    file<<endl;
}
template <typename T>
void writeFileArray(ofstream &file, T values, int n, int width=10){
    //write array of values to file
    for(int i=0; i<n; ++i){
        file<<setw(width)<<left<<values[i];
    }
    file<<endl;
}

//#### reading from file ####
void readFileSkipLines(ifstream &file, int nLines=1);

template <typename T>
void readFileValue(ifstream &file, T &value){
    //read value from first column
    string line;
    getline(file,line);
    istringstream ss(line);
    ss >> value;
}
template <typename T>
void readFileRowVector(ifstream &file, vector<T> &vec, int nCols){
    //read in vector of n columns
    vec.clear();
    string line;
    getline(file,line);
    istringstream ss(line);
    T val;
    for (int i=0; i<nCols; ++i){
        ss >> val;
        vec.push_back(val);
    }
}
template <typename T>
void readFileAll(ifstream &file, vector< vector<T> > &all){
    all.clear();
    string line;
    T value;
    vector<T> row;
    while(getline(file, line)){
        istringstream ss(line);
        row.clear();
        while(ss>>value) row.push_back(value);
        all.push_back(row);
    }
}

#endif //MX2_EASYIO_H
