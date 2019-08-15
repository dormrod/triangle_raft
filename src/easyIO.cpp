#include "easyIO.h"

//#### writing to file ####
void writeFileDashedLine(ofstream &file, int n){
    //write dashed line
    for (int i=0; i<n;++i) file<<"-";
    file<<endl;
}

void writeFileIndent(ofstream &file, int w){
    //write indent of number of spaces
    for (int i=0; i<w;++i) file<<" ";
}

//#### reading from file ####
void readFileSkipLines(ifstream &file, int nLines){
    //skip over line
    string line;
    for(int i=0; i<nLines;++i) getline(file,line);
}