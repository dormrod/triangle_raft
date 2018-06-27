//logs messages from code

#ifndef MX2_LOGFILE_H
#define MX2_LOGFILE_H

#include <iostream>
#include <ctime>
#include <chrono>
#include <vector>
#include "easyIO.h"

using namespace std;

class Logfile {
    //class for logging messages from calculation

private:
    //variables
    string filename;
    chrono::high_resolution_clock::time_point time0, time1;

    //functions
    void initialise();

public:

    //constructors
    Logfile();
    Logfile(string name);
    ~Logfile();

    //logging functions
    string currentTime();
    string currentTimeAndDate();
    int timeElapsedSeconds();
    int timeElapsedMinutes();
    void errorlog(string message, string errorType, int errorCode=1);
    template <typename T>
    void log(string message, T value, string time, int indent, bool linebreak){
        ofstream file(filename,ios::in|ios::app);
        for(int i=0; i<indent; ++i) writeFileIndent(file);
        if(time=="min"){
            int dt=timeElapsedMinutes();
            string timeElapsed=to_string(dt)+" minutes";
            writeFileValue(file,message,value,timeElapsed,1);
        }
        else if(time=="sec"){
            int dt=timeElapsedSeconds();
            string timeElapsed=to_string(dt)+" seconds";
            writeFileValue(file,message,value,timeElapsed,1);
        }
        else writeFileValue(file,message,value);
        if(linebreak) writeFileDashedLine(file);
        file.close();
    };
    template <typename T>
    void log(vector<T> values, int indent, bool linebreak){
        ofstream file(filename,ios::in|ios::app);
        for(int i=0; i<indent; ++i) writeFileIndent(file);
        writeFileVector(file,values);
        if(linebreak) writeFileDashedLine(file);
        file.close();
    };
};


#endif //MX2_LOGFILE_H
