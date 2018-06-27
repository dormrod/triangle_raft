#include "logfile.h"

Logfile::Logfile(){
    //default constructor
    filename="logfile.log";
    initialise();
}

Logfile::Logfile(string name){
    //construct with filename
    filename=name+".log";
    initialise();
}

Logfile::~Logfile() {
    //close with message
    ofstream file(filename,ios::in|ios::app);
    string timeAndDate=currentTimeAndDate();
    string message="Simulation complete: "+timeAndDate;
    writeFileValue(file,message);
    file.close();
}

string Logfile::currentTimeAndDate() {
    //get current time and date
    time_t now = time(0);
    char* dt = ctime(&now);
    string timeAndDate(dt);
    return timeAndDate;
}

string Logfile::currentTime() {
    //get current time
    time_t now = time(0);
    tm *ltm = localtime(&now);
    string time=to_string(ltm->tm_hour)+":"+to_string(ltm->tm_min)+":"+to_string(ltm->tm_sec);
    return time;
}

int Logfile::timeElapsedSeconds() {
    //get time elapsed since last call in seconds
    time1=chrono::high_resolution_clock::now();
    int deltaTime=chrono::duration_cast<chrono::seconds>(time1-time0).count();
    time0=time1;
    return deltaTime;
}

int Logfile::timeElapsedMinutes() {
    //get time elapsed since last call in seconds
    time1=chrono::high_resolution_clock::now();
    int deltaTime=chrono::duration_cast<chrono::seconds>(time1-time0).count();
    deltaTime/=60;
    time0=time1;
    return deltaTime;
}

void Logfile::initialise(){
    //create file with date and time header
    ofstream file(filename,ios::in|ios::trunc);
    string timeAndDate=currentTimeAndDate();
    string message="Simulation run: "+timeAndDate;
    writeFileValue(file,message);
    time0=chrono::high_resolution_clock::now();
    currentTime();
    file.close();
}

void Logfile::errorlog(string message, string errorType, int errorCode){
    //write error log to file and exit if severe
    ofstream file(filename,ios::in|ios::app);
    string errorMessage;
    if(errorType=="critical") errorMessage="Critical error at "+currentTime()+": "+message;
    writeFileValue(file,errorMessage);
    if(errorType=="critical") exit(errorCode);
    writeFileValue(file,message);
    file.close();
}
