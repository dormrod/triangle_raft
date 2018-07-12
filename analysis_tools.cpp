#include "analysis_tools.h"

DiscreteDistribution::DiscreteDistribution() {
    //default constructor
}

DiscreteDistribution::DiscreteDistribution(vector<int> values) {
    //calculate distribution from provided values

    if(values.size()==0){
        n=0;
        return;
    }

    //generate unique list of values
    sort(values.begin(),values.end());
    vector<int> uniqueValues;
    uniqueValues.clear();
    uniqueValues.push_back(values[0]);
    for (int i=1;i<values.size();++i)if (values[i]!=uniqueValues.rbegin()[0]) uniqueValues.push_back(values[i]);

    //set distribution values and size
    n=uniqueValues.size();
    x=col_vector<int>(uniqueValues);
    p_raw=col_vector<double>(n);
    p=col_vector<double>(n);

    //calculate probabilities
    map<int,int> xRef;
    for(int i=0; i<n; ++i) xRef[x[i]]=i;
    int ref;
    for(int i=0; i<values.size(); ++i){
        ref=xRef[values[i]];
        p_raw[ref]+=1.0;
    }
    for(int i=0; i<n; ++i) p=p_raw/values.size();

    //calculate mean
    mean=0.0;
    for(int i=0; i<n; ++i) mean+=x[i]*p[i];
}

vector<int> DiscreteDistribution::getValues() {
    vector<int> values;
    for(int i=0; i<n; ++i) values.push_back(x[i]);
    return values;
}

vector<double> DiscreteDistribution::getProbabilities() {
    vector<double> probabilities;
    for(int i=0; i<n; ++i) probabilities.push_back(p[i]);
    return probabilities;
}

double DiscreteDistribution::getProbability(int xValue) {
    double probability=0.0;
    for(int i=0; i<n; ++i){
        if(xValue==x[i]){
            probability=p[i];
            break;
        }
    }
    return probability;
}

vector<double> leastSquaresLinearRegression(vector<double> x, vector<double> y){
    //perform least squares linear regression and return gradient, intercept and r-squared

    vector<double> coefficients(3); //gradient, intercept and r-squared
    int n=x.size();
    double sumX=0.0, sumY=0.0, sumXY=0.0, sumXX=0.0, sumYY=0.0;
    for(int i=0; i<n; ++i){
        sumX=sumX+x[i];
        sumY=sumY+y[i];
        sumXY=sumXY+x[i]*y[i];
        sumXX=sumXX+x[i]*x[i];
        sumYY=sumYY+y[i]*y[i];
    }
    double sqSumX=sumX*sumX, sqSumY=sumY*sumY;
    double sdX=sqrt(n*sumXX-sqSumX), sdY=sqrt(n*sumYY-sqSumY);
    double r=(n*sumXY-sumX*sumY)/(sdX*sdY);
    coefficients[0]=(r*sdY)/sdX;
    coefficients[1]=(sumY-coefficients[0]*sumX)/n;
    coefficients[2]=r*r;
    return coefficients;
};
