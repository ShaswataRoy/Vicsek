#include<stdio.h>
#include<math.h>
#include <random>
#include <iostream>

using namespace std;

double x[10],y[10];

double dist(int i,int j)
{
    return(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])));
}

double uniform(double a,double b)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(a, b);
    return(dis(gen));
}


void initialize_array(double *array,double a,double b)
{
    for(int i=0;i<10;i++)
    {
        array[i] = uniform(a,b);
    }
}

int main()
{
    initialize_array(x,0,10);
    initialize_array(y,0,10);
    if(dist(1,2)<(10/sqrt(3)))
    {
        cout<<"Yes"<<"\n";
    }
    else
    {
        cout<<"No"<<"\n";
    }
}
