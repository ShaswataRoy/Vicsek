#include<stdio.h>
#include <cstring>
#include<math.h>
#include <sys/time.h>
#include <random>
#include <iostream>
#include <omp.h>
#include "progress.h"

using namespace std;

//Define Variables
const int MAX_SIZE=10000;
const int MAX_TIME=10000;
double v=0.03;
double eta=0.1;
double L=31;
int N=4000;

//double ordall[MAX_TIME];
FILE *fp;


//Initial State
double x[MAX_SIZE],y[MAX_SIZE],
    theta[MAX_SIZE],theta_new[MAX_SIZE];
double order,order_new;

//Initialize
double uniform(double a,double b)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(a, b);
    return(dis(gen));
}

void initialize_array(double *array,double a,double b)
{
    for(int i=0;i<N;i++)
    {
        array[i] = uniform(a,b);
    }
}

void initialize()
{
    initialize_array(x,0,L);
    initialize_array(y,0,L);
    initialize_array(theta,-M_PI,M_PI);
    initialize_array(theta_new,-M_PI,M_PI);

    order=0;
    order_new=1;
}

//Functions

double dist(int i,int j)
{
    return(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])));
}

void update()
{
    int k;
    double theta_near;
    for(int i=0;i<N;i++)
    {
        theta_near = 0;
        k=0;

        for(int j=0;j<N;j++)
        {
            if(dist(i,j)<1.0)
            {
                theta_near+=theta[j];
                k+=1;
            }
        }
        if(k!=0)
        {
            //theta_new[i] = theta_near/k + eta*uniform(-M_PI,M_PI);
            theta_new[i] = theta_near/k + uniform(-eta/2,eta/2);
        }
    }

    for(int i=0;i<N;i++)
    {
        theta[i]=theta_new[i];
        x[i]=fmod((x[i]+v*cos(theta[i])),L);
        y[i]=fmod((y[i]+v*sin(theta[i])),L);
    }
}

double array_mean(double *array,int min,int max)
{
    double sum=0.0;
    for(int i=min;i<max;i++)
    {
        sum+=array[i];
    }
    return(sum/(max-min));
}
//Plot

void compute_order(int no)
{
    initialize();
    int a,b;
    double eps=0.;
    int repeat =1;
    const int no_points=20;
    const double eta_max=2*M_PI;


    double order_arr[no_points],eta_arr[no_points];
    double theta_cos,theta_sin;
    double ordall[MAX_TIME];
    int t=0;
    for(int i=0;i<no_points;i++)
    {
        eta_arr[i] = eta_max*i/(1.0*no_points);
    }

    initialize();

    for(int i=0;i<no_points;i++)
    {
        eta=eta_arr[i];
        for(t=0;t<50;t++)
        {

            theta_cos = 0.;
            theta_sin = 0.;
            update();

            for(int j=0;j<no;j++)
            {
                theta_cos+=cos(theta[j]);
                theta_sin+=sin(theta[j]);
            }
            order_new=sqrt(theta_cos*theta_cos+theta_sin*theta_sin)/no;
            ordall[t]=order_new;
        }

        order_arr[i] += array_mean(ordall,t*3/4,t);
    }





}

int main(int argc,char* argv[])
{

    struct timeval start, end;
    gettimeofday(&start, NULL);
    cout<<"\nN = "<<argv[1]<<"\n";
    N=atoi(argv[1]);L=atof(argv[2]);

    // benchmark code


    compute_order(N);

    gettimeofday(&end, NULL);

    double delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
             end.tv_usec - start.tv_usec) / 1.e6;
    cout<<"\rElasped time is "<< delta<<" seconds.\n";
    return(0);
}
