#include<stdio.h>
#include <cstring>
#include<math.h>
#include <sys/time.h>
#include <random>
#include <iostream>
#include <omp.h>

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

double array_std(double *array,int min,int max)
{
    double sum=0.0;
    double mean=array_mean(array,min,max);
    for(int i=min;i<max;i++)
    {
        sum+=(array[i]-mean)*(array[i]-mean);
    }
    return(sqrt(sum/(max-min)));
}


double array_binder(double *array,int min,int max)
{
    double second = 0.0;
    double fourth = 0.0;
    for(int i=min;i<max;i++)
    {
        second += pow(array[i],2);
        fourth += pow(array[i],4);
    }
    fourth = fourth/(max-min);
    second = second/(max-min);
    return(1-fourth/(second*second*3));
}

//Plot

double compute_order(int no,double eta)
{
    initialize();
    int a,b;
    double eps=0.;
    int repeat =1;

    a=200;b=10;eps=0.01;repeat=1000;


    double order_arr;
    double theta_cos,theta_sin;
    double ordall[MAX_TIME];
    int t=0;


    for(int r=0;r<repeat;r++)
    {
        initialize();

        for(t=0;t<MAX_TIME;t++)
        {
            if(t>a)
            {
                if(fabs(array_std(ordall,t-2*b,t-b)-
                        array_std(ordall,t-b,t))<eps)
                    break;
            }
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

        order_arr += array_std(ordall,t*3/4,t);
    }

    return(order_arr/repeat);

}

double logistic(double x)
{
    return(x*(1-x));


}

void golden(double (*func)(double x),double a,double c,double b,double fa,double fc,double fb,int i)
{
    double phi = (1+sqrt(5))/2;

    if(i==1)
    {
        double d = b - (b-a)/phi;
        double fd = (*func)(d);
        if(fabs(fc-fd)<0.1)
            return;
        fprintf(fp,"%lf\t%lf",d,fd);
        if(fd<fc)
            golden(func,a,c,d,fa,fc,fd,2);
        else
            golden(func,c,d,b,fc,fd,fb,1);
    }

    else if(i==2)
    {
        double d = a + (b-a)/phi;
        double fd = (*func)(d);
        if(fabs(fc-fd)<0.1)
            return;
        fprintf(fp,"%ld\t%ld",d,fd);
        if(fd>fc)
            golden(func,a,d,c,fa,fd,fc,2);

        else
            golden(func,d,c,b,fd,fc,fb,1);
    }
}

int main()
{

    struct timeval start, end;
    gettimeofday(&start, NULL);
    N=60;L=3.873;

    // benchmark code
    string str = "vicsek";
    str+=to_string(N)+"s.txt";

    //fp = fopen(str.c_str(),"w");
    fp = fopen("logistic.txt","w");
    golden(logistic,0,(-1+sqrt(5))/2,1,0,logistic((-1+sqrt(5))/2),0,1);

    //compute_order(N);
    fclose(fp);
    gettimeofday(&end, NULL);

    double delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
             end.tv_usec - start.tv_usec) / 1.e6;
    cout<<"Elasped time is "<< delta<<" seconds.\n";
    return(0);
}
