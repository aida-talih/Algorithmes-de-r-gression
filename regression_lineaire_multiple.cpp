#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <iostream>
#include <math.h>
#define Nmax 2
#define Mmax 32
using namespace std;


class LinearRegression {

public:

    double x[Mmax][2];

    double *y;

    int m;

    double *theta;

    LinearRegression(double x[][Nmax], double y[], int m);

    void train(double alpha);

    double predict(double x[][Nmax], int index);

private:

    static double compute_cost(double x[][Nmax], double y[], double theta[], int m, int index);

    static double h(double x[][Nmax], double theta[],int index,int m);

    static double *calculate_predictions(double x[][Nmax], double theta[], int m );

    static double *gradient_descent(double x[][Nmax], double y[], double alpha,int m);

};

LinearRegression::LinearRegression(double x[][Nmax], double y[], int m) {
    for(int i=0;i<Mmax;i++){
        for(int j=0;j<2;j++){
            this->x[i][j]=x[i][j];

        }
    }

    this->y = y;
    this->m = m;
}

double array_sum(double arr[][Nmax], int len,int index) {
    double s = 0;

    for (int i = 0; i < Mmax; ++i) {
        s += arr[i][index];
    }

    return s;
}



double* array_multiplication(double arr1[][Nmax], double arr2[], int len,int index) {
    double *arr = new double[Mmax];

    for (int i = 0; i < Mmax; ++i) {
        arr[i] = arr1[i][index] * arr2[i];

    }

    return arr;
}

double* array_diff(double arr1[], double arr2[], int len) {
    double *arr = new double[Mmax];

    for (int i = 0; i < Mmax; ++i) {
        arr[i] = arr1[i] - arr2[i];
    }

    return arr;
}




double LinearRegression::h(double x[][Nmax], double theta[],int index,int m) {
    double temp=0;
    for(int i=0;i<Nmax;i++){
        temp=temp+theta[i]*x[index][i];
    }
    return temp;
}


double *LinearRegression::calculate_predictions(double x[][Nmax], double theta[], int m) {
    double *predictions = new double[Mmax];

    // calculate h for each training example

   for (int i = 0; i < Mmax; ++i) {
       predictions[i] = h(x, theta,i,Nmax);
   }

   return predictions;
}

void LinearRegression::train(double alpha) {

    this->theta = gradient_descent(x, y, alpha,m);
    cout << "Theta: " ;
    for(int i=0;i<Nmax;i++){
     cout << theta[i];

    }
}



double *LinearRegression::gradient_descent(double x[][Nmax], double y[], double alpha,int m) {

    double *theta = new double[Nmax]; double d=0; double *somme; double temp; double errors_x2[Mmax][2];

    for(int j=0;j<Nmax;j++){
        theta[j]=1;
    }

    do{
        for(int p=0;p<Nmax;p++){
            somme[p]=0;
        }
        double *predictions = calculate_predictions(x, theta, Mmax);

        double *diff = array_diff(predictions, y, Mmax);

        for(int i=0;i<Nmax;i++){
            double *tp=array_multiplication( x,diff, Mmax,i);

            for(int j=0;j<Mmax;j++){
                     errors_x2[j][i]=tp[j];

                   }
        }

        for(int j=0;j<Nmax;j++){
            somme[j]=(1.0 / Mmax) * array_sum(errors_x2, m,j);
            theta[j] = theta[j] - alpha * somme[j];

        }
        temp=0;
        for(int i=0;i<Nmax;i++){
            temp=temp+somme[i]*somme[i];
        }

 } while(sqrt(temp)>0.01);  //norme

    return theta;
}

int main()
{


double x[32][2]={
{1,4},
{1,4},
{1,7},
{1,7},
{1,8}
,{1,9}
,{1,10}
,{1,10}
,{1,10}
,{1,11}
,{1,11}
,{1,12}
,{1,12}
,{1,12}
,{1,12}
,{1,13}
,{1,13}
,{1,13}
,{1,13}
,{1,14}
,{1,14}
,{1,14}
,{1,14}
,{1,15}
,{1,15}
,{1,15}
,{1,16}
,{1,16}
,{1,17}
,{1,17}
,{1,17}
,{1,18}};




double slope[]={
2,
10
,40
,22
,16
,10
,18
,26
,34
,17
,28
,14
,20
,24
,28
,26
,34
,34
,46
,26
,36
,60
,80
,20
,26
,54
,32
,40
,32
,40
,50
,42
};

   LinearRegression lr(x, slope, 2);
  cout << "Enter learning rate alpha (default: 0.01): ";
    double alpha=0.0001;

    cout << "Training model..." << endl;
    lr.train(alpha);


    return 0;
}

