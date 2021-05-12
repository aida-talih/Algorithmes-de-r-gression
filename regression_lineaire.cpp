#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <iostream>
#include <math.h>

using namespace std;


class LinearRegression {

public:

    // First feature
    double *x;

    // Target feature
    double *y;

    // Number of training examples
    int m;

    // The theta coefficients
    double *theta;


    LinearRegression(double x[], double y[], int m);


    void train(double alpha);


    double predict(double x);

private:


    static double h(double x, double theta[]);


    static double *calculate_predictions(double x[], double theta[], int m);


    static double *gradient_descent(double x[], double y[], double alpha, int m);

};

LinearRegression::LinearRegression(double x[], double y[], int m) {
    this->x = x;
    this->y = y;
    this->m = m;
}


double array_sum(double arr[], int len) {
    double s = 0;

    for (int i = 0; i < len; ++i) {
        s += arr[i];
    }

    return s;
}



double* array_multiplication(double arr1[], double arr2[], int len) {
    double *arr = new double[len];

    for (int i = 0; i < len; ++i) {
        arr[i] = arr1[i] * arr2[i];
    }

    return arr;
}

double* array_diff(double arr1[], double arr2[], int len) {
    double *arr = new double[len];

    for (int i = 0; i < len; ++i) {
        arr[i] = arr1[i] - arr2[i];
    }

    return arr;
}


void LinearRegression::train(double alpha) {

    this->theta = gradient_descent(x, y, alpha, m);

    cout << endl << "Theta: " << theta[0] << " " << theta[1] << endl;
}

double LinearRegression::predict(double x) {
    return h(x, theta);
}



double LinearRegression::h(double x, double theta[]) {
    return theta[0] + theta[1] * x;
}

double *LinearRegression::calculate_predictions(double x[], double theta[], int m) {
    double *predictions = new double[m];


    for (int i = 0; i < m; ++i) {
        predictions[i] = h(x[i], theta);
    }

    return predictions;
}

double Norme(double a,double b) {

    return sqrt( pow(b,2)+pow(a,2));
}

double *LinearRegression::gradient_descent(double x[], double y[], double alpha, int m) {
    double *theta = new double[2];
    theta[0] = 1;
    theta[1] = 1;
    int i=0;
    double d=0;

 do{

        double *predictions = calculate_predictions(x, theta, m);
        double *diff = array_diff(predictions, y, m);

        double *errors_x1 = diff;
        double *errors_x2 = array_multiplication(diff, x, m);



        theta[0] = theta[0] - alpha * (1.0 / m) * array_sum(errors_x1, m);
        theta[1] = theta[1] - alpha * (1.0 / m) * array_sum(errors_x2, m);

        d=sqrt( ((1.0 / m)*array_sum(errors_x1, m)*(1.0 / m)*array_sum(errors_x1, m))+((1.0 / m)*array_sum(errors_x2, m)*(1.0 / m)*array_sum(errors_x2, m)));

       cout<<"la norme"<<d<<endl;


} while(d>0.001);

    return theta;
}

int main()
{

double diameter[]={
4,
4
,7
,7
,8
,9
,10
,10
,10
,11
,11
,12
,12
,12
,12
,13
,13
,13
,13
,14
,14
,14
,14
,15
,15
,15
,16
,16
,17
,17
,17
,18

};


  double alpha=0.0001;


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

    LinearRegression lr(diameter, slope, 32);


    cout << "Training model..." << endl;
    lr.train(alpha);



    return 0;
}

