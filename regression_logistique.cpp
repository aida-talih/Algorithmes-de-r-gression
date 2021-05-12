#include<bits/stdc++.h>

using namespace std;
bool custom_sort(double a, double b)
{
    double  a1=abs(a-0);
    double  b1=abs(b-0);
    return a1<b1;
}


double norme(double x,double y,double z){

return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
}

int main()
{
    int i=0;
    int l=0,k;
    double score[4];

    double x1[] = {2.7810836, 1.465489372, 3.396561688, 1.38807019, 3.06407232,7.627531214,5.332441248,6.922596716,8.675418651 ,7.673756466};
    double x2[] = {2.550537003,2.362125076,4.400293529,1.850220317,3.005305973,2.759262235,2.088626775,1.77106367,-0.2420686549,3.508563011};
    double y[] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};


vector<double>error; // for storing the error values
double err;    // for calculating error on each stage
double b0 = 0; // initializing b0
double b1 = 0; // initializing b1
double b2=  0;
double d;// initializing b2
double alpha = 0.001; // initializing our learning rate
double  e = 2.71828;
double s0,s1,s2;
/*Training Phase*/
double p;
double pred;


do{
s0=0;s1=0;s2=0;

   for(int idx=0;idx<10;idx++){
    p = -(b0 + b1 * x1[idx]+ b2* x2[idx]);
    pred  = 1/(1+ pow(e,p));

    err = pred-y[idx];

     s0=s0+err* 1.0;

     s1=s1+err *x1[idx];

     s2=s2+err* x2[idx];

     error.push_back(err);
     }
     s0=s0*0.1;
     s1=s1*0.1;
     s2=s2*0.1;

      b0 = b0 - (alpha *s0);
      b1 = b1 - (alpha *s1);
      b2 = b2 -(alpha*s2);

    cout<<"Final Values are: "<<"B0="<<b0<<" "<<"B1="<<b1<<" "<<"B2="<<b2<<" error="<<error[0]<<endl;

    cout<<s0<<endl;cout<<s1<<endl;cout<<s2<<endl;

    d=norme(s0,s1,s2);

    cout<<d<<endl;

}while(d>0.001);


sort(error.begin(),error.end(),custom_sort);

return 1;
}
