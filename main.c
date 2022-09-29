#include <stdio.h>
//#include <stdlib.h>
#include <math.h>

#define MMAX 25

double mon6(double x){
    return x*x*x*x*x*x;
}
double circ(double x){
    return sqrt(1.-x*x);
}

long MYexp2(long m){
    long ret=1;
    for(long i=0; i<m; i++){
        ret = 2*ret;
    }
    return ret;
}

void doubleArrToFile(double* arr, long arrlen, char* filename){
    FILE * fp=fopen(filename,"w");
    fwrite(arr, sizeof(double), arrlen, fp);
    fclose(fp);
}

double L_intg(double (*integrand_ptr)(double), double start, double end, long N){
    double I=0.;
    double dx = (end-start)/N;
    for(long i=0; i<N; i++){
        I += integrand_ptr(start+i*dx);
    }
    I=I*dx;
    return I;
}

double R_intg(double (*integrand_ptr)(double), double start, double end, long N){
    double I=0.;
    double dx = (end-start)/N;
    for(long i=1; i<N+1; i++){
        I += integrand_ptr(start+i*dx);
    }
    I=I*dx;
    return I;
}

double trap_intg(double (*integrand_ptr)(double), double start, double end, long N){
    double I=(integrand_ptr(start)+integrand_ptr(end))/2;
    double dx = (end-start)/N;
    for(long i=1; i<N; i++){
        I += integrand_ptr(start+i*dx);
    }
    I=I*dx;
    return I;
}

double mid_intg(double (*integrand_ptr)(double), double start, double end, long N){
    double I=0.;
    double dx = (end-start)/N;
    for(long i=0; i<N; i++){
        I += integrand_ptr(start+i*dx+0.5*dx);
    }
    I=I*dx;
    return I;
}

int main(){
    double I,I_exact;
    double err[MMAX];

    // f(x)=x^6
    I_exact = 1./(6+1);

    for(int M=1; M<MMAX+1; M++){
        I = L_intg(&mon6,0.,1.,MYexp2(M));
        err[M-1] = fabs(I-I_exact)/I_exact;
    }
    doubleArrToFile(err, MMAX, "data/mon6-L-error.data");
    for(int M=1; M<MMAX+1; M++){
        I = R_intg(&mon6,0.,1.,MYexp2(M));
        err[M-1] = fabs(I-I_exact)/I_exact;
    }
    doubleArrToFile(err, MMAX, "data/mon6-R-error.data");
    for(int M=1; M<MMAX+1; M++){
        I = trap_intg(&mon6,0.,1.,MYexp2(M));
        err[M-1] = fabs(I-I_exact)/I_exact;
    }
    doubleArrToFile(err, MMAX, "data/mon6-trap-error.data");
    for(int M=1; M<MMAX+1; M++){
        I = mid_intg(&mon6,0.,1.,MYexp2(M));
        err[M-1] = fabs(I-I_exact)/I_exact;
    }
    doubleArrToFile(err, MMAX, "data/mon6-mid-error.data");

    // f(x)=sqrt(1-x^2)
    I_exact = 3.14159265358979/4.;
    printf("I_exact %f\n", I_exact);

    for(int M=1; M<MMAX+1; M++){
        I = L_intg(&circ,0.,1.,MYexp2(M));
        err[M-1] = fabs(I-I_exact)/I_exact;
    }
    doubleArrToFile(err, MMAX, "data/circ-L-error.data");
    for(int M=1; M<MMAX+1; M++){
        I = R_intg(&circ,0.,1.,MYexp2(M));
        err[M-1] = fabs(I-I_exact)/I_exact;
    }
    doubleArrToFile(err, MMAX, "data/circ-R-error.data");
    for(int M=1; M<MMAX+1; M++){
        I = trap_intg(&circ,0.,1.,MYexp2(M));
        err[M-1] = fabs(I-I_exact)/I_exact;
    }
    doubleArrToFile(err, MMAX, "data/circ-trap-error.data");
    for(int M=1; M<MMAX+1; M++){
        I = mid_intg(&circ,0.,1.,MYexp2(M));
        err[M-1] = fabs(I-I_exact)/I_exact;
    }
    doubleArrToFile(err, MMAX, "data/circ-mid-error.data");

    return 0;
}
