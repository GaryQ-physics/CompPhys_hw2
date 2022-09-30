#include <stdio.h>
// #include <stdlib.h>
#include <math.h>
#include <sys/random.h>

#define MMAX 25

double rand(){
    int r;
    unsigned int max=4294967295;
    unsigned int tmp;
    r=getrandom(&tmp, sizeof(unsigned int), GRND_NONBLOCK);
    if(r!=4){
        printf("random returned something other than 4: %d\n",r);
    }
    return (double)tmp/max;
}

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

double simp_intg(double (*integrand_ptr)(double), double start, double end, long N){
    // only works for even N
    double I=0;
    double dx = (end-start)/N;

    I = (integrand_ptr(start)+integrand_ptr(end));
    for(long i=1; i<N; i++){
        I += ( ((i%2)+1)*2. )*integrand_ptr(start+i*dx);
    }
    I=I*dx/3;
    return I;
}

double mid_ballint_3D(long N){
    // 2D integration for 3D ball
    double x,y,d,I;

    // end and start are 1 and -1 for each
    d = 2./N; // dx and dy

    I=0.;
    for(long i=0; i<N; i++){
        for(long j=0; j<N; j++){
            x = -1.+i*d+0.5*d;
            y = -1.+j*d+0.5*d;
            if(x*x + y*y < 1.){
                I += sqrt(1. - x*x - y*y);
            }
        }
    }
    I=I*d*d;
    return I;
}

double mc_ballint_3D(long N){
    // 2D integration for 3D ball
    double I,x,y;

    I=0.;
    for(long i=0; i<N*N; i++){
        x = 2*rand() - 1.; // random between -1 and 1
        y = 2*rand() - 1.;
        if(x*x + y*y < 1.){
            I += sqrt(1. - x*x - y*y);
        }
    }
    I=4.*I/(N*N);
    return I;
}

void doball(){
    double I,I_exact;
    long newMAX = 13;
    double err[newMAX];
    // n=3 ball
    I_exact = 2*3.14159265358979/3.;
    printf("I_exact %f\n", I_exact);

    for(int M=1; M<newMAX+1; M++){
        I = mid_ballint_3D(MYexp2(M));
        err[M-1] = fabs(I-I_exact)/I_exact;
    }
    doubleArrToFile(err, newMAX, "data/3-ball-mid-error.data");

    for(int M=1; M<newMAX+1; M++){
        I = mc_ballint_3D(MYexp2(M));
        err[M-1] = fabs(I-I_exact)/I_exact;
    }
    doubleArrToFile(err, newMAX, "data/3-ball-mc-error.data");
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
    for(int M=1; M<MMAX+1; M++){
        I = simp_intg(&mon6,0.,1.,MYexp2(M));
        err[M-1] = fabs(I-I_exact)/I_exact;
    }
    doubleArrToFile(err, MMAX, "data/mon6-simp-error.data");

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
    for(int M=1; M<MMAX+1; M++){
        I = simp_intg(&circ,0.,1.,MYexp2(M));
        err[M-1] = fabs(I-I_exact)/I_exact;
    }
    doubleArrToFile(err, MMAX, "data/circ-simp-error.data");

    doball();
    return 0;
}
