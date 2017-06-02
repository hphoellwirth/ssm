// ----------------------------------------------------------------------
// Information
// ----------------------------------------------------------------------

// Continuous Sequential Importance Sampling Algorithm
//
// (Author) Hans-Peter Höllwirth
// (Date)   06.2017

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double *array;

// compare function for index sorting
int cmp(const void *a, const void *b){
    int ia = *(int *)a;
    int ib = *(int *)b;
    return array[ia] < array[ib] ? -1 : array[ia] > array[ib];
}
 
// Output:  s   - particle filtered (continuous version)
// Input 1: p   - predictive density
// Input 2: w   - normal pdf evaluated at y[t]
// Input 3: u   - sorted uniformly random sampled vector (rejection sampling)
// Input 4: len - length of vectors s, p, w, u
void csir(double *s, double *p, double *w, double *u, int *len) {
    
    int P = len[0];
    int i; 

    // standardize w
    double sum_w = 0.0;    
    for (i = 0; i < P; i++)
        sum_w += w[i];
        
    for (i = 0; i < P; i++) 
        w[i] = w[i] / sum_w;
    
    // sort p index
    int p_idx[P];
    for (i = 0; i < P; i++) {
        p_idx[i] = i;
    } 
    array = p;
    qsort(p_idx, P, sizeof(*p_idx), cmp);    
    
    // compute cumulated sum of w
    double w_cum[(P+1)];
    w_cum[0] = 0.0;
    for (i = 1; i <= P; i++) {
        w_cum[i] = w_cum[i-1] + w[p_idx[i-1]];       
    }  
    
    // duplicate first element of p 
    double * p_new = malloc(P * sizeof(p[0]));
    memcpy(p_new, p, (P+1) * sizeof(p[0]));
    p_new[0] = p[p_idx[0]];
    for (i = 0; i < P; i++)
        p_new[i+1] = p[p_idx[i]];    

    // compute s
    int j = 0;
    for (i = 0; i < P; i++) {
        while((w_cum[i] < u[j]) && (u[j] <= w_cum[i+1])){
            s[j] = p_new[i] + ((p_new[i+1]-p_new[i])/(w_cum[i+1]-w_cum[i])) * (u[j]-w_cum[i]);
            if (j < P) {
                j++;
            }
            else break;
        }
    }  
}

/*
// test function csir
int main(void){
    int P[] = {5};
    double p[] = {1.737899, 2.247636, 2.813836, 2.521015, 2.959267};
    double w[] = {0.02438808, 0.04078347, 0.05066587, 0.04691185, 0.05149138};
    double u[]        = {0.01631192, 0.04410280, 0.05565297, 0.62613028, 0.79035955};

    double s[P[0]];    
    //double *s;
    csir (s, p, w, u, P);
    
    printf("s\n");
    for (int i=0; i<P[0]; i++) { 
        printf("%lf\n", s[i]);
    }    
    
    return 0;
}
*/





