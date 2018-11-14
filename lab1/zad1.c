#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define size 10000000

float sum_iter(float *tab, int len){
    float sum = 0;
    for(int i=0; i<len; i++) sum += tab[i];
    return sum;
}

float sum_rec(float *tab, int l, int r){
    if(l > r) return 0;
    else if(l == r) return tab[l];
    else{
        int m = (l+r)/2;
        return sum_rec(tab, l, m) + sum_rec(tab, m+1,r);
    }
}

float sum_kahan(float* tab, int len){
    float sum = tab[0];
    float c = 0.0;             
    for(int i = 1; i<len; i++){
        float y = tab[i] - c;
        float t = sum + y;           
        c = (t - sum) - y;    
        sum = t;                   
    }                  
    return sum;
}

float my_abs(float v){
    if(v >= 0) return v;
    else return -v;
}

float calcError(float v1, float v2) {
    return my_abs(v1-v2);
}

float calcAbsError(float v1, float v2) {
    return my_abs(v1-v2)/v1;
}

int main(){
    float *tab = malloc(sizeof(float) * size);
    float sum1, sum2, sum3, err1, err2, expect, e;
    clock_t start, end;
    double cpu_time_used;
    e = 0.135791;
    for(int i=0; i<size; i++) tab[i] = e;
    int len = 1;
    for(int i=1; i<=7; i++){
        len *= 10;
	    expect = e*len;
        
        start = clock();
	    sum1 = sum_iter(tab, len);
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	    
        err2 = calcAbsError(expect, sum1);
        printf("%f %f\n", err2, cpu_time_used);	

        start = clock();
	    sum1 = sum_rec(tab, 0, len);
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

        err2 = calcAbsError(expect, sum1);
        printf("%f %f\n", err2, cpu_time_used);

        start = clock();
	    sum1 = sum_kahan(tab, len);
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

        err2 = calcAbsError(expect, sum1);
        printf("%f %f\n\n", err2, cpu_time_used);	
        	    
/*printf("SUM ITER\n%dValue: %f\nExpected: %f\nResult: %f\nError: %f\nAbsError: %f\n\n", len, e, expect, sum1, err1, err2);*/
    }
    return 0;
}
/*
    for(int i=1; i<8; i++){
    e = 0.135791 * i;
    expect = e*size;

    for(int i=0; i<size; i++) tab[i] = e;
    
    start = clock();
    sum1 = sum_iter(tab);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    err1 = calcError(expect, sum1);
    err2 = calcAbsError(expect, sum1);
    printf("SUM ITER\nError: %f\nTime: %f\n\n", err2, cpu_time_used);


    start = clock();
    sum2 = sum_rec(tab, 0, size);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    err1 = calcError(expect, sum2);
    err2 = calcAbsError(expect, sum2);
    printf("SUM REC\nError: %f\nTime: %f\n\n", err2, cpu_time_used);


    start = clock();
    sum3 = sum_kahan(tab);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    err1 = calcError(expect, sum3);
    err2 = calcAbsError(expect, sum3);
printf("SUM KAHAN\nError: %f\nTime: %f\n\n", err2, cpu_time_used);
    }*/

/*   
    start = clock();
    sum1 = sum_iter(tab);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    diff = calcError(sum1, expect); 
    error = diff/sum1;
    printf("SUM ITER\nExpected: %f\nResult: %f\nDiffrenece: %f\nError: %f\nTime: %f\n\n", expect, sum1, diff, error, cpu_time_used);
    
    start = clock();
    sum2 = sum_rec(tab, 0, size);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    diff = sum2 - expect; 
    error = diff/sum2;
    printf("SUM REC\nExpected: %f\nResult: %f\nDiffrenece: %f\nError: %f\nTime: %f\n\n", expect, sum2, diff, error, cpu_time_used);

    start = clock();
    sum3 = sum_kahan(tab);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    diff = sum3 - expect; 
    error = diff/sum3;
    printf("SUM KAHAN\nExpected: %f\nResult: %f\nDiffrenece: %f\nError: %f\nTime: %f\n\n", expect, sum3, diff, error, cpu_time_used);
*/
