#include <stdio.h>
#include <stdlib.h>

void MultMatrix(double* , double* , double* , int, int, int);
//Первая матрица, Вторая матрица, Пустая матрица нужного размера, Количество строк первой, Количество столбцов первой, Количество столбцов второй
void lineswap(double* , int, int, int);
//Матрица, Длина строки, Номер первой строки, Номер второй строки
void lineaddition(double* , int, int, double, int);
//Матрица, Длина строки, Номер изменяемой строки, Коэффициент, Номер прибавляемой строки
void linemultp(double* , int, int, double);
//Матрица, Длина строки, Номер строки, Коэффициент
void columnswap(double* , int, int, int, int);
//Матрица, Длина столбца, Длина строки, Номер первого столбца, Номер второго столбца
void columnaddition (double* , int, int, int, double, int);
//Матрица, Длина столбца, Длина строки, Номер изменяемого столбца, Коэффициент, Номер прибавляемого столбца
void columnmultp(double* , int , int , int , double );
//Матрица, Длина столбца, Длина строки, Номер столбца, Коэффициент
long long determinator(double* , int);
//Матрица(квадратная), Размер
void obr(double* , double* , int);
//Матрица(квадратная), Пустая матрица того же размера, Размер
void trmatrix(double* , double* , int, int);
//Матрица, Пустая матрица того же размера, Количество строк, Количество столбцов
int rank(double* , int, int);
//Матрица, Количество строк, Количество столбцов
int main(){
    
}

void MultMatrix(double* d1, double *d2, double* c, int N, int M, int K){
    int i, j, l;

    for (i = 0; i < N; i++){
        for (l = 0; l < K; l++){
            c[i*K + l] = 0;
            for (j = 0; j < M; j++){
                c[i*K + l] += d1[i*M + j] * d2[j*K + l];
                }
        }

    }
}

void lineswap(double* a, int m, int k1, int k2){
    double c[m];
    int j;
    k1--;
    k2--;
    for (j = 0; j < m; j++){
            c[j] = a[k1*m + j];
            a[k1*m + j] = a[k2*m + j];
            a[k2*m + j] = c[j];
    }
}

void lineaddition(double* a, int m, int k1, double p, int k2){
    int j;
    k1--;
    k2--;
    for(j = 0; j < m; j++){
        a[k1*m + j] += p*a[k2*m + j]; 
    }
}

void linemultp(double* a, int m, int k, double p){
    int j;
    k--;
    for (j = 0; j < m; j++){
        a[k*m + j] *= p;
    }
}

void columnswap(double* a, int n, int m, int w1, int w2){
    double c[m];
    int j;
    w1--;
    w2--;
    for (j = 0; j < n; j++){
        c[j] = a[w1 + j*m];
        a[w1 + j*m] = a[w2 + j*m];
        a[w2 + j*m] = c[j];
    }
}

void columnaddition(double* a, int n, int m, int w1, double p, int w2){
    int j;
    w1--;
    w2--;
    for(j = 0; j < n; j++){
        a[w1 + j*m] += p*a[w2 + j*m]; 
    }
}

void columnmultp(double* a, int n, int m, int k, double p){
    int j;
    k--;
    for (j = 0; j < n; j++){
        a[k + j*m] *= p;
    }
}

long long determinator(double* z, int n){
    int i, j, k;
    int stepen = 0;
    int p = 1;
    int s = -1;
    long long det = 1;
    double a[n*n];
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            a[i*n + j] = z[i*n + j];
        }
    }
    for(j = 0; j < n; j++){
        for (i = j + 1; i < n;i++){
                if(a[(j)*n + j] != 0){
                    lineaddition(a, n, i + 1, -(a[(i)*n + j])/(a[(j)*n + j]), j + 1);
                } else {
                    for(k = j + 1; k < n; k++){
                        if(a[k*n + j] != 0){
                            lineswap(a, n, j + 1, k + 1);
                            stepen += 2*(k - j) -1;
                            break; 
                        }
                    }
                    lineaddition(a, n, i + 1, -(a[(i)*n + j])/(a[(j)*n + j]), j + 1);
                }
            }
    }
    for(i = 1; i <= stepen; i++) p = p*s;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            if (i == j){
                det *= a[i*n + j];
            }
        }
    }
    det *= p;
    return det;
}

void obr(double* z, double* a, int n){
    int i, j, k;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            a[i*n + j] = z[i*n + j];
        }
    }
    double ao[n*n];
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i==j){
                ao[i*n + j] = 1;
            } else ao[i*n + j] = 0;
        }
    }
    for(j = 0; j < n; j++){
        for (i = 0; i < n;i++){
                if(a[(j)*n + j] != 0){
                    linemultp(ao, n, j + 1, 1/a[j*n + j]);
                    linemultp(a, n, j + 1, 1/a[j*n + j]);
                    if(i != j){
                        lineaddition(ao, n, i + 1, -(a[(i)*n + j]), j + 1);
                        lineaddition(a, n, i + 1, -(a[(i)*n + j]), j + 1);
                    }
                } else {
                    for(k = j + 1; k < n; k++){
                        if(a[k*n + j] != 0){
                            lineswap(ao, n, j + 1, k + 1);
                            lineswap(a, n, j + 1, k + 1);
                            break; 
                        }
                    }
                    linemultp(ao, n, j + 1, 1/a[j*n + j]);
                    linemultp(a, n, j + 1, 1/a[j*n + j]);
                    if(i == j){
                        lineaddition(ao, n, i + 1, -(a[(i)*n + j]), j + 1);
                        lineaddition(a, n, i + 1, -(a[(i)*n + j]), j + 1);
                    }
                }
            }
    }
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            a[i*n + j] = ao[i*n + j];
        }
    }
}

void trmatrix(double* c, double* tc, int n, int m){
    int i, j;
    for (i = 0; i < n; i++){
        for (j = 0; j < m; j++){
            tc[i + j*n] = c[i*m + j];
        }
    }   
}

int rank(double* z, int n, int m){
    int i, j, k;
    int flag = 0;
    int r = n;
    double a[n*m];
    for (i = 0; i < n; i++){
        for (j = 0; j < m; j++){
            a[i*m + j] = z[i*m + j];
        }
    }
    for(j = 0; j < m; j++){
        for (i = j + 1; i < n;i++){
                if(a[(j)*m + j] != 0){
                    lineaddition(a, m, i + 1, -(a[(i)*m + j])/(a[(j)*m + j]), j + 1);
                    flag = 0;
                } else {
                    for(k = j + 1; k < m; k++){
                        if(a[k*m + j] != 0){
                            lineswap(a, m, j + 1, k + 1);
                            flag = 1;
                            break; 
                        }
                    }
                    if (flag != 1){
                        flag = 0;
                        break;
                    }
                    lineaddition(a, m, i + 1, -(a[(i)*m + j])/(a[(j)*m + j]), j + 1);
                    flag = 0;
                }
            }
    }
    for (i = 0; i < n; i++){
        for (j = 0; j < m; j++){
            if(a[i*m + j] != 0) flag = 1;
        }
        if (flag == 0) r--;
        flag = 0;
    }
    return r;
}