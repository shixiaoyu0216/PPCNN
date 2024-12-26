#include "Binomialdistribution.h"

Binomialdistribution::Binomialdistribution()
{

}

void Binomialdistribution::Binomial_distribution_gen(ZZX& s, int n, int ita){
    struct timeb timeSeed;
    ftime(&timeSeed);
    srand(timeSeed.time*1000+timeSeed.millitm);
    for (int i = 0; i < n; i++) {
        int *a = new int[n]();
        for (int j = 0; j < n; j++) {
            a[j] = rand() % 2;
        }
        int *b = new int[n]();
        for (int j = 0; j < n; j++) {
            b[j] = rand() % 2;
        }
        int s_temp = 0;
        for (int j = 0; j < ita; j++) {
            s_temp += a[j] - b[j];
        }
        SetCoeff(s,i,s_temp);
        delete a;
        delete b;
    }
}
