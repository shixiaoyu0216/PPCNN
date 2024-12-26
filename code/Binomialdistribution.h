#ifndef BINOMIALDISTRIBUTION_H
#define BINOMIALDISTRIBUTION_H
#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/RR.h>
// #include <util.h>
#include <ctime>
#include <sys/timeb.h>
using namespace std;
using namespace NTL;

class Binomialdistribution
{
public:
    Binomialdistribution();
    void Binomial_distribution_gen(ZZX& s, int n,int ita);
};

#endif // BINOMIALDISTRIBUTION_H
