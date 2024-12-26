#ifndef ENC_H
#define ENC_H
#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/RR.h>
#include <Binomialdistribution.h>
#include <util.h>
#include <ctime>
using namespace std;
using namespace NTL;

class Enc
{
public:
    Enc();
    void encrypt(ZZX m,ZZX u[],ZZX& v,int n,ZZ q,int k,ZZ q2d,ZZX r[],ZZX **A,ZZX e1[],ZZX e2,ZZ_pX p_temp1,ZZ_pX p_temp2,ZZ_pX p_temp3,ZZ_pX p_temp4,ZZ_pX pm,ZZX t[],int dt,int du ,int dv);
    void decrypt(ZZX& m,ZZX u[],ZZX& v,int n,ZZ q,int k,ZZX s[], ZZ_pX p_temp1,ZZ_pX p_temp2,ZZ_pX p_temp3,ZZ_pX p_temp4,ZZ_pX pm,int du ,int dv,int dp);
};

#endif // ENC_H
