#ifndef PARAMGEN_H
#define PARAMGEN_H
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

class KeyparaGen
{
public:
    KeyparaGen();
    void keygen();
};    

#endif //// KEYPARAGEN_H
