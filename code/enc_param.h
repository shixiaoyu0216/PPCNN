#ifndef ENC_PARAM_H
#define ENC_PARAM_H
#include <NTL/ZZX.h>
#include <Binomialdistribution.h>
using namespace std;
using namespace NTL;
struct enc_param{
    public:
        Binomialdistribution Bin_dis;
        ZZX r1[2];
        ZZX e11[2];
        ZZX e12;
};
#endif//ENC_PARAM_H