#ifndef CLOUD_CIPHER_H
#define CLOUD_CIPHER_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <math.h>
#include <NTL/RR.h>
#include <util.h>
#include <param.h>
#include <enc.h>
#include <cloud_util.h>
#include <NTL/ZZX.h>
#include <util.h>
#include <helib/helib.h>

using namespace std;
using namespace NTL;
class Cloud_cipher{
    public:
    Cloud_cipher();
    void weight_plainGen(int kerN,int kerD,int kerH, int kerW, string orignalweight,string orignalbias,string plainfile);
    void dense_weight_plainGen(int den_num1,int fc1_w,string w_dense1,string b_dense1,string fc_weightpliantex);
    void homMultHaveAdd(string outputcipherfile,string weightplainfile, const char* inputfile,int kerN,int kerD,int kerH,int kerW,int inH,int inW,int strid_H,int strid_W,string padding); 
    void homMultHaveAdd(string outputfile,string weightfile, string inputfile,int dennum,int fc1_w);
    void homMultConv_Conv1(int packedRowCount, int packedColCount, int packedRowUnit, int packedColUnit, string outputcipherfile, string weightplainfile, const char *inputfile, int kerN, int kerD, int kerH, int kerW, int inH, int inW, int strid_H, int strid_W, string padding);
    void homMultConv_Conv2(string outputcipherfile, string weightplainfile, const char *inputfile, int kerN, int kerD, int kerH, int kerW, int inH, int inW, int strid_H, int strid_W, int pool_strid_W, int pool_strid_H, string padding, ZZ_pX pm,  ZZ_pX p_temp1,  ZZ_pX p_temp2,  ZZ_pX p_temp3,  ZZ_pX p_temp4);
    void cipherAdd_Noise(int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W, int pool_stridH, int pool_stridW, string outcipherfile1, string cipherfile1,string cipherfile2);
    void input_cipherGen_Noslide(int***C,int inH,int inW,int kerN,int kerH,int kerW, int packedNum_row, int packedNum_col,const char *cipherfile,string padding);  
    void divNoise(ZZX *d, int kerN, string outputcipherfile, string inputcipher, helib:: EncryptedArray ea);
    void divNoise(ZZX *d,  string outputcipherfile, string inputcipher, helib:: EncryptedArray ea);
    void input_cipherGen_Noslide_twoPart(int***C_1,int***C_2,int inH,int inW,int kerN,int kerH,int kerW, int packedNum_row, int packedNum_col,const char *cipherfile,string padding);
    void AbsNoise_MulNoiseFC(ZZX *d, ZZX **u1, ZZX *v1, int outputlen,string outcipherfile,string cipherfile,helib:: EncryptedArray ea);
void AbsNoise_MulNoise(int**r, ZZX *d, int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W, int pool_stridH, int pool_stridW, string outcipherfile,string cipherfile,helib:: EncryptedArray ea,ZZX *t, ZZX **A, ZZ_pX pm, ZZ_pX p_temp1, ZZ_pX p_temp2, ZZ_pX p_temp3, ZZ_pX p_temp4, ZZ q2d, enc_param ***usp);
void plainEnc_Add(ZZX**R, int**C,int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W, int pool_stridH, int pool_stridW, const char *outcipherfile1, const char *cipherfile1,string padding,ZZX *t, ZZX **A,ZZ_pX pm, ZZ_pX p_temp1, ZZ_pX p_temp2, ZZ_pX p_temp3, ZZ_pX p_temp4, ZZ q2d, enc_param ***usp);

// void plainIPPacking_Enc(ZZX**R, int**C,int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W, int pool_stridH, int pool_stridW, const char *cipherfile1,string padding,ZZX *t, ZZX **A,ZZ_pX pm, ZZ_pX p_temp1, ZZ_pX p_temp2, ZZ_pX p_temp3, ZZ_pX p_temp4, ZZ q2d, enc_param ***usp);
    // void dense_input_cipherGen(int* C, int lenth_C,string fc_inputcipherfile);

};
#endif // CLIENT_UTIL_H