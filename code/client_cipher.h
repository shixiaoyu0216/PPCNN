#ifndef CLIENT_CIPHER_H
#define CLIENT_CIPHER_H
#include <iostream>
// #include <streambuf>
#include <fstream>
#include <string>
#include <vector>
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <math.h>
#include <NTL/RR.h>
#include <util.h>
#include <param.h>
#include <client_util.h>
#include <string.h>
#include <helib/helib.h>
// #include <NTL/BasicThreadPool.h>

using namespace std;
using namespace NTL;

class Client_cipher{
    public:
    Client_cipher();
    void input_cipherGen(int *** A,int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W,const char *cipherfile,string padding);
    void dense_input_cipherGen(int* C, int lenth_C,string fc_inputcipherfile, helib::EncryptedArray ea);
    void input_cipherGen_Conv1(int ***AA, int inD, int inH, int inW, int kerN, int kerH, int kerW, int strid_H, int strid_W, const char *cipherfile, string padding);
    void input_cipherGen_Conv2(int***AA,int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W,int pool_strid_H, int pool_strid_W, const char *cipherfile,string padding,ZZX *t, ZZX** A,  ZZ_pX pm,  ZZ_pX p_temp1,  ZZ_pX p_temp2,  ZZ_pX p_temp3,  ZZ_pX p_temp4, ZZ q2d,enc_param ***usp);
    int*** DecCNN2D(string evalfile,string skfile,int kerN,int kerH,int kerW,int inH, int inW,int strid_H,int strid_W,string padding);
    int ***DecCNN2D_Conv1(int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W, int pool_stridH, int pool_stridW,  string evalfile, string skfile, string padding);
    int ***DecCNN3D_Conv2(int packedRowCount, int packedColCount, int packedRowUnit, int packedColUnit,  string evalfile, string skfile, int kerN, int kerH, int kerW, int inH, int inW, int strid_H, int strid_W,string padding);
    int*** DecCNN3DHaveAdd(string evalfile,string skfile,int kerN,int kerD,int kerH,int kerW,int inH, int inW,int strid_H,int strid_W,string padding);
    void dense_input_cipherGen(int* C, int lenth_C,string fc_inputcipherfile, ZZX *t, ZZX** A,  ZZ_pX pm,  ZZ_pX p_temp1,  ZZ_pX p_temp2,  ZZ_pX p_temp3,  ZZ_pX p_temp4, ZZ q2d, enc_param *usp);
    int* DecDenseHaveAdd_relu(int den_num1,int fc1_w,string resfile);
    double DecDenseHaveAdd_softmax(int den_num1,int fc1_w,string resfile);

};
#endif // CLIENT_CIPHER_H