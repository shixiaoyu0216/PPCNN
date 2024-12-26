#ifndef CLIENT_UTIL_H
#define CLIENT_UTIL_H
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
#include <util.h>
#include <helib/helib.h>
#include <sys/time.h>



using namespace std;
using namespace NTL;

class Client_util
{
public:
    Client_util();
    int** plain_input(int***AA, int inD, int inH, int inW, int kerH, int kerW, int strid_H, int strid_W);
    void read_data_Images(string filepath, vector<vector<float> >&a);
    void read_data_Label(string filename, vector<double>&labels);
    void read_Mnist_Label(string filename, vector<double>&labels);
    void read_Mnist_Images(string filename, vector<vector<double>>&images);   
    ZZ dec_one_poly(string skfile, ZZX h1, ZZX h2, ZZX h3);
    ZZX dec_onepackedPloy(string skfile, ZZX h1, ZZX h2, ZZX h3);
    ZZ* GetpakedPCoeff_bias(ZZX m1m2, int kerH, int kerW, int packed_num1);
    ZZ* GetpakedPCoeff_Nobias(ZZX m1m2, int kerH, int kerW, int packed_num2);
    int relu(float x);
    int relu(double x);
    ZZ relu2(float x);
    void matPacking(int **dat, int num_datrows, int num_datcols, int **packedMat, int num_rows, int  num_cols);
    void matPackingNoOne(int **dat, int num_datrows, int num_datcols, int **packedMat, int num_rows, int  num_cols);
    void slideToMat_SAME(int ***C, int ***A, int inD,int inH, int inW,int kerH,int kerW,int strid_H,int strid_W);
    void slideToMat_VALID(int ***C, int ***A, int inD,int inH, int inW,int kerH,int kerW,int strid_H,int strid_W);   
    void matPacking_by_poolsize(ZZX** polys, int **dat, int datrows, int datcols, int pcol);
    void matPacking_by_poolsize1(ZZX** polys, int **dat, int datrows, int datcols, int pcol);
    void matPacking_by_poolsize11(ZZX** polys, int **dat, int datrows, int datcols, int pcol);
    void matPacking_by_poolsize12(ZZX** polys, int **dat, int datrows, int datcols, int pcol);
    void packingbyPoolWin(ZZX *poly, int **dat, int datrows, int datcols, int pcol, int pool_stridH, int pool_stridW, int max_pool_num);

    int ***maxPool(int ***input,int kerN,int outH,int outW,int pool_sizeH,int pool_sizeW,int pool_stridH, int pool_stridW);
    int* flatten(int ***P, int kerN,int in2_H,int in2_W);
    void softmax(double *x,int num, double *a);
    void Encode_Enc(int ***AA, int inD, int inH, int inW, string cipherfile,helib:: EncryptedArray ea);
    void Encode_Enc(int *AA, int AA_len , string cipherfile,helib:: EncryptedArray ea);
    void Genl_Dec(ZZX *h5,int kerN,string cipherfile,string skfile,int slots,long p);
    int*** GetCoef_Relu(ZZX *h5,int kerN, int outH,int outW);
    int*** Decode_Relu(ZZX *h,int kerN, int outH,int outW,helib::EncryptedArray ea, long p);
    int*** Decode(ZZX *h,int kerN, int outH,int outW,helib::EncryptedArray ea, long p);
    int* Decode(ZZX h,int C_len,helib::EncryptedArray ea,long p);
    int* Decode_Relu(ZZX h,int C_len,helib::EncryptedArray ea, long p);
    int* dense_input(int* C, int lenth_C);
    void Encode_Enc(int ***AA, int inD, int inH, int inW, string cipherfile,helib:: EncryptedArray ea, ZZX *t, ZZX**A, ZZ_pX pm, ZZ_pX p_temp1, ZZ_pX p_temp2, ZZ_pX p_temp3, ZZ_pX p_temp4,ZZ q2d,  enc_param  ***usp);

    inline void file_to_string(vector<string> &record, const string& line, char delimiter)
    {
    int linepos=0;
    char c;
    int linemax=line.length();
    string curstring;
    record.clear();
    while(linepos<linemax)
    {
        c = line[linepos];
        // cout<<c;       //输出文件每行内容，分析后进行提取
        if(isdigit(c)||c=='.'){
            curstring+=c;
        }
        else if(c==delimiter&&curstring.size()){
            record.push_back(curstring);
            curstring="";
        }
        ++linepos;
    }
    cout<<'\n';
    if(curstring.size())
        record.push_back(curstring);
    return;
}
 
inline float string_to_float(string str){
    int i=0,len=str.length();
    float sum=0;
    while(i<len){
        if(str[i]=='.') break;
        sum=sum*10+str[i]-'0';
        ++i;
    }
    ++i;
    float t=1,d=1;
    while(i<len){
        d*=0.1;
        t=str[i]-'0';
        sum+=t*d;
        ++i;
    }
    return sum;
}


};
#endif // CLIENT_UTIL_H