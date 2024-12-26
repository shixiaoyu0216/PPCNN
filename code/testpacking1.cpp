
#include <iostream>
#include <math.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/RR.h>

using namespace std;
using namespace NTL;
int main(){
    int kerN = 6, kerH=5, kerW=5,  pcol=32,pool_stridH=2,pool_stridW=2;
    int inH = 32, strid_H = 1,inW =32, strid_W = 1;
    int n = 2048;
    int outH = ceil((inH*1.0)/strid_H);   
    int outW = ceil((inW * 1.0)/strid_W);      
    int datrows = kerH * kerW;
    int datcols = outH * outW; 
    int real_nrow = datrows + 1;            
    int max_packed_num = n/(real_nrow);
    int max_pool_num = max_packed_num / (pool_stridH * pool_stridW);
    int packed_num = max_pool_num * pool_stridH * pool_stridH; 
    int packed_col = ceil(double(1.0*datcols/packed_num));
    
    int part = ceil(double(datcols) / double(pcol)); //part是将dat的列要截断为32个部分
    int prow = part ;  
    int packedNum_row = pool_stridH; 
    int packedNum_col = pool_stridW;
    int packedColCount = ceil(pcol / (packedNum_col * 1.0)); 
    int packedRowCount = ceil(prow / (packedNum_row * 1.0));   

  
    int data[prow * pcol];
    for(int i = 0; i < prow * pcol; i++){
       data[i] = rand() % 100;
    } 
    //每col列一截断
    int Pdat1[prow][pcol];
    for(int i = 0; i < prow; i++){       
        for(int j = 0; j < pcol; j++){ 
            Pdat1[i][j] = data[i * pcol + j];
            cout<<Pdat1[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    //是否需要填0判断
    int pp_row = packedRowCount * packedNum_row;
    int pp_col = packedColCount * packedNum_col;
    int p_dat2[pp_row][pp_col];

    for(int i = 0; i < pp_row; i++){
        for(int j = 0; j < pp_col; j++){
            p_dat2[i][j] = 0;
        }
    }

    for(int i = 0; i < prow; i++){
        for(int j = 0; j < pcol; j++){
            p_dat2[i][j] = Pdat1[i][j];
        }
    }


    int p_dat3[packedRowCount*packedColCount][pool_stridW*pool_stridH];
    for (int row = 0; row < packedRowCount; row++){
        int startRowId = row * packedNum_row;
        int endRowId = (row + 1) * packedNum_row - 1;
        for (int col = 0; col < packedColCount; col++){
            int startColId = col * packedNum_col;
            int endColId = (col + 1) * packedNum_col - 1;
            int index1 = 0;
            for(int i = startRowId; i<=endRowId; i++){
                for(int j = startColId; j<=endColId; j++){
                    p_dat3[row*packedColCount+col][index1] = p_dat2[i][j];
                    index1++;
                }
            }
 
        }
    }
     
    int p_dat4[packed_col*max_pool_num][pool_stridW*pool_stridH];   
    for(int i = 0; i < packedRowCount * packedColCount; i++){
        for(int j = 0; j < pool_stridW * pool_stridH; j++){
            p_dat4[i][j] = p_dat3[i][j];        
        }
    }

    ZZX polys[packed_col];
    for(int z = 0; z < packed_col; z++){
        polys[z].SetLength(n);
        int index = 0;
        for(int i = 0; i < max_pool_num; i++){      
            for(int j = 0; j < pool_stridW * pool_stridH; j++){
                SetCoeff(polys[z], kerH*kerW+(kerH*kerW+1)*index, p_dat4[z*max_pool_num+i][j]);
                index++;
            }
        }
        cout<<polys[z]<<endl;
    }
    return 0;
}