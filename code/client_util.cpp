#include "client_util.h"
using namespace std;

Client_util::Client_util()
{
}

int **Client_util::plain_input(int ***AA, int inD, int inH, int inW, int kerH, int kerW, int strid_H, int strid_W)
{
    util U;
    Client_util cu;
    int outH = ceil(inH / strid_H);
    int outW = ceil(inW / strid_W);

    int ***C = U.newArray3(inD, kerH * kerW, outH * outW);
 

    cu.slideToMat_SAME(C, AA, inD, inH, inW, kerH, kerW, strid_H, strid_W); //test right

    int **CC = U.newArray2(inD * kerH * kerW, outH * outW);

    int dim3 = 0;
    int tmp = 0;
    for (int i = 0; i < inD * kerH * kerW; i++)
    {
        if (i > 0 && i % (kerH * kerW) == 0)
        {
            dim3 = dim3 + 1;
            tmp = tmp + 1;
        }
        int slides = tmp * kerH * kerW;
        for (int j = 0; j < inH * inW; j++)
        {
            CC[i][j] = C[dim3][i - slides][j];
        }
    }
  
    U.deleteArray3(C, inD, kerH * kerW, outH * outW);
    return CC;
}

void Client_util::matPacking_by_poolsize(ZZX **polys, int **dat, int datrows, int datcols, int pcol)
{
    util U;
    //dat:datrow*datcols; 加上bias 行数多一变为 tempMat:real_nrow*（pcol*part）i.e.26*784; Pmat:prow*pcol

    //  int datrows=25,datcols=784;//输入矩阵dat 的行和列

    int part = ceil(double(datcols) / double(pcol)); //part是将dat的列要截断为几个部分
    int real_nrow = datrows + 1;
    int prow = part * (datrows + 1); //将dat28列一截断生成的矩阵Pdat的行和列

    int **Pdat = U.newArray2(prow, pcol);
    U.initArray2(Pdat, prow, pcol, 0);

    //tempMat是第一行全为1，其余值dat矩阵的值大小为26*784
    int **tempMat = new int *[real_nrow];
    for (int i = 0; i < real_nrow; i++)
    {
        tempMat[i] = new int[datcols];
    }

    for (int i = 0; i < real_nrow; i++)
    {
        for (int j = 0; j < datcols; j++)
        {
            tempMat[i][j] = 0;
        }
    }
    for (int j = 0; j < datcols; j++)
    {
        tempMat[0][j] = 1;
    }

    for (int i = 1; i < real_nrow; i++)
    {
        for (int j = 0; j < datcols; j++)
        {
            tempMat[i][j] = dat[i - 1][j];
        }
    }

    // data每列倒序排列
    int value = 0;
    // cout <<value<<endl;
    int iters = real_nrow / 2;
    for (int i = 0; i < iters; i++)
    {
        for (int j = 0; j < datcols; j++)
        {
            value = tempMat[i][j];
            tempMat[i][j] = tempMat[real_nrow - 1 - i][j];
            tempMat[real_nrow - 1 - i][j] = value;
        }
    }

    int eleCeil = datrows + 1;       // 1:bias
    int packedNum_row = 2 * eleCeil; //26*2
    int packedNum_col = 4;
    // int unitNum = 2 * 4;
    int packedColCount = ceil(pcol / (packedNum_col * 1.0)); //7
    int packedRowCount = ceil(prow / (packedNum_row * 1.0)); //14

    int tempid = 0;
    for (int i = 0; i < part; i++)
    {
        int rowId = i * real_nrow;
        for (int j = 0; j < pcol; j++)
        {
            for (int k = 0; k < real_nrow; k++)
            {
                Pdat[rowId + k][j] = tempMat[k][tempid];
            }
            tempid++;
        }
    }

    for (int row = 0; row < packedRowCount; row++)
    {
        int startRowId = row * packedNum_row;
        int endRowId = (row + 1) * packedNum_row - 1;
        for (int col = 0; col < packedColCount; col++)
        {
            int startColId = col * packedNum_col;
            int endColId = (col + 1) * packedNum_col - 1;
            int index1 = 0;
            for (int j = startColId; j <= endColId; j++)
            {
                for (int i = startRowId; i < startRowId + eleCeil; i++)
                {
                    SetCoeff(polys[row][col], index1, Pdat[i][j]);
                    index1++;
                }
            }
            for (int j = startColId; j <= endColId; j++)
            {
                for (int i = startRowId + eleCeil; i <= endRowId; i++)
                {
                    SetCoeff(polys[row][col], index1, Pdat[i][j]);
                    index1++;
                }
            }
        }
    }

    // cout<<"packing result***************"<<endl;
    // for(int i=0; i<packedRowCount; i++){
    //     for(int j=0; j<packedColCount; j++){
    //         cout<<i<<"--"<<j<<":"<<polys[i][j]<<endl;
    //     }
    //     cout<<endl;
    // }
    //  cout<<polys[3][6]<<endl;

    for (int i = 0; i < real_nrow; i++)
    {
        delete[] tempMat[i];
    }
    delete[] tempMat;

    U.deleteArray2(Pdat, 26 * 28, 28);

    // return polys;
}
void Client_util::matPacking_by_poolsize1(ZZX **polys, int **dat, int datrows, int datcols, int pcol)
{
    util U;
    //dat:datrow*datcols; 加上bias 行数多一变为 tempMat:real_nrow*（pcol*part）i.e.26*784; Pmat:prow*pcol

    //  int datrows=36,datcols=49;//输入矩阵dat 的行和列

    int part = ceil(double(datcols) / double(pcol)); //part是将dat的列要截断为几个部分
    int real_nrow = datrows + 1;
    int prow = part * (datrows + 1); //将dat 7 列一截断生成的矩阵Pdat的行和列

    int **Pdat = U.newArray2(prow, pcol);
    U.initArray2(Pdat, prow, pcol, 0);

    //tempMat是第一行全为1，其余值dat矩阵的值大小为36*49
    int **tempMat = new int *[real_nrow];
    for (int i = 0; i < real_nrow; i++)
    {
        tempMat[i] = new int[datcols];
    }

    for (int i = 0; i < real_nrow; i++)
    {
        for (int j = 0; j < datcols; j++)
        {
            tempMat[i][j] = 0;
        }
    }
    for (int j = 0; j < datcols; j++)
    {
        tempMat[0][j] = 1;
    }

    for (int i = 1; i < real_nrow; i++)
    {
        for (int j = 0; j < datcols; j++)
        {
            tempMat[i][j] = dat[i - 1][j];
        }
    }

    // data每列倒序排列
    int value = 0;

    int iters = real_nrow / 2;
    for (int i = 0; i < iters; i++)
    {
        for (int j = 0; j < datcols; j++)
        {
            value = tempMat[i][j];
            tempMat[i][j] = tempMat[real_nrow - 1 - i][j];
            tempMat[real_nrow - 1 - i][j] = value;
        }
    }

    int eleCeil = datrows + 1;       // 1:bias
    int packedNum_row = 2 * eleCeil; //37*2
    int packedNum_col = 2;
    // int unitNum = 2 * 4;
    int packedColCount = ceil(pcol / (packedNum_col * 1.0)); //4
    int packedRowCount = ceil(prow / (packedNum_row * 1.0)); //4

    int tempid = 0;
    for (int i = 0; i < part; i++)
    {
        int rowId = i * real_nrow;
        for (int j = 0; j < pcol; j++)
        {
            for (int k = 0; k < real_nrow; k++)
            {
                Pdat[rowId + k][j] = tempMat[k][tempid];
            }
            tempid++;
        }
    }

    //是否需要填0 判断
    int pp_row = packedRowCount * packedNum_row;
    int pp_col = packedColCount * packedNum_col;
    int **p_dat2 = U.newArray2(pp_row, pp_col);
    U.initArray2(p_dat2, pp_row, pp_col, 0);
    for (int i = 0; i < prow; i++)
    {
        for (int j = 0; j < pcol; j++)
        {
            p_dat2[i][j] = Pdat[i][j];
        }
    }

    for (int row = 0; row < packedRowCount; row++)
    {
        int startRowId = row * packedNum_row;
        int endRowId = (row + 1) * packedNum_row - 1;
        for (int col = 0; col < packedColCount; col++)
        {
            int startColId = col * packedNum_col;
            int endColId = (col + 1) * packedNum_col - 1;
            int index1 = 0;

            for (int j = startColId; j <= endColId; j++)
            {
                for (int i = startRowId; i < startRowId + eleCeil; i++)
                {
                    SetCoeff(polys[row][col], index1, p_dat2[i][j]);
                    index1++;
                }
            }

            for (int j = startColId; j <= endColId; j++)
            {
                for (int i = startRowId + eleCeil; i <= endRowId; i++)
                {
                    SetCoeff(polys[row][col], index1, p_dat2[i][j]);
                    index1++;
                }
            }
        }
    }

    // cout<<"packing result***************"<<endl;
    // for(int i=0; i<packedRowCount; i++){
    //     for(int j=0; j<packedColCount; j++){
    //         cout<<i<<"--"<<j<<":"<<polys[i][j]<<endl;
    //     }
    //     cout<<endl;
    // }
    //  cout<<polys[3][6]<<endl;

    for (int i = 0; i < real_nrow; i++)
    {
        delete[] tempMat[i];
    }
    delete[] tempMat;

    U.deleteArray2(Pdat, prow, pcol);
    U.deleteArray2(p_dat2, pp_row, pp_col);
    // return polys;
}

void Client_util::packingbyPoolWin(ZZX *poly, int **dat, int datrows, int datcols, int pcol, int pool_stridH, int pool_stridW, int max_pool_num)
{
    int part = ceil(double(datcols) / double(pcol)); //part是将dat的列要截断为32个部分
    int real_nrow = datrows + 1;
    int prow = part *  real_nrow;   //将dat 32 列一截断生成的矩阵Pdat的行和列   

    //tempMat是第一行全为1
    int **tempMat = new int *[real_nrow];
    for (int i = 0; i < real_nrow; i++){
        tempMat[i] = new int[datcols];
    }  
    for (int i = 0; i < real_nrow; i++){
        for (int j = 0; j < datcols; j++){
            tempMat[i][j] = 0;
        }
    }
    for (int j = 0; j < datcols; j++){
        tempMat[0][j] = 1;
    }  
    for (int i = 1; i < real_nrow; i++){
        for (int j = 0; j < datcols; j++){
            tempMat[i][j] =  dat[i - 1][j];            
        }
    }

    // data每列倒序排列
    int value = 0;
    int iters = real_nrow / 2;
    for (int i = 0; i < iters; i++)
    {
        for (int j = 0; j < datcols; j++)
        {
            value = tempMat[i][j];
            tempMat[i][j] = tempMat[real_nrow - 1 - i][j];
            tempMat[real_nrow - 1 - i][j] = value;
        }
    }
    
    //每col列一截断
    int tempid = 0;
    int Pdat1[prow][pcol];
    for(int i = 0; i < part; i++){
        int rowId = i*real_nrow;
        for(int j = 0; j < pcol; j++){
            for(int k = 0; k < real_nrow; k++){
                Pdat1[rowId+k][j] = tempMat[k][tempid];
            }
            tempid++;           
        }
    }


      
    int packedNum_row = pool_stridH * real_nrow; 
    int packedNum_col = pool_stridW;
    int packedColCount = ceil(pcol / (packedNum_col * 1.0)); 
    int packedRowCount = ceil(prow / (packedNum_row * 1.0));   

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


    int p_dat3[packedRowCount*packedColCount][real_nrow*pool_stridW*pool_stridH]={};
    for (int row = 0; row < packedRowCount; row++){
        int startRowId = row * packedNum_row;
        int endRowId = (row + 1) * packedNum_row - 1;
        for (int col = 0; col < packedColCount; col++){
            int startColId = col * packedNum_col;
            int endColId = (col + 1) * packedNum_col - 1;
            int index1 = 0;
        
            for (int j = startColId; j <= endColId; j++){
                for (int i = startRowId; i <startRowId+real_nrow; i++){
                    p_dat3[row*packedColCount+col][index1] = p_dat2[i][j];
                    index1++;
                }
            }
          
            for(int j = startColId; j <= endColId; j++){
                for(int i = startRowId + real_nrow; i <= endRowId; i++){
                    p_dat3[row*packedColCount+col][index1] = p_dat2[i][j];
                    index1++;
                }
            }
        }
    } 

    
    for (int i = 0; i < real_nrow; i++)
    {
        delete[] tempMat[i];
    }
    delete[] tempMat;


    int poly_num = ceil(1.0 * packedRowCount * packedColCount /max_pool_num);//2
    int p_dat4[poly_num*max_pool_num][pool_stridW*pool_stridH*real_nrow]={};

   
    for(int i = 0; i < packedRowCount*packedColCount; i++){
        for(int j = 0; j < pool_stridW*pool_stridH*real_nrow; j++){
            p_dat4[i][j] = p_dat3[i][j];        
        }

    }


    for(int z = 0; z < poly_num; z++){
        poly[z].SetLength(n);
        int index = 0;
        for(int i = 0; i < max_pool_num; i++){      
            for(int j = 0; j < pool_stridW*pool_stridH*real_nrow; j++){
                SetCoeff(poly[z], index, p_dat4[z*max_pool_num+i][j]);
                index++;
            }
        }
    }
}








/*
扩展矩阵的每4列和2行打包为一个多项式，且扩展矩阵中包含1*/
void Client_util::matPacking_by_poolsize11(ZZX **polys, int **dat, int datrows, int datcols, int pcol)
{
    util U;
    //dat:datrow*datcols; 加上bias 行数多一变为 tempMat:real_nrow*（pcol*part）i.e.26*784; Pmat:prow*pcol

    //  int datrows=36,datcols=49;//输入矩阵dat 的行和列

    int part = ceil(double(datcols) / double(pcol)); //part是将dat的列要截断为几个部分
    int real_nrow = datrows + 1;
    int prow = part * (datrows + 1); //将dat 7 列一截断生成的矩阵Pdat的行和列

    int **Pdat = U.newArray2(prow, pcol);
    U.initArray2(Pdat, prow, pcol, 0);
    //tempMat是第一行全为1，其余值dat矩阵的值大小为36*49
    int **tempMat = new int *[real_nrow];
    for (int i = 0; i < real_nrow; i++)
    {
        tempMat[i] = new int[datcols];
    }

    for (int i = 0; i < real_nrow; i++)
    {
        for (int j = 0; j < datcols; j++)
        {
            tempMat[i][j] = 0;
        }
    }
    for (int j = 0; j < datcols; j++)
    {
        tempMat[0][j] = 1;
    }

    for (int i = 1; i < real_nrow; i++)
    {
        for (int j = 0; j < datcols; j++)
        {
            tempMat[i][j] = dat[i - 1][j];
        }
    }

    // data每列倒序排列
    int value = 0;

    int iters = real_nrow / 2;
    for (int i = 0; i < iters; i++)
    {
        for (int j = 0; j < datcols; j++)
        {
            value = tempMat[i][j];
            tempMat[i][j] = tempMat[real_nrow - 1 - i][j];
            tempMat[real_nrow - 1 - i][j] = value;
        }
    }

    int eleCeil = datrows + 1; // 1:bias
    int packedNum_row = 2 * eleCeil;
    int packedNum_col = 4;
    // int unitNum = 2 * 4;
    int packedColCount = ceil(pcol / (packedNum_col * 1.0));
    int packedRowCount = ceil(prow / (packedNum_row * 1.0));

    int tempid = 0;
    for (int i = 0; i < part; i++)
    {
        int rowId = i * real_nrow;
        for (int j = 0; j < pcol; j++)
        {
            for (int k = 0; k < real_nrow; k++)
            {
                Pdat[rowId + k][j] = tempMat[k][tempid];
            }
            tempid++;
        }
    }

    //是否需要填0 判断
    int pp_row = packedRowCount * packedNum_row;
    int pp_col = packedColCount * packedNum_col;
    int **p_dat2 = U.newArray2(pp_row, pp_col);
    U.initArray2(p_dat2, pp_row, pp_col, 0);
    for (int i = 0; i < prow; i++)
    {
        for (int j = 0; j < pcol; j++)
        {
            p_dat2[i][j] = Pdat[i][j];
        }
    }

    for (int row = 0; row < packedRowCount; row++)
    {
        int startRowId = row * packedNum_row; //row*2
        int endRowId = (row + 1) * packedNum_row - 1;
        for (int col = 0; col < packedColCount; col++)
        {
            int startColId = col * packedNum_col;
            int endColId = (col + 1) * packedNum_col - 1;
            int index1 = 0;

            for (int j = startColId; j <= endColId; j++)
            {
                for (int i = startRowId; i < startRowId + eleCeil; i++)
                {
                    SetCoeff(polys[row][col], index1, p_dat2[i][j]);
                    index1++;
                }
            }

            for (int j = startColId; j <= endColId; j++)
            {
                for (int i = startRowId + eleCeil; i <= endRowId; i++)
                {
                    SetCoeff(polys[row][col], index1, p_dat2[i][j]);
                    index1++;
                }
            }
        }
    }

    // cout<<"packing result***************"<<endl;
    // for(int i=0; i<packedRowCount; i++){
    //     for(int j=0; j<packedColCount; j++){
    //         cout<<i<<"--"<<j<<":"<<polys[i][j]<<endl;
    //     }
    //     cout<<endl;
    // }

    for (int i = 0; i < real_nrow; i++)
    {
        delete[] tempMat[i];
    }
    delete[] tempMat;

    U.deleteArray2(Pdat, prow, pcol);
    U.deleteArray2(p_dat2, pp_row, pp_col);
    // return polys;
}






/*
扩展矩阵的每4列打包为一个多项式，且扩展矩阵中不包含1*/
void Client_util::matPacking_by_poolsize12(ZZX **polys, int **dat, int datrows, int datcols, int pcol)
{
    util U;
    int part = ceil(double(datcols) / double(pcol)); //part是将dat的列要截断为几个部分
    int real_nrow = datrows;
    int prow = part * datrows; //将dat 7 列一截断生成的矩阵Pdat的行和列
    int **Pdat = U.newArray2(prow, pcol);
    U.initArray2(Pdat, prow, pcol, 0);
    // data每列倒序排列
    int value = 0;
    int iters = real_nrow / 2;
    for (int i = 0; i < iters; i++)
    {
        for (int j = 0; j < datcols; j++)
        {
            value = dat[i][j];
            dat[i][j] = dat[real_nrow - 1 - i][j];
            dat[real_nrow - 1 - i][j] = value;
        }
    }
    // for (int i = 0; i < real_nrow; i++){
    //     for (int j = 0; j < datcols; j++){
    //         cout<<dat[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    // cout<<"每part=4个部分一截断"<<endl;
    int eleCeil = datrows;
    int packedNum_row = 2 * eleCeil; //37*2
    int packedNum_col = 4;
    int packedColCount = ceil(pcol / (packedNum_col * 1.0)); //4
    int packedRowCount = ceil(prow / (packedNum_row * 1.0)); //4
    // ZZX **polys = U.newZZXArray2(packedRowCount, packedColCount);//6*4*4

    int tempid = 0;
    for (int i = 0; i < part; i++)
    {
        int rowId = i * real_nrow;
        for (int j = 0; j < pcol; j++)
        {
            for (int k = 0; k < real_nrow; k++)
            {
                Pdat[rowId + k][j] = dat[k][tempid];
            }
            tempid++;
        }
    }
    //  cout<<"Pdat="<<endl;
    //     for(int i=0;i<prow;i++){
    //         for(int j=0;j<pcol;j++){
    //             cout<<Pdat[i][j]<<" ";
    //         }
    //         cout<<endl;
    //     }
    //     cout<<endl;

    //是否需要填0 判断
    int pp_row = packedRowCount * packedNum_row;
    int pp_col = packedColCount * packedNum_col;
    int **p_dat2 = U.newArray2(pp_row, pp_col);
    U.initArray2(p_dat2, pp_row, pp_col, 0);
    for (int i = 0; i < prow; i++)
    {
        for (int j = 0; j < pcol; j++)
        {
            p_dat2[i][j] = Pdat[i][j];
        }
    }

    for (int row = 0; row < packedRowCount; row++)
    {
        int startRowId = row * packedNum_row;
        int endRowId = (row + 1) * packedNum_row - 1;
        for (int col = 0; col < packedColCount; col++)
        {
            int startColId = col * packedNum_col;
            int endColId = (col + 1) * packedNum_col - 1;
            int index1 = 0;
            for (int j = startColId; j <= endColId; j++)
            {
                for (int i = startRowId; i < startRowId + eleCeil; i++)
                {
                    SetCoeff(polys[row][col], index1, p_dat2[i][j]);
                    index1++;
                }
            }

            for (int j = startColId; j <= endColId; j++)
            {
                for (int i = startRowId + eleCeil; i <= endRowId; i++)
                {
                    SetCoeff(polys[row][col], index1, p_dat2[i][j]);
                    index1++;
                }
            }
        }
    }

    // cout<<"packing result***************"<<endl;
    // for(int i=0; i<packedRowCount; i++){
    //     for(int j=0; j<packedColCount; j++){
    //         cout<<i<<"--"<<j<<":"<<polys[i][j]<<endl;
    //     }
    //     cout<<endl;
    // }
    //  cout<<polys[3][6]<<endl;
    U.deleteArray2(Pdat, prow, pcol);
    U.deleteArray2(p_dat2, pp_row, pp_col);
}

// client_enc_param Client_util::randnum_gen(){
//     client_enc_param usp;
// 	for(int i=0;i<k;i++){
//         usp.r1[i].SetLength(n);
//         usp.Bin_dis.Binomial_distribution_gen(usp.r1[i],n,ita);
//         }
//     for(int i=0;i<k;i++){
//        usp.e11[i].SetLength(n);
//        usp.Bin_dis.Binomial_distribution_gen(usp.e11[i],n,ita);
//     }
//     usp.Bin_dis.Binomial_distribution_gen(usp.e12,n,ita);
//     return usp;
// }

void Client_util::read_Mnist_Label(string filename, vector<double> &labels)
{
    ifstream file(filename, ios::binary);
    util U;
    if (file.is_open())
    {
        int magic_number = 0;
        int number_of_images = 0;
        file.read((char *)&magic_number, sizeof(magic_number));
        file.read((char *)&number_of_images, sizeof(number_of_images));
        magic_number = U.ReverseInt(magic_number);
        number_of_images = U.ReverseInt(number_of_images);
        //cout << "magic number = " << magic_number << endl;
        //cout << "number of images = " << number_of_images << endl;
        for (int i = 0; i < number_of_images; i++)
        {
            unsigned char label = 0;
            file.read((char *)&label, sizeof(label));
            labels.push_back((double)label);
        }
    }
}

// vector< vector< long  > > A( k, vector< vector< long long > >( k, vector< long long>( n, 0 )));
//     int temperature=0;
//     ifstream in_file_pk_B(file_B_path,ios::binary);
//     for(int i = 0;i<k;i++){
//         for(int j = 0;j<k;j++){
//             for(int g=0;g<n;g++){
//                 in_file_pk_B.read((char*)&temperature,sizeof(int));
//                 A[i][j][g]=temperature;
//             }
//         }
//     }

void Client_util::read_data_Images(string filepath, vector<vector<float>> &a)
{
    vector<string> row;
    string line;
    string filename;
    ifstream in(filepath); //读取工程文件下的.csv文件
    if (in.fail())
    {
        cout << "File not found" << endl;
        exit(0);
    } //读取失败，返回为假

    vector<float> b;
    while (getline(in, line) && in.good())
    {
        file_to_string(row, line, ','); //把line里的单元格数字字符提取出来，“,”为单元格分隔符
        for (int i = 0, leng = row.size(); i < leng; i++)
        {
            // cout << row[i] << " ";
            b.push_back(string_to_float(row[i]));
        }
        a.push_back(b); //二维数组a存储读入变量
        b.clear();
    }
    in.close();
    //  ofstream out("file1.csv");  //写入到文件file1
    //  for(int i=0;i<a.size();++i){
    //     for(int j=0;j<a[0].size();++j){
    //         cout<<a[i][j]<<',';
    //     }
    //     cout<<'\n';
    // }
    // out.close();
}

void Client_util::read_Mnist_Images(string filename, vector<vector<double>> &images)
{
    util U;
    ifstream file(filename, ios::binary);
    if (file.is_open())
    {
        int magic_number = 0;
        int number_of_images = 0;
        int n_rows = 0;
        int n_cols = 0;
        unsigned char label;
        file.read((char *)&magic_number, sizeof(magic_number));
        file.read((char *)&number_of_images, sizeof(number_of_images));
        file.read((char *)&n_rows, sizeof(n_rows));
        file.read((char *)&n_cols, sizeof(n_cols));
        magic_number = U.ReverseInt(magic_number);
        number_of_images = U.ReverseInt(number_of_images);
        n_rows = U.ReverseInt(n_rows);
        n_cols = U.ReverseInt(n_cols);
        //cout << "magic number = " << magic_number << endl;
        //cout << "number of images = " << number_of_images << endl;
        //cout << "rows = " << n_rows << endl;
        //cout << "cols = " << n_cols << endl;

        for (int i = 0; i < number_of_images; i++)
        {
            vector<double> tp;
            for (int r = 0; r < n_rows; r++)
            {
                for (int c = 0; c < n_cols; c++)
                {
                    unsigned char image = 0;
                    file.read((char *)&image, sizeof(image));
                    tp.push_back(image);
                }
            }
            images.push_back(tp);
        }
    }
}

ZZX Client_util::dec_onepackedPloy(string skfile, ZZX h1, ZZX h2, ZZX h3)
{
    util U;
    //read sk.cvs
    ZZX s[k];
    ifstream inFile2;
    inFile2.open(skfile, ios::in);
    if (!inFile2)
    {
        cout << "Cannot open sk in dec one ploy!" << endl;
        exit(0);
    }
    for (int i = 0; i < k; i++)
    {
        inFile2 >> s[i];
    }
    //cout<<endl;
    inFile2.close();

    ZZX s0u10w, s1u11w;
    ZZ_pX p_temp1, p_temp2, p_temp3, p_temp4;

    ZZ_p::init(q2);
    ZZ_pX ppm;
    SetCoeff(ppm, 0, 1);
    SetCoeff(ppm, n, 1);
    ZZ_pXModulus pm(ppm);

    ZZ_pX pp_temp1, pp_temp2, pp_temp3, pp_temp4;
    U.Mult(s0u10w, pm, s[0], h1, n, q2, pp_temp1, pp_temp2, pp_temp3, pp_temp4);
    U.Mult(s1u11w, pm, s[1], h2, n, q2, pp_temp1, pp_temp2, pp_temp3, pp_temp4);
    ZZX h4 = h3 - s0u10w - s1u11w;

    ZZ two_d;
    RR two_d_rr;
    power2(two_d, dp);
    two_d_rr = to_RR(two_d);
    RR q2_rr;
    q2_rr = to_RR(q);

    RR two_dp_rr;
    power2(two_dp_rr, dp);
 
    ZZ two_d2;
    ZZ two_d3;
    two_d2 = two_d / 2;
    two_d3 = two_d2 * (-1);
  

    for (int i = 0; i < n; i++)
    {
        RR coeff_rr;
        ZZ coeff_zz;
        GetCoeff(coeff_zz, h4, i);
        coeff_rr = to_RR(coeff_zz);
        coeff_zz = RoundToZZ(coeff_rr * two_dp_rr / q2_rr) % two_d;
        if (coeff_zz >= two_d2)
        {
            coeff_zz = coeff_zz - two_d;
        }
        else if (coeff_zz <= two_d3)
        {
            coeff_zz = coeff_zz + two_d;
        }
        SetCoeff(h4, i, coeff_zz);
    }
    return h4;
}

ZZ Client_util::dec_one_poly(string skfile, ZZX h1, ZZX h2, ZZX h3)
{
    util U;
    //read sk.cvs
    ZZX s[k];
    ifstream inFile2;
    inFile2.open(skfile, ios::in);
    for (int i = 0; i < k; i++)
    {
        inFile2 >> s[i];
        //cout<<s[i];
    }
    //cout<<endl;
    inFile2.close();

    ZZX s0u10w, s1u11w;
    ZZ_pX p_temp1, p_temp2, p_temp3, p_temp4;

    ZZ_p::init(q2);
    ZZ_pX ppm;
    SetCoeff(ppm, 0, 1);
    SetCoeff(ppm, n, 1);
    ZZ_pXModulus pm(ppm);

    ZZ_pX pp_temp1, pp_temp2, pp_temp3, pp_temp4;
    U.Mult(s0u10w, pm, s[0], h2, n, q2, pp_temp1, pp_temp2, pp_temp3, pp_temp4);
    U.Mult(s1u11w, pm, s[1], h3, n, q2, pp_temp1, pp_temp2, pp_temp3, pp_temp4);
    ZZX h4 = h1 - s0u10w - s1u11w;

    ZZ two_d;
    RR two_d_rr;
    power2(two_d, dp);
    two_d_rr = to_RR(two_d);
    RR q2_rr;
    q2_rr = to_RR(q);

    RR two_dp_rr;
    power2(two_dp_rr, dp);
 
    ZZ two_d2;
    ZZ two_d3;
    two_d2 = two_d / 2;
    two_d3 = two_d2 * (-1);
   
    for (int i = 0; i < n; i++)
    {
        RR coeff_rr;
        ZZ coeff_zz;
        GetCoeff(coeff_zz, h4, i);
        coeff_rr = to_RR(coeff_zz);
        coeff_zz = RoundToZZ(coeff_rr * two_dp_rr / q2_rr) % two_d;
        if (coeff_zz >= two_d2)
        {
            coeff_zz = coeff_zz - two_d;
        }
        else if (coeff_zz <= two_d3)
        {
            coeff_zz = coeff_zz + two_d;
        }
        SetCoeff(h4, i, coeff_zz);
    }
    ZZ inner;
    GetCoeff(inner, h4, 0);
    return inner;
    //cout<<inner<<" ";
}

ZZ *Client_util::GetpakedPCoeff_bias(ZZX m1m2, int kerH, int kerW, int packed_num1)
{
    util U;
    ZZ *inner = new ZZ[packed_num1];
    for (int z = 0; z < packed_num1; z++)
    {
        GetCoeff(inner[z], m1m2, kerH * kerW + (kerH * kerW + 1) * z);
    }

    return inner;
}

ZZ *Client_util::GetpakedPCoeff_Nobias(ZZX m1m2, int kerH, int kerW, int packed_num2)
{

    ZZ *inner = new ZZ[packed_num2];
    for (int z = 0; z < packed_num2; z++)
        GetCoeff(inner[z], m1m2, kerH * kerW - 1 + kerH * kerW * z);
    return inner;
}

int Client_util::relu(float x)
{
    float y;
    int z;
    if (x < 0)
        y = 0;
    else
        y = x;
    z = (int)y;
    return z;
}
int Client_util::relu(double x)
{
    util U;
    double y;
    int z;

    if (x < 0)
        y = 0;
    else
        y = x;
    z = U.doubleToInt(y);
    return z;
}

ZZ Client_util::relu2(float x)
{

    ZZ y;
    float z;
    if (x < 0)
        z = 0;
    else
        z = x;
    y = conv<ZZ>(z);
    return y;
}

void Client_util::softmax(double *x, int num, double *a)
{
    // for(int i=0;i<num;i++)
    //     cout<<x[i]<<" ";
    // cout<<endl;
    double sum = 0.0;
    for (int i = 0; i < num; i++)
    {
        x[i] = exp(x[i]);
        sum += x[i];
    }
    // cout<<"sum value  "<<sum<<endl;
    for (int i = 0; i < num; i++)
    {
        //  cout<<i<<"   "<<x[i];
        if (sum == 0)
        {
            a[i] = 0;
        }
        else
        {
            a[i] = x[i] / sum;
        }
        //  cout<<"  "<< a[i]<<endl;
    }
    // cout<<endl;
}

/**
 * dat: original matrix
 * packedMat: the return matrix with sizeof (num_rows, num_cols), the inital value is 0
 * */
void Client_util::matPacking(int **dat, int num_datrows, int num_datcols, int **packedMat, int num_rows, int num_cols)
{
    util U;
    int real_nrow = num_datrows + 1;
    int **tempMat = new int *[real_nrow];
    for (int i = 0; i < real_nrow; i++)
    {
        tempMat[i] = new int[num_datcols];
    }
    for (int j = 0; j < num_datcols; j++)
    {
        tempMat[0][j] = 1;
    }
    // cout<<"dat"<<endl;
    //     for(int i=0; i<num_datrows; i++){
    //     // for(int j=0; j<num_datcols; j++){
    //         cout<<dat[i][783]<<" ";
    //     // }
    //     // cout<<endl;
    // }
    // cout<<endl;
    for (int i = 1; i < real_nrow; i++)
    {
        for (int j = 0; j < num_datcols; j++)
        {
            tempMat[i][j] = dat[i - 1][j];
        }
    }
    //         cout <<"tempMat"<<endl;
    // for(int i=0; i<real_nrow; i++){
    //     // for(int j=0; j<num_datcols; j++){
    //         cout<<tempMat[i][num_datcols-1]<<" ";
    //     // }
    //     // cout<<endl;
    // }
    // cout<<endl;

    // for(int i=0; i<num_rows; i++){
    //     // for(int j=0; j<num_datcols; j++){
    //         cout<<packedMat[i][num_cols-1]<<" ";
    //     // }
    //     // cout<<endl;
    // }
    // cout<<endl;
    int value = 999;
    // cout <<value<<endl;
    int iters = real_nrow / 2;
    for (int i = 0; i < iters; i++)
    {
        for (int j = 0; j < num_datcols; j++)
        {
            value = tempMat[i][j];
            tempMat[i][j] = tempMat[real_nrow - 1 - i][j];
            tempMat[real_nrow - 1 - i][j] = value;
        }
    }
    // packing
    int iter_cols = num_rows / real_nrow;
    int iter_rows = num_cols - 1;
    int remain_cols = num_datcols % iter_cols;
    int id = 0;
    int count = 0;
    //  cout<<"iter_cols;"<<iter_cols<<"  iter_rows;"<< iter_rows<<"  reian;"<<remain_cols<<endl;

    for (int j = 0; j < iter_rows; j++)
    {
        count = 0;
        for (int i = 0; i < iter_cols; i++)
        {
            id = j * iter_cols + i;
            for (int idr = 0; idr < real_nrow; idr++)
            {
                packedMat[count][j] = tempMat[idr][id];
                count++;
            }
        }
    }
    count = 0;
    for (int i = remain_cols; i > 0; i--)
    {
        for (int j = 0; j < real_nrow; j++)
        {
            packedMat[count][num_cols - 1] = tempMat[j][num_datcols - i];
            count++;
        }
    }
    // cout<<"packMat"<<endl;
    // cout<<"count value:"<<count<<endl;
    //  for(int i=0; i<num_rows; i++){
    //     // for(int j=0; j<num_cols; j++){
    //         cout<<packedMat[i][num_cols-1]<<" ";
    //     // }
    //     // cout<<endl;
    // }
    //cout<<endl;
    for (int i = 0; i < real_nrow; i++)
    {
        delete[] tempMat[i];
    }
    delete[] tempMat;
}
/**
 * dat: original matrix
 * packedMat: the return matrix with sizeof (num_rows, num_cols), the inital value is 0
 * */
void Client_util::matPackingNoOne(int **dat, int num_datrows, int num_datcols, int **packedMat, int num_rows, int num_cols)
{
    // util U;
    int real_nrow = num_datrows;
    int **tempMat = new int *[real_nrow];
    for (int i = 0; i < real_nrow; i++)
    {
        tempMat[i] = new int[num_datcols];
    }

    // cout<<"dat"<<endl;
    //     for(int i=0; i<num_datrows; i++){
    //     // for(int j=0; j<num_datcols; j++){
    //         cout<<dat[i][783]<<" ";
    //     // }
    //     // cout<<endl;
    // }
    // cout<<endl;
    for (int i = 0; i < real_nrow; i++)
    {
        for (int j = 0; j < num_datcols; j++)
        {
            tempMat[i][j] = dat[i][j];
        }
    }
    //         cout <<"tempMat"<<endl;
    // for(int i=0; i<real_nrow; i++){
    //     // for(int j=0; j<num_datcols; j++){
    //         cout<<tempMat[i][num_datcols-1]<<" ";
    //     // }
    //     // cout<<endl;
    // }
    // cout<<endl;

    // for(int i=0; i<num_rows; i++){
    //     // for(int j=0; j<num_datcols; j++){
    //         cout<<packedMat[i][num_cols-1]<<" ";
    //     // }
    //     // cout<<endl;
    // }
    // cout<<endl;
    int value = 999;
    // cout <<value<<endl;
    int iters = real_nrow / 2;
    for (int i = 0; i < iters; i++)
    {
        for (int j = 0; j < num_datcols; j++)
        {
            value = tempMat[i][j];
            tempMat[i][j] = tempMat[real_nrow - 1 - i][j];
            tempMat[real_nrow - 1 - i][j] = value;
        }
    }
    // packing
    int iter_cols = num_rows / real_nrow;
    int iter_rows = num_cols - 1;
    int remain_cols = num_datcols % iter_cols;
    int id = 0;
    int count = 0;
    //  cout<<"iter_cols;"<<iter_cols<<"  iter_rows;"<< iter_rows<<"  reian;"<<remain_cols<<endl;

    for (int j = 0; j < iter_rows; j++)
    {
        count = 0;
        for (int i = 0; i < iter_cols; i++)
        {
            id = j * iter_cols + i;
            for (int idr = 0; idr < real_nrow; idr++)
            {
                packedMat[count][j] = tempMat[idr][id];
                count++;
            }
        }
    }
    count = 0;
    for (int i = remain_cols; i > 0; i--)
    {
        for (int j = 0; j < real_nrow; j++)
        {
            packedMat[count][num_cols - 1] = tempMat[j][num_datcols - i];
            count++;
        }
    }
    // cout<<"packMat"<<endl;
    // cout<<"count value:"<<count<<endl;
    //  for(int i=0; i<num_rows; i++){
    //     // for(int j=0; j<num_cols; j++){
    //         cout<<packedMat[i][num_cols-1]<<" ";
    //     // }
    //     // cout<<endl;
    // }
    //cout<<endl;
    for (int i = 0; i < real_nrow; i++)
    {
        delete[] tempMat[i];
    }
    delete[] tempMat;
}

void Client_util::slideToMat_SAME(int ***C, int ***A, int inD, int inH, int inW, int kerH, int kerW, int strid_H, int strid_W)
{
    util U;
    const int outH = ceil((inH * 1.0) / strid_H);
    const int outW = ceil((inW * 1.0) / strid_W);
    const int pad_W = max((outW - 1) * strid_W + kerW - inW, 0);
    const int pad_H = max((outH - 1) * strid_H + kerH - inH, 0);
    const int pad_top = int((pad_H) / 2);
    const int pad_left = int((pad_W) / 2);

    int ***B = U.newArray3(inD, inH + pad_H, inW + pad_W);

    U.initArray3(B, inD, inH + pad_H, inW + pad_W, 0);
    //pading zeros to A become B
    for (int z = 0; z < inD; z++)
    {
        int i = 0;
        for (int u = 0; u < inH; u++)
        {
            int j = 0;
            for (int v = 0; v < inW; v++)
            {
                B[z][i + pad_top][j + pad_left] = A[z][u][v];
                j = j + 1;
            }
            i = i + 1;
        }
    }

    // //cout A and B
    // cout<<"A=";
    // for (int i=0;i<inH;i++){
    //     for(int j=0;j<inW;j++){
    //         cout<<A[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    // cout<<"B=";
    // for (int i=0; i<inH+pad_H; i++){
    //     for(int j=0; j<inW+pad_W; j++){
    //         cout << B[0][i][j] <<" ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    //after padding,slide window on B and transform to C
    for (int z = 0; z < inD; z++)
    {
        int y = 0;
        for (int i = 0; i < outH; i++)
        {
            int start_u = strid_H * i;
            for (int j = 0; j < outW; j++)
            {
                int x = 0;
                int start_v = strid_W * j;
                for (int u = start_u; u < start_u + kerH; u++)
                {
                    for (int v = start_v; v < start_v + kerW; v++)
                    {
                        C[z][x][y] = B[z][u][v];
                        x = x + 1;
                    }
                }
                y = y + 1;
            }
        }
    }
    U.deleteArray3(B, inD, inH + pad_H, inW + pad_W);
}

void Client_util::slideToMat_VALID(int ***C, int ***A, int inD, int inH, int inW, int kerH, int kerW, int strid_H, int strid_W)
{
    // util U;
    int outH = ceil(float(inH - kerH + 1) / float(strid_H));
    int outW = ceil(float(inW - kerW + 1) / float(strid_W));
    //after padding,slide window on A and transform to C
    for (int z = 0; z < inD; z++)
    {
        int y = 0;
        for (int i = 0; i < outH; i++)
        {
            for (int j = 0; j < outW; j++)
            {
                int x = 0;
                for (int u = 0; u < kerH; u++)
                    for (int v = 0; v < kerW; v++)
                    {
                        C[z][x][y] = A[z][u + i][v + j];
                        x = x + 1;
                    }
                y = y + 1;
            }
        }
    }
}

int ***Client_util::maxPool(int ***input, int kerN, int outH, int outW,
                            int pool_sizeH, int pool_sizeW, int pool_stridH, int pool_stridW)
{
    util U;
    //padding
    int in2_H = ceil(outH * 1.0 / pool_stridH);
    int in2_W = ceil(outW * 1.0 / pool_stridW);
    int pad2_H = max((in2_H - 1) * pool_stridH + pool_sizeH - outH, 0);
    int pad2_W = max((in2_W - 1) * pool_stridW + pool_sizeW - outW, 0);
    int pad2_top = int((pad2_H) / 2);
    int pad2_left = int((pad2_W) / 2);

    //B is padding 0 in P
    int B[kerN][outH + pad2_H][outW + pad2_W];
    for (int z = 0; z < kerN; z++)
    {
        for (int i = 0; i < outH + pad2_H; i++)
        {
            for (int j = 0; j < outW + pad2_W; j++)
            {
                B[z][i][j] = 0;
            }
        }
    }
    for (int z = 0; z < kerN; z++)
    {
        int i = 0;
        for (int u = 0; u < outH; u++)
        {
            int j = 0;
            for (int v = 0; v < outW; v++)
            {
                B[z][i + pad2_top][j + pad2_left] = input[z][u][v];
                j = j + 1;
            }
            i = i + 1;
        }
    }
    /*      
    cout<<"Test padding..."<<endl;
    cout<<"B[31][i][j]="<<endl;
    for (int i=0;i<outH+pad2_H;i++){
        for (int j=0;j<outW+pad2_W;j++){
            cout<<B[31][i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;   */

    //maxpool in every fearture map
    int ***P_res;

    P_res = U.newArray3(kerN, in2_H, in2_W);
    // cout<<kerN<<","<<in2_H<<","<<in2_W<<endl;
    int max;
    for (int z = 0; z < kerN; z++)
    {
        int u1 = 0, u2 = 0;
        for (int k1 = 0; k1 <= outH + pad2_H - pool_sizeH; k1 = k1 + pool_stridH)
        {
            for (int k2 = 0; k2 <= outW + pad2_W - pool_sizeW; k2 = k2 + pool_stridW)
            {
                max = B[z][k1][k2];
                for (int i = k1; i < k1 + pool_sizeH; i++)
                {
                    for (int j = k2; j < k2 + pool_sizeW; j++)
                    {
                        if (B[z][i][j] > max)
                        {
                            max = B[z][i][j];
                        }
                    }
                }
                P_res[z][u1][u2] = max;

                u2++;
                if (u2 % in2_W == 0)
                {
                    u1++;
                    u2 = 0;
                }
            }
        }
    }
    /*    
    cout<<"test maxpool.."<<endl;
    cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
    cout<<"P_res[31][i][j]="<<endl;
    for (int i=0;i<in2_H;i++){
        for (int j=0;j<in2_W;j++){
            cout<<P_res[31][i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;*/
    return P_res;
}

int *Client_util::flatten(int ***P, int kerN, int in2_H, int in2_W)
{
    // cout<<"//Test flatten..."<<endl;
    //cout<<"C[3137]="<<endl;
    int num = kerN * in2_H * in2_W;
    int *C = new int[num];
    for (int i = 1; i < num; i++)
        C[i] = {0}; //flatten 7*7*64
    int u = 0;
    for (int z = 0; z < kerN; z++)
    {
        for (int i = 0; i < in2_H; i++)
        {
            for (int j = 0; j < in2_W; j++)
            {
                C[u] = P[z][i][j];
                u = u + 1;
            }
        }
    }
    //cout<<endl;
    return C;
}

void Client_util::Encode_Enc(int ***AA, int inD, int inH, int inW, string cipherfile, helib::EncryptedArray ea, 
ZZX *t, ZZX **A, ZZ_pX pm, ZZ_pX p_temp1, ZZ_pX p_temp2, ZZ_pX p_temp3, ZZ_pX p_temp4, ZZ q2d, enc_param ***usp)
{    
    timeval t1,t2;
    gettimeofday(&t1,NULL);
    util U;
    long slots = ea.size();
    std::vector<vector<long>> p1(inD, vector<long>(slots, 0));
    //一个输入通道加密为一条密文,所以p是两维
    for (int i = 0; i < inD; i++)
    {
        int z1 = 0;
        for (int j = 0; j < inH; j++)
        {
            for (int z = 0; z < inW; z++)
            {
                p1[i][z1] = AA[i][j][z];
                z1++;
            }
        }
    }
    ZZX encode1_poly[inD];
    for (int i = 0; i < inD; i++)
    {
        ea.encode(encode1_poly[i], p1[i]);
    }
   

    ZZX u[inD][k];
    ZZX v[inD];
 
    Enc enc;
    for (int z = 0; z < inD; z++)
    {
        enc.encrypt(encode1_poly[z], u[z], v[z], n, q, k, q2d, usp[0][0][z].r1, A, usp[0][0][z].e11, usp[0][0][z].e12, p_temp1, p_temp2, p_temp3, p_temp4, pm, t, dt, du, dv);
    }
    gettimeofday(&t2,NULL);
    double deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    cout << "编码+加密的计算总时间为:" << deltaT1/1000<<"ms"<< endl;

    ofstream outFile1;
    outFile1.open(cipherfile, ios::ate);
    for (int z = 0; z < inD; z++)
    {
        for (int j = 0; j < k; j++)
            outFile1 << u[z][j] << endl;
        outFile1 << v[z] << endl;
    }
    outFile1.close();
}

void Client_util::Encode_Enc(int *AA, int AA_len, string cipherfile, helib::EncryptedArray ea)
{
    util U;
    long slots = ea.size();
    std::vector<long> p(slots, 0);
    for (int i = 0; i < AA_len; i++)
    {
        p[i] = AA[i];
    }
    ZZX encode_poly;
    ea.encode(encode_poly, p);

    ZZX t[k];
    ifstream inFile1;
    inFile1.open("./file/key/pk-t.csv", ios::in);
    if (!inFile1)
    {
        cout << "Cannot open pk-t in conv1!\n";
        exit(0);
    }
    for (int i = 0; i < k; i++)
    {
        inFile1 >> t[i];
    }
    inFile1.close();
    for (int i = 0; i < k; i++)
    {
        U.decompress(t[i], n, q, dt);
    }
    //read A in pk-A
    ZZX **A = U.newZZXArray2(k, k);
    ifstream inFile2;
    inFile2.open("./file/key/pk-A.csv", ios::in);
    if (!inFile2)
    {
        cout << "Cannon open pk-A in conv1!\n";
        exit(0);
    }
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            inFile2 >> A[i][j];
        }
    }
    inFile2.close();
    ZZ_p::init(q);
    ZZ_pX ppm;
    SetCoeff(ppm, 0, 1);
    SetCoeff(ppm, n, 1);
    ZZ_pXModulus pm(ppm);
    ZZ_pX p_temp1, p_temp2, p_temp3, p_temp4;
    //Enc
    Enc enc;
    //q/2d
    ZZ two_d;
    RR two_d_rr;
    power2(two_d, dp); //two_d=2^dp
    two_d_rr = to_RR(two_d);
    RR q_rr;
    q_rr = to_RR(q);
    ZZ q2d;
    q2d = RoundToZZ(q_rr / two_d_rr); //q2d=round(q/2^dp)

    //u,vand c, k=2;
    ZZX u[k];
    ZZX v;
    enc_param usp;
    Client_util cu;
    usp = U.randnum_gen();

    enc.encrypt(encode_poly, u, v, n, q, k, q2d, usp.r1, A, usp.e11, usp.e12, p_temp1, p_temp2, p_temp3, p_temp4, pm, t, dt, du, dv);

    ofstream outFile1;
    outFile1.open(cipherfile, ios::ate);
    for (int j = 0; j < k; j++)
        outFile1 << u[j] << endl;
    outFile1 << v << endl;

    outFile1.close();
    
}

void Client_util::Genl_Dec(ZZX *h5, int kerN, string cipherfile, string skfile, int slots, long p)
{
    util U;
    ZZX m_dec1[kerN];
    ZZX u[kerN][k];
    ZZX v[kerN];
    ifstream inFile1;
    inFile1.open(cipherfile, ios::in);
    if (!inFile1)
    {
        cout << "Cannot open cipher1 file!\n";
        exit(0);
    }
    for (int z1 = 0; z1 < kerN; z1++)
    {
        inFile1 >> u[z1][0];
        inFile1 >> u[z1][1];
        inFile1 >> v[z1];
    }

    inFile1.close();

    //read sk.cvs
    ZZX s[k];
    ifstream inFile2;
    inFile2.open(skfile, ios::in);
    if (!inFile2)
    {
        cout << "Cannot open sk in dec one ploy!" << endl;
        exit(0);
    }
    for (int i = 0; i < k; i++)
    {
        inFile2 >> s[i];
    }
    inFile2.close();

 timeval t1,t2;
    gettimeofday(&t1,NULL);
    ZZX s0u10w1[kerN], s1u11w1[kerN], h6[kerN];
    ZZ_p::init(q2);
    ZZ_pX ppm;
    SetCoeff(ppm, 0, 1);
    SetCoeff(ppm, slots, 1);
    ZZ_pX pp_temp1, pp_temp2, pp_temp3, pp_temp4;

    ZZ two_d;
    RR two_d_rr;
    power2(two_d, dp);
    two_d_rr = to_RR(two_d);
    RR q_rr;
    q_rr = to_RR(q);
    ZZ q2d;
    q2d = RoundToZZ(q_rr / two_d_rr);

    RR two_dp_rr;
    power2(two_dp_rr, dp);
    RR q2_rr;
    ZZ two_d2;
    ZZ two_d3;
    two_d2 = two_d / 2;
    two_d3 = two_d2 * (-1);
    q2_rr = q_rr;

    for (int i = 0; i < kerN; i++)
    {
        U.Mult(s0u10w1[i], ppm, s[0], u[i][0], slots, q2, pp_temp1, pp_temp2, pp_temp3, pp_temp4);
        U.Mult(s1u11w1[i], ppm, s[1], u[i][1], slots, q2, pp_temp1, pp_temp2, pp_temp3, pp_temp4);
        h5[i] = v[i] - s0u10w1[i] - s1u11w1[i];
        for (int j = 0; j < n; j++)
        {
            RR coeff_rr;
            ZZ coeff_zz;
            GetCoeff(coeff_zz, h5[i], j);
            coeff_rr = to_RR(coeff_zz);
            coeff_zz = RoundToZZ(coeff_rr * two_dp_rr / q2_rr) % two_d;
            if (coeff_zz >= two_d2)
            {
                coeff_zz = coeff_zz - two_d;
            }
            else if (coeff_zz <= two_d3)
            {
                coeff_zz = coeff_zz + two_d;
            }
            SetCoeff(h5[i], j, coeff_zz);
        }
        U.Module(h5[i], n, conv<ZZ>(p));
    }
        gettimeofday(&t2,NULL);
    double deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    cout << "一般的解密的计算时间为:" << deltaT1/1000<<"ms"<< endl;

}

int ***Client_util::Decode_Relu(ZZX *h, int kerN, int outH, int outW, helib::EncryptedArray ea, long p)
{
    timeval t1,t2;
    gettimeofday(&t1,NULL);
    util U;
    int slots = ea.size();
    std::vector<vector<long>> pp(kerN, vector<long>(slots, 0));
    for (int i = 0; i < kerN; i++)
        ea.decode(pp[i], h[i]);

    for (int i = 0; i < kerN; i++)
    {
        for (int j = 0; j < slots; j++)
        {
            pp[i][j] = pp[i][j] % p;
            if (pp[i][j] >= p / 2)
            {
                pp[i][j] = pp[i][j] - p;
            }
            else if (pp[i][j] <= -p / 2)
            {
                pp[i][j] = pp[i][j] + p;
            }
        }
    }
    int ***mat = U.newArray3(kerN, outH, outW);
    for (int i = 0; i < kerN; i++)
    {
        for (int j = 0; j < outH; j++)
        {
            for (int z = 0; z < outW; z++)
            {
                mat[i][j][z] = pp[i][j * outW + z];
                if(mat[i][j][z] < 0){
                    mat[i][j][z] = 0;
                }
            }
        }
    }
    gettimeofday(&t2,NULL);
    double deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    cout << "模块E中解码和relu的计算总时间为:" << deltaT1/1000<<"ms"<< endl;

    return mat;
}

int *Client_util::Decode_Relu(ZZX h, int C_len, helib::EncryptedArray ea, long p)
{
    util U;
    int slots = ea.size();
    std::vector<long> pp(slots, 0);
    ea.decode(pp, h);
    for (int i = 0; i < slots; i++)
    {
        pp[i] = pp[i] % p;
        if ((pp[i] <= p / 2) && (pp[i] >= -p / 2))
        {
        }
        else if (pp[i] > p / 2)
        {
            pp[i] = pp[i] - p;
        }
        else if (pp[i] < -p / 2)
        {
            pp[i] = pp[i] + p;
        }
    }
    int *mat;
    for (int i = 0; i < C_len; i++)
    {
        mat[i] = pp[i];
        if (mat[i] < 0)
        {
            mat[i] = 0;
        }
    }
    return mat;
}

int ***Client_util::Decode(ZZX *h, int kerN, int outH, int outW, helib::EncryptedArray ea, long p)
{
    util U;
    int slots = ea.size();
    std::vector<vector<long>> pp(kerN, vector<long>(slots, 0));
    for (int i = 0; i < kerN; i++)
        ea.decode(pp[i], h[i]);
    for (int i = 0; i < kerN; i++)
    {
        for (int j = 0; j < slots; j++)
        {
            pp[i][j] = pp[i][j] % p;
            if ((pp[i][j] <= (p - 1) / 2) && (pp[i][j] >= -(p - 1) / 2))
            {
            }
            else if (pp[i][j] > 0)
            {
                pp[i][j] = pp[i][j] - p;
            }
            else if (pp[i][j] < 0)
            {
                pp[i][j] = pp[i][j] + p;
            }
        }
    }

    int ***mat = U.newArray3(kerN, outH, outW);
    for (int i = 0; i < kerN; i++)
    {
        for (int j = 0; j < outH; j++)
        {
            for (int z = 0; z < outW; z++)
            {
                mat[i][j][z] = pp[i][j * outW + z];
            }
        }
    }
    return mat;
}

int *Client_util::Decode(ZZX h, int C_len, helib::EncryptedArray ea, long p)
{

    int slots = ea.size();
    std::vector<long> pp(slots, 0);
    ea.decode(pp, h);
    for (int i = 0; i < slots; i++)
    {
        pp[i] = pp[i] % p;
        if ((pp[i] <= (p - 1) / 2) && (pp[i] >= -(p - 1) / 2))
        {
        }
        else if (pp[i] > 0)
        {
            pp[i] = pp[i] - p;
        }
        else if (pp[i] < 0)
        {
            pp[i] = pp[i] + p;
        }
    }

    int *mat;
    for (int i = 0; i < C_len; i++)
    {
        mat[i] = pp[i];
    }
    return mat;
}

int *Client_util::dense_input(int *C, int lenth_C)
{
    int *CC = new int[lenth_C + 1];
    CC[0] = 1;
    for (int i = 1; i < lenth_C + 1; i++)
        CC[i] = C[i - 1];
    // cout<<"The X in fc ="<<endl;
    //  for(int i=0;i<lenth_C+1;i++)
    //     cout<<CC[i]<<" ";
    //     cout<<endl;
    return CC;
}
