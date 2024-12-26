#include <client_cipher.h>
#include <sys/time.h>


using namespace std;
using namespace NTL;

Client_cipher::Client_cipher(){

}


int ***Client_cipher::DecCNN2D_Conv1(int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W, int pool_stridH, int pool_stridW,  string evalfile, string skfile, string padding){
  
    Client_util cu;
    util U;
    int outW = ceil((inW*1.0)/strid_W);   
    int outH = ceil((inH*1.0)/strid_H);   
    int datrows = kerH * kerW;
    int datcols = outH * outW; 
    int real_nrow = datrows + 1;            
    int max_packed_num = n / (real_nrow);
    int max_pool_num = max_packed_num / (pool_stridH * pool_stridW);
    int packed_num = max_pool_num * pool_stridH * pool_stridH; 
    int packed_col = ceil(double(1.0 * datcols/packed_num));
    
    ZZX **h1 = U.newZZXArray2(kerN, packed_col);
    ZZX **h2 = U.newZZXArray2(kerN, packed_col);
    ZZX **h3 = U.newZZXArray2(kerN, packed_col);
    ifstream inFile1;
    inFile1.open(evalfile, ios::in);
    if (!inFile1){
        cout << "Cannot open homomorphic eval result file!" << endl;
        exit(0);
    }

    for (int z = 0; z < kerN; z++){
        for(int i = 0; i < packed_col; i++){         
            inFile1 >> h1[z][i];            
            inFile1 >> h2[z][i];
            inFile1 >> h3[z][i];                       
        } 
    }    
    inFile1.close();

    timeval t1, t2;
    gettimeofday(&t1, NULL);
    ZZX WX1;                                        
    ZZ *p1 = new ZZ[packed_num]; 
    int ** inner_res = U.newArray2(kerN, packed_col * packed_num);  
   
    Enc enc;
    for(int z = 0; z < kerN; z++){
        for(int i = 0; i < packed_col; i++){ 
            WX1 = cu.dec_onepackedPloy(skfile, h1[z][i], h2[z][i], h3[z][i]);                  
            p1 = cu.GetpakedPCoeff_bias(WX1, kerH, kerW, packed_num);
            for(int j=0; j < packed_num; j++){
                inner_res[z][i * packed_num + j] = conv<int>(p1[j]);                     
            }
        }
    }
    gettimeofday(&t2, NULL);
    double deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    cout << "CONV dec time is " << (double)deltaT / 1000.0 << "ms" << endl;
 
    int poolCol = packed_num * packed_col/(pool_stridW * pool_stridH);//266
    int** pool_res1 = U.newArray2(kerN, poolCol);
    int stepCol = pool_stridW * pool_stridH;
    timeval t5, t6;
    gettimeofday(&t5, NULL);
    int max = 0;
    for(int z = 0; z < kerN; z++){      
        for(int i = 0; i < poolCol; i++){
            max = inner_res[z][stepCol * i];
            for(int j = 0; j < stepCol; j++){
                if(inner_res[z][stepCol * i + j] > max)
                   max = inner_res[z][stepCol * i + j];
            }
            pool_res1[z][i] =  max * (double(1.0 / (w_precision)));
        }
    }

    int*** pool_res2 = U.newArray3(kerN, outH / pool_stridH, outW / pool_stridW);
    for(int z = 0; z < kerN; z++){
        for(int i = 0; i < outH / pool_stridH; i++){
            for(int j = 0; j < outW / pool_stridW; j++){
                pool_res2[z][i][j] = pool_res1[z][i*outW/pool_stridW+j];
            }
        }
    }  
    cout<<"end maxpool"<<endl;
    gettimeofday(&t6, NULL);
    double deltaT3 = (t6.tv_sec - t5.tv_sec) * 1000000 + t6.tv_usec - t5.tv_usec;
    cout << "conv1 maxpool time is " << (double)deltaT3 / 1000.0 << "ms" << endl;
 

    // gettimeofday(&t6, NULL);
    // double deltaT3 = (t6.tv_sec - t5.tv_sec) * 1000000 + t6.tv_usec - t5.tv_usec;
    // cout << "conv1 maxpool time is " << (double)deltaT3 / 1000.0 << "ms" << endl;
    U.deleteArray2(inner_res,kerN, packed_num*packed_col);  
    U.deleteZZXArray2(h1, kerN, packed_col);
    U.deleteZZXArray2(h2, kerN, packed_col);
    U.deleteZZXArray2(h3, kerN, packed_col);
    delete []p1;
    return pool_res2;
}

int ***Client_cipher::DecCNN3D_Conv2(int packedRowCount, int packedColCount, int packedRowUnit, int packedColUnit,  string evalfile, string skfile, int kerN, int kerH, int kerW, int inH, int inW, int strid_H,int strid_W,string padding)
{   util U;
    int outH = ceil((inH*1.0) / strid_H);
    int outW = ceil((inW*1.0) / strid_W);
  
    Client_util cu;
    clock_t startTime11;
    clock_t endTime11;
    startTime11 = clock();
    ZZX ***h1 = U.newZZXArray3(kerN, packedRowCount,  packedColCount);//16*4*4
    ZZX ***h2 = U.newZZXArray3(kerN, packedRowCount,  packedColCount);
    ZZX ***h3 = U.newZZXArray3(kerN, packedRowCount,  packedColCount);
    
    ZZX ***CsumMat1=U.newZZXArray3(kerN, packedRowCount,  packedColCount);
    ZZX ***CsumMat2=U.newZZXArray3(kerN, packedRowCount,  packedColCount);
    ZZX ***CsumMat3=U.newZZXArray3(kerN, packedRowCount,  packedColCount);

    ifstream inFile1;
    inFile1.open(evalfile, ios::in);
    if (!inFile1){
        cout << "Cannot open homomorphic eval result file!" << endl;
        exit(0);
    }

    for (int z = 0; z < kerN; z++){
        for(int i=0; i<packedRowCount; i++){
            for (int j=0; j <packedColCount; j++){
                inFile1 >> h1[z][i][j];
                inFile1 >> h2[z][i][j];
                inFile1 >> h3[z][i][j];
            }
        } 
    }
     for (int z = 0; z < kerN; z++){
        for(int i=0; i<packedRowCount; i++){
            for (int j = 0; j<packedColCount; j++){
                inFile1 >> CsumMat1[z][i][j];
                inFile1 >> CsumMat2[z][i][j];
                inFile1 >> CsumMat3[z][i][j];
            }
        } 
    }
    inFile1.close();
    // endTime11 = clock();
    // cout << "Homomorphic res read time is " << (double)(endTime11 - startTime11) / CLOCKS_PER_SEC << "s" << endl;
  
    timeval t1, t2;
    gettimeofday(&t1, NULL);
    // 用户解密第一通道的结果
    ZZX WX1;                                        //第一通道所有多项式相乘计算结果
    ZZ *p1 = new ZZ[packedRowUnit * packedColUnit]; //packed_num 4
    int*** inner_result = U.newArray3(kerN,outH+1,outW+1);  //16*8*8
 
    // cout<<"=====test decrypt result polys[0][4][2]========"<<endl;
    // WX1 = cu.dec_onepackedPloy(skfile, h1[0][4][2], h2[0][4][2], h3[0][4][2]);
    // p1 = cu.GetpakedPCoeff_bias(WX1, kerH, kerW, packedRowUnit * packedColUnit);
    // for(int i=0;i<packedRowUnit * packedColUnit;i++)
    //     cout<<p1[i]<<" ";
    // cout<<endl;
    
    for(int z=0; z<kerN; z++){//16
        for(int i=0; i<packedRowCount; i++){ //4
            for(int j=0; j<packedColCount; j++){ //4
                WX1 = cu.dec_onepackedPloy(skfile, h1[z][i][j], h2[z][i][j], h3[z][i][j]);
                p1 = cu.GetpakedPCoeff_bias(WX1, kerH, kerW, packedRowUnit * packedColUnit);//每个密文多项式结出2*2=4个内积结果
                int cid1 = 0;
                int cid2 = 0;
                for(int u=0; u<packedRowUnit*packedColUnit; u++){
                    // cout<<p1[u]<<" ";
                    int value = conv<int>(p1[u]);
                    //  double yValue = conv<double>(p1[u]);
                    // int value = cu.relu(yValue * (double(1.0 / w_precision)));
                   
                    if(u<2){
                        inner_result[z][i*packedRowUnit][j*packedColUnit+cid1] = value;
                        cid1++;
                    }
                    else{
                        inner_result[z][i*packedRowUnit+1][j*packedColUnit+cid2] = value;
                        cid2++;
                    }
                }
                // cout<<endl;
            }//end for j
            
        }//end for i       
    }//end for z  
    // cout<<endl; 
    // cout<<"test decryption result of inner product in conv2 "<<endl;
    // for(int i=0;i<outH+1;i++){
    //     for(int j=0;j<outW+1;j++){
    //         cout<<inner_result[0][i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;

     //解密其余通道的密文
    ZZX WX2;                                        
    ZZ *p2 = new ZZ[packedRowUnit * packedColUnit]; 
    int*** inner_result2 = U.newArray3(kerN,outH+1, outW+1);
    for(int z=0; z<kerN; z++){
        for(int i=0; i<packedRowCount; i++){ //4
            for(int j=0; j<packedColCount; j++){ //4
                WX2 = cu.dec_onepackedPloy(skfile, CsumMat1[z][i][j],  CsumMat2[z][i][j],  CsumMat3[z][i][j]);
                p2 = cu.GetpakedPCoeff_Nobias(WX2, kerH, kerW, packedRowUnit * packedColUnit);//每个密文多项式结出2*2=4个内积结果
                int cid1 = 0;
                int cid2 = 0;
                for(int u=0; u<packedRowUnit*packedColUnit; u++){
                    int value = conv<int>(p2[u]);
                    // int value = cu.relu(yValue * (double(1.0 / w_precision)));
                    if(u<2){
                        inner_result2[z][i*packedRowUnit][j*packedColUnit+cid1] = value;
                        cid1++;
                    }
                    else{
                        inner_result2[z][i*packedRowUnit+1][j*packedColUnit+cid2] = value;
                        cid2++;
                    }
                }
            }//end for j
        }//end for i
    }//end for z   

    gettimeofday(&t2, NULL);
    double deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
   
    cout << "CONV2 dec time is " << (double)deltaT / 1000.0 << "ms" << endl;
    delete[] p1,p2;

    // 每个输出通道内，两个输入通道的对应值相加。并求relu
    timeval t3, t4;
    gettimeofday(&t3, NULL);
    int dec_Mat[kerN][outH][outW];
    int ***relu_result=U.newArray3(kerN,outH,outW);
        for (int i = 0; i < kerN; i++)
        {
            for (int j = 0; j < outH; j++)
            {
                for(int z=0;z< outW;z++)
                {
                    dec_Mat[i][j][z] = inner_result[i][j][z] + inner_result2[i][j][z];
                    relu_result[i][j][z] = cu.relu(dec_Mat[i][j][z] * (double(1.0 / w_precision)));
                }
            }
        }
    gettimeofday(&t4, NULL);
    double deltaT2 = (t4.tv_sec - t3.tv_sec) * 1000000 + t4.tv_usec - t3.tv_usec;
    cout << "conv2 homomorphic add and relu time is " << (double)deltaT2 / 1000.0 << "ms" << endl;

    U.deleteArray3(inner_result, kerN, outH, outW);
    U.deleteZZXArray3(h1, kerN, packedRowCount, packedColCount);
    U.deleteZZXArray3(h2, kerN, packedRowCount, packedColCount);
    U.deleteZZXArray3(h3, kerN, packedRowCount, packedColCount); 
   
    U.deleteZZXArray3(CsumMat1, kerN, packedRowCount, packedColCount);
    U.deleteZZXArray3(CsumMat2, kerN, packedRowCount, packedColCount);
    U.deleteZZXArray3(CsumMat3, kerN, packedRowCount, packedColCount);

    return relu_result;
}




void Client_cipher::input_cipherGen_Conv1(int ***AA, int inD, int inH, int inW, int kerN, int kerH, int kerW, int strid_H, int strid_W, const char *cipherfile, string padding)
{
  util U;
    Client_util cu;
    int outH = ceil(inH / strid_H);
    int outW = ceil(inW / strid_W);
    int ***C = U.newArray3(inD, kerH*kerW, outH*outW);
    U.initArray3(C,inD, kerH*kerW, outH*outW,0);
    // padding to 32*32


    int datrows = kerH * kerW;
    int datcols = outH * outW;             //输入矩阵dat 的行:25和列:784
    int pcol = outW; 
    int part = ceil(double(datcols) / double(pcol)); //pcol=28. part是将dat的列要截断为28个部分
    int prow = part * (datrows + 1);                              //将dat 每28列一截断生成的矩阵Pdat的行和列1.e.28*26
       
    int eleCeil = datrows + 1;       // 1:bias  total:26
    int packedNum_row = 2 * eleCeil; //26*2
    int packedNum_col = 4;
    int packedColCount = ceil(pcol / (packedNum_col * 1.0));
    int packedRowCount = ceil(prow / (packedNum_row * 1.0));


        
    ZZX t[k];
    for (int i = 0; i < k; i++){
        t[i] = ZZX();
    }
    
    ifstream inFile1;
    inFile1.open("./file/key/pk-t.csv", ios::in);
    if (!inFile1){
        cout << "Cannot open pk-t in conv1!\n";
        exit(1);
    }
    for (int i = 0; i < k; i++){
        inFile1 >> t[i];
    }
    inFile1.close();   

    ZZX **A = U.newZZXArray2(k, k);
    ifstream inFile2;
    inFile2.open("./file/key/pk-A.csv", ios::in);
    if (!inFile2) {
        cout << "Cannon open pk-A in conv1!\n";
        exit(0);
    }
    for (int i = 0; i < k; i++){
        for (int j = 0; j < k; j++){
            inFile2 >> A[i][j];
        }
    }
    inFile2.close();


    ZZ_p::init(q);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,n,1);
    ZZ_pXModulus pm(ppm);
    ZZ_pX p_temp1,p_temp2,p_temp3,p_temp4;
    for (int i = 0; i < k; i++){
        U.decompress(t[i], n, q, dt);
    }
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
    enc_param usp[packedRowCount][packedColCount];
    for(int i=0;i<packedRowCount;i++)
        for(int j=0;j<packedColCount;j++)
            usp[i][j] = U.randnum_gen(); 

    
    timeval p1, p2;
    gettimeofday(&p1, NULL);
    cu.slideToMat_SAME(C,AA,inD,inH,inW,kerH,kerW,strid_H,strid_W);
    ZZX ***polys = U.newZZXArray3(inD, packedRowCount, packedColCount);//14*7
     for (int i = 0; i < inD; i++){
        cu.matPacking_by_poolsize(polys[i], C[i], datrows, datcols, pcol);
    }
    ZZX ***u2 = U.newZZXArray3(packedRowCount, packedColCount, k);
    ZZX **v2 = U.newZZXArray2(packedRowCount, packedColCount);
    
    for (int i = 0; i < packedRowCount; i++){
        for (int j = 0; j < packedColCount; j++){       
            enc.encrypt(polys[0][i][j], u2[i][j], v2[i][j], n, q, k, q2d, usp[i][j].r1, A, usp[i][j].e11, usp[i][j].e12, p_temp1, p_temp2, p_temp3, p_temp4, pm, t, dt, du, dv);          
        }
    }
        gettimeofday(&p2, NULL);
    double deltaT1 = (p2.tv_sec - p1.tv_sec) * 1000000 + p2.tv_usec - p1.tv_usec;
    cout << "con1编码和加密总时间为:" << deltaT1/1000<<"ms"<< endl;


    fstream fw(cipherfile);
    if (fw.is_open()){
        for (int i = 0; i < packedRowCount; i++){
            for (int j = 0; j < packedColCount; j++){
                for (int z = 0; z < k; z++){
                    fw << u2[i][j][z] << "\n";
                }
                    fw << v2[i][j] << "\n";
            }
        }
    }
    fw.close();

    U.deleteZZXArray3(u2, packedRowCount, packedColCount, k);
    U.deleteZZXArray2(v2, packedRowCount, packedColCount);
    U.deleteZZXArray3(polys, inD, packedRowCount, packedColCount);    
    U.deleteArray3(C, inD, kerH * kerW, outH * outW);

}


void Client_cipher::input_cipherGen(int***AA,int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W,const char *cipherfile,string padding)
{  
    util U;
    Client_util cu;
    if (padding=="SAME"){  
        int outH=ceil(inH*1.0/strid_H);
        int outW=ceil(inW*1.0/strid_W);

        int packed_num1 = n/(kerH*kerW+1);
        int packed_col1 =ceil(double(outH*outW)/packed_num1);
        int packed_num2 = n/(kerH*kerW);
        int packed_col2 =ceil(double(outH*outW)/packed_num2);

        int*** C=U.newArray3(inD,kerH*kerW,outH*outW);      
        cu.slideToMat_SAME(C,AA,inD,inH,inW,kerH,kerW,strid_H,strid_W);
        cout<<inD<<";"<<kerH*kerW<<";"<<outH*outW<<endl;
    
    int*** packedC; 
    packedC=U.newArray3(inD,n,packed_col1);
    U.initArray3(packedC,inD,n,packed_col1,0);
     cout<<inD<<";"<<n<<";"<<packed_col1<<endl;
    //第一通道加上一行1,并且倒置填充矩阵得到packedC
    cu.matPacking(C[0],kerH*kerW,outH*outW,packedC[0],n,packed_col1);
    //其余通道不加一行1,倒置填充矩阵得到packedC
    for(int z=1;z<inD;z++)
        cu.matPackingNoOne(C[z],kerH*kerW,outH*outW,packedC[z],n,packed_col2);
    //将packedC装配为多项式M
    ZZX M[inD][packed_col1];//packed_col ploy and each poly's degree is 256	
    for(int u=0;u<inD;u++){  
        for(int j=0;j<packed_col1;j++){
            M[u][j].SetLength(n);
            for(int i=0;i<n;i++){
                SetCoeff(M[u][j],i,packedC[u][i][j]);
                }
            }
    }
    //释放内存
    U.deleteArray3(packedC,inD,n,packed_col1);
    U.deleteArray3(C,inD,kerH*kerW,outH*outW);  
  

    // 加密M 
    // read t in pk.csv 
	ZZX t[k];
    ifstream inFile1;
    inFile1.open("./file/key/pk-t.csv",ios::in); 
    if(!inFile1){
        cout<<"Cannot open pk-t in conv1!\n";
        exit(0);
    }   
    for(int i=0;i<k;i++){
        inFile1>>t[i];       
    }    
    inFile1.close();	
 
	for(int i=0;i<k;i++){
        U.decompress(t[i],n,q,dt);
    }

	//read A in pk-A
    ZZX **A=U.newZZXArray2(k,k);
    ifstream inFile2;
    inFile2.open("./file/key/pk-A.csv",ios::in);
    if(!inFile2){
        cout<<"Cannon open pk-A in conv1!\n";
        exit(0);
    }
    for(int i=0;i<k;i++){
        for(int j=0;j<k;j++){
            inFile2>>A[i][j];     
        }
	}    
    inFile2.close();
	
    ZZ_p::init(q);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,n,1);
    ZZ_pXModulus pm(ppm);
    ZZ_pX p_temp1,p_temp2,p_temp3,p_temp4;
    //Enc
    Enc enc;
    //q/2d
    ZZ two_d;
    RR two_d_rr;
    power2(two_d,dp);//two_d=2^dp
    two_d_rr=to_RR(two_d);
    RR q_rr;
    q_rr = to_RR(q);
    ZZ q2d;
    q2d=RoundToZZ(q_rr/two_d_rr);//q2d=round(q/2^dp)
	
    //u,vand c, k=2;
    ZZX ***u2=U.newZZXArray3(inD,packed_col1,k);
    ZZX **v2=U.newZZXArray2(inD,packed_col1);
    //ZZX v2[inD][packed_col];

    
   // client_enc_param usp;
    //usp=cu.randnum_gen(); 
    
    enc_param usp;
    clock_t startTime;
    clock_t endTime;
    startTime=clock();
   
    for(int z=0;z<inD;z++){
        for(int i=0;i<packed_col1;i++){           
            usp=U.randnum_gen();         
            enc.encrypt(M[z][i],u2[z][i],v2[z][i],n,q,k,q2d,usp.r1,A,usp.e11,usp.e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv);
        }
    }
    // endTime=clock();
    // cout<<"CONV enc input time is "<<(double)(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;
    
    //  store inputcipher file
    ofstream outFile1;
    outFile1.open(cipherfile,ios::ate);
    for(int z=0;z<inD;z++){
        for(int i=0;i<packed_col1;i++){
            for (int j=0;j<k;j++){
                outFile1<<u2[z][i][j]<<endl;
            }     
            outFile1<<v2[z][i]<<endl; 
        }  
    }
    outFile1.close(); 

    // long cu2[inD][packed_col1][k][n];
    // long cv2[inD][packed_col1][n];
    // for(int z=0;z<inD;z++){
    //     for(int i=0;i<packed_col1;i++){
    //         for (int j=0;j<k;j++){
    //             for(int u=0;u<n;u++){
    //                 cu2[z][i][j][u]=conv<long>(coeff(u2[z][i][j],u));
    //             }                    
    //         }
    //         for(int v=0;v<n;v++){
    //             cv2[z][i][v]=conv<long>(coeff(v2[z][i],v));
    //         }
    //     }
    // }

    // FILE *pFile;
    // pFile=fopen(cipherfile,"w");    
    // fwrite(cu2,sizeof(cu2[0][0][0][0]),inD*packed_col1*k*n,pFile);
    // fwrite(cv2,sizeof(cv2[0][0][0]),inD*packed_col1*n,pFile);
    // fclose(pFile);

    U.deleteZZXArray3(u2,inD,packed_col1,k); 
    U.deleteZZXArray2(A,k,k); 
    U.deleteZZXArray2(v2,inD,packed_col1);     
    }


    else if(padding=="VALID"){
               
        int outH=ceil(float(inH-kerH+1)/float(strid_H));
        int outW=ceil(float(inW-kerW+1)/float(strid_W));
        cout<<"outH="<<outH<<endl;
        int packed_num1 = n/(kerH*kerW+1);
        int packed_col1 =ceil(double(outH*outW)/packed_num1);
        int packed_num2 = n/(kerH*kerW);
        int packed_col2 =ceil(double(outH*outW)/packed_num2);
        cout<<"********************************"<<endl;
        cout<<"packed_col1="<<packed_col1<<endl;
        cout<<"packed_col2="<<packed_col2<<endl;
        cout<<"********************************"<<endl;

        int*** C=U.newArray3(inD,kerH*kerW,outH*outW);       
        cu.slideToMat_VALID(C,AA,inD,inH,inW,kerH,kerW,strid_H,strid_W);        
        ZZX **N=U.newZZXArray2(inD,outH*outW);
    
        //输入图像的第一个通道要装配上1,其余通道则不需要。
        for(int j=0;j<outH*outW;j++){
        N[0][j].SetLength(kerH*kerW+1);
            SetCoeff(N[0][j],0,1);
            for(int i=0;i<kerH*kerW;i++){
                SetCoeff(N[0][j],i+1,C[0][i][j]);
            }
        }
        
        //输入图像的第一个通道要装配上1,其余通道则不需要。
        for(int z=1;z<inD;z++){
            for(int j=0;j<outH*outW;j++){
                N[z][j].SetLength(kerH*kerW);
                    for(int i=0;i<kerH*kerW;i++){
                        SetCoeff(N[z][j],i,C[z][i][j]);
                    }
                }
        }
        // //存储多项式
        // ofstream outFile2;
        // outFile2.open(plainfile,ios::ate);
        // for(int z=0;z<inD;z++)
        //     for(int i=0;i<outH*outW;i++){
        //         outFile2<<N[z][i]<<endl;
        //     }
        // outFile2.close();
        U.deleteZZXArray2(N,inD,outH*outW);
        
        int*** packedC; 
        packedC=U.newArray3(inD,n,packed_col1);
        U.initArray3(packedC,inD,n,packed_col1,0);
        //第一通道加上一行1,并且倒置填充矩阵packedC
        cu.matPacking(C[0],kerH*kerW,outH*outW,packedC[0],n,packed_col1);
        //其余通道不加一行1,倒置填充矩阵packedC
        for(int z=1;z<inD;z++)
            cu.matPackingNoOne(C[z],kerH*kerW,outH*outW,packedC[z],n,packed_col2);
        //将packedC装配为多项式M
        ZZX M[inD][packed_col1];//packed_col ploy and each poly's degree is 256	
        for(int u=0;u<inD;u++){  
            for(int j=0;j<packed_col1;j++){
                M[u][j].SetLength(n);
                for(int i=0;i<n;i++){
                    SetCoeff(M[u][j],i,packedC[u][i][j]);
                    }
                }
        }
        //释放内存
        U.deleteArray3(packedC,inD,n,packed_col1);
        U.deleteArray3(C,inD,kerH*kerW,outH*outW);  
    

        // 加密M 
        // read t in pk.csv 
        ZZX t[k];
        ifstream inFile1;
        inFile1.open("./file/key/pk-t.csv",ios::in); 
        if(!inFile1){
            cout<<"Cannot open pk-t in conv1!\n";
            exit(0);
        }   
        for(int i=0;i<k;i++){
            inFile1>>t[i];       
        }    
        inFile1.close();	
    
        for(int i=0;i<k;i++){
            U.decompress(t[i],n,q,dt);
        }

        //read A in pk-A
        ZZX **A=U.newZZXArray2(k,k);
        ifstream inFile2;
        inFile2.open("./file/key/pk-A.csv",ios::in);
        if(!inFile2){
            cout<<"Cannon open pk-A in conv1!\n";
            exit(0);
        }
        for(int i=0;i<k;i++){
            for(int j=0;j<k;j++){
                inFile2>>A[i][j];     
            }
        }    
        inFile2.close();
        
        ZZ_p::init(q);
        ZZ_pX ppm;
        SetCoeff(ppm,0,1);
        SetCoeff(ppm,n,1);
        ZZ_pXModulus pm(ppm);
        ZZ_pX p_temp1,p_temp2,p_temp3,p_temp4;
        //Enc
        Enc enc;
        //q/2d
        ZZ two_d;
        RR two_d_rr;
        power2(two_d,dp);//two_d=2^dp
        two_d_rr=to_RR(two_d);
        RR q_rr;
        q_rr = to_RR(q);
        ZZ q2d;
        q2d=RoundToZZ(q_rr/two_d_rr);//q2d=round(q/2^dp)
        
        //u,vand c, k=2;
        ZZX ***u2=U.newZZXArray3(inD,packed_col1,k);
        ZZX **v2=U.newZZXArray2(inD,packed_col1);
        //ZZX v2[inD][packed_col];

        
        // client_enc_param usp;
        //usp=cu.randnum_gen(); 
        
        enc_param usp;
        clock_t startTime;
        clock_t endTime;
        startTime=clock();
    
        for(int z=0;z<inD;z++){
            for(int i=0;i<packed_col1;i++){           
                usp=U.randnum_gen();  
                enc.encrypt(M[z][i],u2[z][i],v2[z][i],n,q,k,q2d,usp.r1,A,usp.e11,usp.e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv);
            }
        }
        endTime=clock();
        cout<<"CONV enc input time is "<<(double)(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;
        //store inputcipher file
        // ofstream outFile1;
        // outFile1.open(cipherfile,ios::ate);
        // for(int z=0;z<inD;z++){
        //     for(int i=0;i<packed_col1;i++){
        //         for (int j=0;j<k;j++){
        //             outFile1<<u2[z][i][j]<<endl;
        //         }     
        //         outFile1<<v2[z][i]<<endl; 
        //     }  
        // }
        // outFile1.close(); 

        long cu2[inD][packed_col1][k][n];
        long cv2[inD][packed_col1][n];
        for(int z=0;z<inD;z++){
            for(int i=0;i<packed_col1;i++){
                for (int j=0;j<k;j++){
                    for(int u=0;u<n;u++){
                        cu2[z][i][j][u]=conv<long>(coeff(u2[z][i][j],u));
                    }                    
                }
                for(int v=0;v<n;v++){
                    cv2[z][i][v]=conv<long>(coeff(v2[z][i],v));
                }
            }
        }
        FILE *pFile;
        pFile=fopen(cipherfile,"wb");    
        fwrite(cu2,sizeof(cu2[0][0][0][0]),inD*packed_col1*k*n,pFile);
        fwrite(cv2,sizeof(cv2[0][0][0]),inD*packed_col1*n,pFile);
        fclose(pFile);
        
        U.deleteZZXArray3(u2,inD,packed_col1,k); 
        U.deleteZZXArray2(A,k,k); 
        U.deleteZZXArray2(v2,inD,packed_col1);     
        
        }  
   
}

void Client_cipher::input_cipherGen_Conv2(int***AA,int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W,int pool_strid_H, int pool_strid_W, const char *cipherfile,string padding,ZZX *t, ZZX** A,  
ZZ_pX pm,  ZZ_pX p_temp1,  ZZ_pX p_temp2,  ZZ_pX p_temp3,  ZZ_pX p_temp4, ZZ q2d,enc_param ***usp){  
    util U;
    Client_util cu;
    if(padding=="SAME")
    {  
        int outH = ceil((inH * 1.0)/strid_H);
        int outW = ceil((inW * 1.0)/strid_W);      
        int datrows = kerH * kerW;
        int datcols = outH * outW; 
        int real_nrow = datrows + 1;            
        int max_packed_num = n / (real_nrow);
        int max_pool_num = max_packed_num / (pool_strid_H * pool_strid_W);
        int packed_num = max_pool_num * pool_strid_H * pool_strid_H; 
        int packed_col = ceil(double(1.0 * datcols / packed_num));
        // cout << "max_pool-num=" << max_pool_num << endl;
        // cout << "packed_num=" << packed_num << endl;
        // cout << "paked_col=" << packed_col << endl;
        
        int*** C = U.newArray3(inD, kerH * kerW, outH * outW);         
   
        //right 21.11.25
        cu.slideToMat_SAME(C, AA, inD, inH, inW, kerH, kerW, strid_H, strid_W);
    
       //所有通道按照方法2装配为多项式 
        ZZX **poly = new ZZX*[inD];
        for(int i = 0; i < inD; i++)
            poly[i] = new ZZX[packed_col];
        int pcol = outW;  
        
        for (int i = 0; i < inD; i++){
            //outH*ouW列每个列都添加了1，再倒排,再装配为多项式
            cu.packingbyPoolWin(poly[i], C[i], datrows, datcols, pcol, pool_strid_H, pool_strid_W, max_pool_num);          
        } 
        U.deleteArray3(C, inD, kerH*kerW, outH*outW);    
       

        ZZX ***u2 = U.newZZXArray3(inD, packed_col, k);
        ZZX **v2 = U.newZZXArray2(inD, packed_col);

        Enc enc;      
        for(int z = 0; z < inD; z++){
            for(int i = 0; i < packed_col; i++){                        
                enc.encrypt(poly[z][i],u2[z][i],v2[z][i],n,q,k,q2d,usp[0][z][i].r1, A, usp[0][z][i].e11, usp[0][z][i].e12, 
                p_temp1, p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv);
            }
        } 
     
        //  store inputcipher file
        ofstream outFile1;
        outFile1.open(cipherfile, ios::ate);
        for(int z = 0; z < inD; z++){
            for(int i = 0; i < packed_col; i++){                
                for(int j = 0; j < k; j++){
                    outFile1 << u2[z][i][j] << endl;
                }
                outFile1 << v2[z][i] << endl; 
            }
        }        
        outFile1.close();      
        U.deleteZZXArray3(u2, inD, packed_col, k);   
        U.deleteZZXArray2(v2, inD, packed_col);  
    }    
}



void Client_cipher::dense_input_cipherGen(int* C, int lenth_C,string fc_inputcipherfile, ZZX *t, ZZX** A,  ZZ_pX pm,  ZZ_pX p_temp1,  ZZ_pX p_temp2,  ZZ_pX p_temp3,  ZZ_pX p_temp4, ZZ q2d, enc_param *usp)
{  
    util U;
    //timeval p1, p2;
    // gettimeofday(&p1, NULL);   
    int packednum=ceil(1.0*lenth_C/(n-1));       
    // int slots=ea.size();
    // ZZX encode_M[packednum];//for f1:13 ploy and each poly's degree n is 256
    // std::vector <vector<long>>p(packednum,vector<long>(slots,0));

	// for(int i=0;i<packednum;i++){
    //     p[i][0]=1;
    //     for(int j=0;j<slots-1;j++){
    //         p[i][j+1]=C[i*(slots-1)+j];
    //     }       
    //     // ea.encode(encode_M[i],p[i]);
    // }
    ZZX M[packednum];
    int CC[packednum*n];
    for(int i=0;i<lenth_C;i++)
        CC[i]=C[i];
    for(int i=lenth_C; i<n*packednum;i++)
        CC[i]=0;
    for(int i=0;i<packednum;i++){
        M[i].SetLength(n);
        SetCoeff(M[i],0,1);
        for(int z=0;z<n-1;z++){
		    SetCoeff(M[i],z+1,CC[i*(n-1)+z]);           
        }
    }  

    ZZX MM[packednum]; 
     
    //u,vand c, k=2;
	ZZX** u1= new ZZX*[packednum]();//input M has 13(outH*outW) cipher,and each cipher is vector and each element is a poly with n degree
    for (int i = 0; i < packednum; i++) {
        u1[i] = new ZZX[k]();
    }
    ZZX v1[packednum];   
    
    timeval t1,t2;
    gettimeofday(&t1,NULL);
     //Enc
     cout<<packednum<<endl;
    Enc enc;
    for(int i=0;i<packednum;i++){        
        MM[i]=U.reangecoeff(M[i],n);  
	    enc.encrypt(MM[i],u1[i],v1[i],n,q,k,q2d,usp[i].r1,A,usp[i].e11,usp[i].e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv);
    }
    gettimeofday(&t2,NULL);
    double deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    // cout << "FC层编码和加密时间为:" << deltaT1/1000<<"ms"<< endl;
 
    //store inputcipher file
    ofstream outFile1;
    outFile1.open(fc_inputcipherfile,ios::ate);
    for(int i=0;i<packednum;i++){
        for (int j=0;j<k;j++){
            outFile1<<u1[i][j]<<endl;
            }     
        outFile1<<v1[i]<<endl; 
	}  
    outFile1.close();         
    U.deleteZZXArray2(u1,packednum,k);
}

void Client_cipher::dense_input_cipherGen(int* C, int lenth_C,string fc_inputcipherfile, helib::EncryptedArray ea)
{  
    util U;
    //timeval p1, p2;
    // gettimeofday(&p1, NULL);   
    int packednum=ceil(1.0*lenth_C/(n-1));       
    // int slots=ea.size();
    // ZZX encode_M[packednum];//for f1:13 ploy and each poly's degree n is 256
    // std::vector <vector<long>>p(packednum,vector<long>(slots,0));

	// for(int i=0;i<packednum;i++){
    //     p[i][0]=1;
    //     for(int j=0;j<slots-1;j++){
    //         p[i][j+1]=C[i*(slots-1)+j];
    //     }       
    //     // ea.encode(encode_M[i],p[i]);
    // }
     ZZX M[packednum];
    int CC[packednum*n];
    for(int i=0;i<lenth_C;i++)
        CC[i]=C[i];
    for(int i=lenth_C; i<n*packednum;i++)
        CC[i]=0;
    for(int i=0;i<packednum;i++){
        M[i].SetLength(n);
        SetCoeff(M[i],0,1);
        for(int z=0;z<n-1;z++){
		    SetCoeff(M[i],z+1,CC[i*(n-1)+z]);           
        }
    }

   
    
   
    //加密M加密前没有编码
    //read t in pk.csv 

	ZZX t[k];
    ifstream inFilepk;
    inFilepk.open("./file/key/pk-t.csv",ios::in);  
    if(!inFilepk){
        cout<<"Cannot open pk-t!"<<endl;
        exit(0);
         }  
        for(int i=0;i<k;i++){
        inFilepk>>t[i];       
    }    
    inFilepk.close();	


	for(int i=0;i<k;i++){
        U.decompress(t[i],n,q,dt);
    }

	//read A in pk-A
    ZZX** A = new ZZX*[k]();
    for (int i = 0; i < k; i++) {
        A[i] = new ZZX[k]();
    }
   
	//ZZX A[k][k];
    ifstream inFile2;
    inFile2.open("./file/key/pk-A.csv",ios::in);
     if(!inFile2){
        cout<<"Cannot open pk-A!"<<endl;
        exit(0);
         } 
    for(int i=0;i<k;i++){
        for(int j=0;j<k;j++){
            inFile2>>A[i][j];     
            //cout<<A[i][j]<<" ";  
        }
	}    
    inFile2.close(); 
	
    ZZ_p::init(q);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,n,1);
    ZZ_pXModulus pm(ppm);
    ZZ_pX p_temp1,p_temp2,p_temp3,p_temp4;
    //Enc
    Enc enc;
    //q/2d
    ZZ two_d;
    RR two_d_rr;
    power2(two_d,dp);//two_d=2^dp
    two_d_rr=to_RR(two_d);
    RR q_rr;
    q_rr = to_RR(q);
    ZZ q2d;
    q2d=RoundToZZ(q_rr/two_d_rr);//q2d=round(q/2^dp)

     ZZX MM[packednum];
    Client_util cu;
    enc_param usp[packednum];
    for(int i=0;i<packednum;i++)
        usp[i]=U.randnum_gen();
	
     
    //u,vand c, k=2;
	ZZX** u1= new ZZX*[packednum]();//input M has 13(outH*outW) cipher,and each cipher is vector and each element is a poly with n degree
    for (int i = 0; i < packednum; i++) {
        u1[i] = new ZZX[k]();
    }
    ZZX v1[packednum];   
    
    timeval t1,t2;
    gettimeofday(&t1,NULL);

    for(int i=0;i<packednum;i++){        
        MM[i]=U.reangecoeff(M[i],n);  
	    enc.encrypt(MM[i],u1[i],v1[i],n,q,k,q2d,usp[i].r1,A,usp[i].e11,usp[i].e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv);
    }
    gettimeofday(&t2,NULL);
    double deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    // cout << "FC层编码和加密时间为:" << deltaT1/1000<<"ms"<< endl;
 
         //store inputcipher file
    ofstream outFile1;
    outFile1.open(fc_inputcipherfile,ios::ate);
    for(int i=0;i<packednum;i++){
        for (int j=0;j<k;j++){
            outFile1<<u1[i][j]<<endl;
            }     
        outFile1<<v1[i]<<endl; 
	}  
        outFile1.close();  
       
        U.deleteZZXArray2(A,k,k);
        U.deleteZZXArray2(u1,packednum,k);
}

/*
对同态计算的结果进行解密，解密结果是[kerN][outH][outW]的fature map
*/
int*** Client_cipher::DecCNN2D(string evalfile,string skfile,int kerN,int kerH,int kerW,int inH, int inW,int strid_H,int strid_W,string padding)
{   
    util U;  
    int outH=ceil(inH/strid_H);
    int outW=ceil(inW/strid_W);
    int packed_num1 = n/(kerH*kerW+1);
    int packed_col1 =ceil(float(outH*outW)/packed_num1); 
    int ***relu_inner=U.newArray3(kerN,outH,outW);
    if (padding=="SAME")
    {    
        int homcol=packed_col1;
        clock_t startTime11;
        clock_t endTime11;
        startTime11=clock();
        //read homomorphic eval result
        ifstream inFile1;
        inFile1.open(evalfile,ios::in);
        if(!inFile1){
            cout<<"Cannot open homomorphic eval result file!"<<endl;
            exit(0);
        }
        
        ZZX h1[kerN][homcol],h2[kerN][homcol],h3[kerN][homcol];    

        for(int i=0;i<kerN;i++){
        for (int j=0;j<homcol;j++){   
            inFile1>>h1[i][j];
            inFile1>>h2[i][j];
            inFile1>>h3[i][j];              
            }         
        }
        inFile1.close();
    
        // endTime11=clock();
        // cout<<"Homomorphic res read time is "<<(double)(endTime11-startTime11)/CLOCKS_PER_SEC<<"s"<<endl;

    
        //用户解密第一通道计算结果。
        clock_t startTime1;
        clock_t endTime1;

        startTime1=clock();
        Client_util cu;
        ZZX WX1[kerN][packed_col1];  //第一通道所有多项式相乘计算结果  
        ZZ* p1=new ZZ[packed_num1];//packed_num1=9
        ZZ dec_Mat1[kerN][packed_col1*packed_num1];
        
        //第一通道解包后的内积结果 
        for(int i=0;i<kerN;i++){   
            for (int j=0;j<packed_col1;j++){ 
                WX1[i][j]=cu.dec_onepackedPloy(skfile,h1[i][j],h2[i][j],h3[i][j]);
                p1=cu.GetpakedPCoeff_bias(WX1[i][j],kerH,kerW,packed_num1);
                for(int z=0;z<packed_num1;z++){
                    dec_Mat1[i][packed_num1*j+z]=p1[z];              
                }     
            }
        }
        delete []p1;
        // endTime1=clock();
        // cout<<"CONV1 dec time is "<<(double)(endTime1-startTime1)/CLOCKS_PER_SEC<<"s"<<endl;

        //取outH*outW列得到一通道的矩阵
        ZZ dec_finalMat[kerN][outH*outW];   
    
        for(int i=0;i<kerN;i++){
            for(int j=0;j<outH*outW;j++){
                dec_finalMat[i][j]=dec_Mat1[i][j];
            }
        }
    // innerfile.close();

        //解密的二维矩阵dec_finalMat：kerN×(outH*outW)转化为三维dec_3DfinalMat
        // kerN×outH×outW
        ZZ dec_3DfinalMat[kerN][outH][outW];
        for(int i=0;i<kerN;i++){
            int j=0;
            for(int u1=0;u1<outH;u1++){
            for(int u2=0;u2<outW;u2++){        
                    dec_3DfinalMat[i][u1][u2]=dec_finalMat[i][j];
                    j++;
                }     
            }
        }
    /*
        cout<<"// 测试二维变三维..."<<endl;
        cout<<"dec_finalMat[10][j]=";
        for(int j=0;j<outH*outW;j++){
                cout<<dec_finalMat[10][j]<<" ";
            }
        cout<<endl;   

        cout<<"dec_3DfinalMat[10][j]="<<endl;
        for(int u1=0;u1<outH;u1++){
            for(int u2=0;u2<outW;u2++){        
                    cout<<dec_3DfinalMat[10][u1][u2]<<" "; 
                }  
                cout<<endl;   
        }
        cout<<endl;*/

        //释放内存
        // U.deleteZZXArray2(v1v2,kerN,homcol);
        // U.deleteZZXArray2(v2u10,kerN,homcol);
        // U.deleteZZXArray2(v1u20,kerN,homcol);
        // U.deleteZZXArray2(v2u11,kerN,homcol);
        // U.deleteZZXArray2(v1u21,kerN,homcol);
        // U.deleteZZXArray2(u10u21,kerN,homcol);
        // U.deleteZZXArray2(u11u20,kerN,homcol);
        // U.deleteZZXArray2(u10u20,kerN,homcol);
        // U.deleteZZXArray2(u11u21,kerN,homcol);  
        
        timeval relu1,relu2;
        gettimeofday(&relu1,NULL);
        double double_inner[kerN][outH][outW];      
        
        for(int i=0;i<kerN;i++){
            for(int u1=0;u1<outH;u1++){
            for(int u2=0;u2<outW;u2++){ 
                double_inner[i][u1][u2]=conv<double>(dec_3DfinalMat[i][u1][u2]);
                relu_inner[i][u1][u2]=cu.relu(double_inner[i][u1][u2]*double(1.0/w_precision));
            }
            }
        }
        // gettimeofday(&relu2,NULL);
        // double deltaT = (relu2.tv_sec - relu1.tv_sec) * 1000000 + relu2.tv_usec - relu1.tv_usec;
        // cout<<"relu time is "<<deltaT/1000<<"ms"<<endl;    
    }
   return relu_inner;
}
 
    


int*** Client_cipher::DecCNN3DHaveAdd(string evalfile,string skfile,int kerN,int kerD,int kerH,int kerW,int inH, int inW,int strid_H,int strid_W,string padding)
{   
    util U;
    timeval t1, t2;
    int outH=ceil(inH/strid_H);
    int outW=ceil(inW/strid_W);
    int packed_num1 = n/(kerH*kerW+1);
    int packed_col1 =ceil(double(outH*outW)/packed_num1);
    int packed_num2 = n/(kerH*kerW);
    int packed_col2 =ceil(double(outH*outW)/packed_num2);
    int ***relu_inner=U.newArray3(kerN,outH,outW); 
    double double_inner[kerN][outH][outW];  
    if(padding=="SAME")
    {          
        ZZX h1[kerN][packed_col1],h2[kerN][packed_col1],h3[kerN][packed_col1];    
        ZZX CsumMat1[kerN][packed_col2];
        ZZX CsumMat2[kerN][packed_col2];
        ZZX CsumMat3[kerN][packed_col2];

        gettimeofday(&t1,NULL);
        ifstream inFile1;
        inFile1.open(evalfile,ios::in);
        if(!inFile1){
            cout<<"Cannot open homomorphic eval result file!---"<<k<<endl;
            exit(0);
        }   
        for(int i=0; i<kerN; i++){
            for(int j=0; j<packed_col1; j++){
                inFile1>>h1[i][j];
                inFile1>>h2[i][j];
                inFile1>>h3[i][j];
            }
        }

        for(int i=0; i<kerN; i++){
            for(int j=0; j<packed_col2; j++){
                inFile1>>CsumMat1[i][j];
                inFile1>>CsumMat2[i][j];
                inFile1>>CsumMat3[i][j];
            }
        }             
    
        gettimeofday(&t2, NULL);
        // double deltaT2 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
        // cout<<"Homomorphic res read time is "<<deltaT2/1000.0<<"ms"<<endl;      
    
        //用户解密第一通道计算结果。
        gettimeofday(&t1, NULL);
        Client_util cu;
        ZZX WX1[kerN][packed_col1];  //第一通道所有多项式相乘计算结果  
        ZZ* p1=new ZZ[packed_num1];//packed_num1=9
        ZZ dec_Mat1[kerN][packed_col1*packed_num1];
        
        //第一通道解包后的内积结果 
        for(int i=0;i<kerN;i++){   
            for (int j=0;j<packed_col1;j++){ 
                WX1[i][j]=cu.dec_onepackedPloy(skfile,h1[i][j],h2[i][j],h3[i][j]);
                p1=cu.GetpakedPCoeff_bias(WX1[i][j],kerH,kerW,packed_num1);
                for(int z=0;z<packed_num1;z++){
                    dec_Mat1[i][packed_num1*j+z]=p1[z];              
                }     
            }
        }
        delete []p1;      
        
        //用户解密其余kerD-1个通道的计算结果。
        ZZX WX2[kerN][packed_col2];
        ZZ*p=new ZZ[packed_num2];
        int col2=packed_col2*packed_num2;     
        ZZ **dec_Mat2=new ZZ*[kerN];
        for(int i=0;i<kerN;i++){
            dec_Mat2[i]=new ZZ[col2];
            for(int j=0;j<col2;j++)
            dec_Mat2[i][j]=ZZ(0);
        }   
       
        for(int i=0;i<kerN;i++){  
            for (int j=0;j<packed_col2;j++){           
                WX2[i][j]=cu.dec_onepackedPloy(skfile,CsumMat1[i][j],CsumMat2[i][j],CsumMat3[i][j]);
                p=cu.GetpakedPCoeff_Nobias(WX2[i][j],kerH,kerW,packed_num2);
                for(int z=0;z<packed_num2;z++){
                    dec_Mat2[i][packed_num2*j+z]=p[z];                      
                }                 
            }
        }
        gettimeofday(&t2, NULL);
        // double deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
        
        // cout<<"CONV dec time is "<<(double)deltaT/1000.0<<"ms"<<endl;
        delete[] p; 
        
        int cols_Num = outH*outW; 
        ZZ   dec_feartureMat[kerN][cols_Num];  
        //2个矩阵取outH*outW列，对应位置相加得到矩阵dec_sumMat:kerN×(outH*outW)

        for(int i=0; i<kerN; i++){
            for(int j=0; j<cols_Num; j++){
                dec_feartureMat[i][j]=dec_Mat1[i][j]+dec_Mat2[i][j];
            }
        }
        // ofstream innerfile;
        // innerfile.open(decplainfile,ios::ate);
        // for(int i=0;i<kerN;i++){
        //     for (int j=0;j<cols_Num;j++){
        //         innerfile<<dec_sumMat[i][j]<<" ";
        //     } 
        // }
        // innerfile.close();
    
        //释放内存
        for (int i = 0; i < kerN; i++) {
            delete[] dec_Mat2[i];
        }
        delete[]  dec_Mat2;          
   
        timeval relu1,relu2;
        gettimeofday(&relu1,NULL);  
        for(int i=0;i<kerN;i++){
            int j=0;
            for(int u1=0;u1<outH;u1++){
                for(int u2=0;u2<outW;u2++){ 
                    double_inner[i][u1][u2]=conv<double>(dec_feartureMat[i][j]);
                    relu_inner[i][u1][u2]=cu.relu(double_inner[i][u1][u2]*double(1.0/w_precision));
                    j=j+1;
                }
            }
        }
        // gettimeofday(&relu2,NULL);
        // double deltaT1 = (relu2.tv_sec - relu1.tv_sec) * 1000000 + relu2.tv_usec - relu1.tv_usec;
        // cout<<"relu time is "<<deltaT1/1000<<"ms"<<endl;         
      
    }  
    return relu_inner;       
}

/*
对同态计算的结果进行解密，解密结果是[kerN][outH][outW]的fature map
*/
int* Client_cipher::DecDenseHaveAdd_relu(int den_num1,int fc1_w,string resfile)
{   util U;
    int a=ceil(double(fc1_w+1)/n);
    // cout<<"In dec packednum="<<a<<endl;
    //int a=13;
    // int den_num1=1024;
    ZZX* sum1=new ZZX[den_num1];  
    ZZX* sum2=new ZZX[den_num1];  
    ZZX* sum3=new ZZX[den_num1];
    // ZZX* sum4=new ZZX[den_num1];
    // ZZX* sum5=new ZZX[den_num1];
    // ZZX* sum6=new ZZX[den_num1];
    // ZZX* sum7=new ZZX[den_num1];
    // ZZX* sum8=new ZZX[den_num1];
    // ZZX* sum9=new ZZX[den_num1];

     //read fc2_res.csv
    ifstream inFile1;
	inFile1.open(resfile,ios::in); 
    if(!inFile1){
        cout<<"Cannot open resfile in Dec!"<<endl;
        exit(0);
    }  
    for(int i=0;i<den_num1;i++){
        inFile1>>sum1[i];
        inFile1>>sum2[i];
        inFile1>>sum3[i];
        // inFile1>>sum4[i];
        // inFile1>>sum5[i];
        // inFile1>>sum6[i];
        // inFile1>>sum7[i];
        // inFile1>>sum8[i];    
        // inFile1>>sum9[i];           
    }
    inFile1.close();
 
    Client_util cu;
    ZZ inner[den_num1];    
    //用户解密全连接1计算结果。
    timeval t1, t2;
   
    gettimeofday(&t1,NULL);
    for(int i=0;i<den_num1;i++){
        inner[i]=cu.dec_one_poly("./file/key/sk.csv",sum1[i],sum2[i],sum3[i]);
   // cout<<inner[i]<<" ";
    }  
   // cout<<endl; 
    gettimeofday(&t2,NULL);
    double deltaT4 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    cout<<"FC dec time is:"<<deltaT4/1000<<"ms"<<endl<<endl;

    timeval relu1,relu2;
    gettimeofday(&relu1,NULL);
    double double_inner[den_num1];
    int *relu_inner=new int[den_num1];

    for (int j=0;j<den_num1;j++){
        double_inner[j]=conv<double>(inner[j]);    
        relu_inner[j]=cu.relu(double_inner[j]*double(1.0/w_precision));
    } 
    gettimeofday(&relu2,NULL);
    double deltaT3 = (relu2.tv_sec - relu1.tv_sec) * 1000000 + relu2.tv_usec - relu1.tv_usec;
    cout<<"relu time is "<<deltaT3/1000<<"ms"<<endl;   
   
   delete[] sum1,sum2,sum3;
   return relu_inner;
   
}

double Client_cipher::DecDenseHaveAdd_softmax(int den_num1,int fc1_w,string resfile)
{   
    util U;
    int a=ceil(double(fc1_w+1)/n);
    // cout<<"In dec packednum="<<a<<endl;
    //int a=13;
  
    ZZX* sum1=new ZZX[den_num1];  
    ZZX* sum2=new ZZX[den_num1];  
    ZZX* sum3=new ZZX[den_num1];

     //read fc2_res.csv
    ifstream inFile1;
	inFile1.open(resfile,ios::in); 
    if(!inFile1){
        cout<<"Cannot open resfile in Dec!"<<endl;
        exit(0);
    }  
    for(int i=0;i<den_num1;i++){
        inFile1>>sum1[i];
        inFile1>>sum2[i];
        inFile1>>sum3[i];
             
     }
    inFile1.close();
 
    Client_util cu;
    ZZ inner[den_num1];    
    //用户解密全连接1计算结果。
  
   timeval t1,t2;
    gettimeofday(&t1,NULL);

    for(int i=0;i<den_num1;i++){
        inner[i]=cu.dec_one_poly("./file/key/sk.csv",sum1[i],sum2[i],sum3[i]);
   // cout<<inner[i]<<" ";
    }  
    //cout<<endl;

    double double_inner[den_num1];
    
    for (int j=0;j<den_num1;j++){
        double_inner[j]=(conv<double>(inner[j]));      
    } 
   
            gettimeofday(&t2,NULL);
    double deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    cout << "FC解密时间为:" << deltaT1/1000<<"ms"<< endl;
   
    // cout<<"===================================="<<endl;
    double *softmax_inner=new double[den_num1];
    for(int i=0;i<den_num1;i++){
        softmax_inner[i] = 0;
        // cout<<double_inner[i]<<" ";
    }
    cout<<endl;

    timeval relu1,relu2;
    gettimeofday(&relu1,NULL);
    cu.softmax(double_inner, den_num1, softmax_inner);
    gettimeofday(&relu2, NULL);
    double deltaT = (relu2.tv_sec - relu1.tv_sec) * 1000000 + relu2.tv_usec - relu1.tv_usec;
    cout<<"softmax time is "<<deltaT/1000<<"ms"<<endl;
   
    double maxId = 0.0;
    double maxValue = 0.0;
    for(int i=0; i<den_num1; i++){
        if(maxValue < softmax_inner[i]){
            maxValue = softmax_inner[i];
            maxId = i;
        }
    }    
    // cout<<"prediction result:"<<maxId<<endl;
    // cout<<" probability="<<maxValue<<endl;     
    // cout<<"===================================="<<endl;       
    delete[] sum1,sum2,sum3;
    delete[] softmax_inner;
   return maxId;
}



