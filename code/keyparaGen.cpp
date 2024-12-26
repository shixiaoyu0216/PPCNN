#include <Binomialdistribution.h>
#include <param.h>
#include <keyparaGen.h>
#include <iostream>
#include<fstream>

KeyparaGen::KeyparaGen()
{

}
void KeyparaGen::keygen(){   
   //skgen   
    Binomialdistribution Bin_dis;
    srand((unsigned) time(NULL));
    ZZX s[k];
    for(int i=0;i<k;i++){
        s[i].SetLength(n);
        Bin_dis.Binomial_distribution_gen(s[i],n,ita);
        //        cout<<s[i]<<endl;
    }
    ofstream outFile1;
    outFile1.open("./file/key/sk.csv",ios::ate);    
         for(int i=0;i<k;i++){
         outFile1<<s[i]<<endl;    
   }    
    outFile1.close();
    
    //pkgen
    util U;
    ZZX e[k];
    for(int i=0;i<k;i++){
        e[i].SetLength(n);
        Bin_dis.Binomial_distribution_gen(e[i],n,ita);
    }

    ZZX** A = new ZZX*[k]();
    for (int i = 0; i < k; i++) {
        A[i] = new ZZX[k]();
    }
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            A[i][j].SetLength(n);
            U.rand(A[i][j],n,q);
            //           cout<<A[i][j]<<endl;
        }
    }
  //Multinit
    ZZ_p::init(q);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,n,1);
    ZZ_pXModulus pm(ppm);
    ZZ_pX p_temp1,p_temp2,p_temp3,p_temp4;

    //gen t
    ZZX *t = new ZZX[k]();
    for(int i = 0; i < k; i++){
        ZZX polytemp1;
        for (int j = 0; j < k; j++) {
            U.Mult(polytemp1,pm,A[i][j],s[j],n,q,p_temp1,p_temp2,p_temp3,p_temp4);
            t[i]=t[i]+polytemp1;
        }
        t[i]=t[i]+e[i];
        U.Module(t[i],n,q);
        U.compress(t[i],n,q,dt);
    }
    ofstream outFile2;
    outFile2.open("./file/key/pk-A.csv",ios::ate);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
    outFile2<<A[i][j]<<" "<<endl;}
    }
     outFile2.close();

     ofstream outFile3;
    outFile3.open("./file/key/pk-t.csv",ios::ate);
    for (int z=0;z<k;z++){
        outFile3<<t[z]<<" "<<endl;
    }
    outFile3.close();
    U.deleteZZXArray2(A,k,k);
    delete []t;
}

    


   