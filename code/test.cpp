#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/RR.h>
#include <Binomialdistribution.h>
#include <util.h>
#include <enc.h>
#include <ctime>
using namespace std;
using namespace NTL;

int main()
{
    //1.params
    int n = 256;
    int k = 2;
    int dt = 79;
    int du = 79;
    int dv = 79;
    int dp = 29;
    int ita = 5;
    ZZ q = conv<ZZ>("4835703278458516698824713");
    ZZ q2 = conv<ZZ>("11972621413014756705924586149611790497021399392059469");
    /*int dt = 65;
    int du = 65;
    int dv = 65;
    int dp = 23;
    int ita = 3;
    ZZ q = conv<ZZ>("73786976294838206473");
    ZZ q2 = conv<ZZ>("4835703278458516698824713");*/
    ZZX m1;
    m1.SetLength(n);
    for (int i = 0; i < n; i++){
        SetCoeff(m1,i,RandomBnd(1000));
    }

    cout<<"随机生成的第一个整数向量为:"<<endl;
    cout<<m1<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;

    

    ZZX m2;
    m2.SetLength(n); 
    for (int i = 0; i < n; i++){
        SetCoeff(m2,i,RandomBnd(1000));
    }
    cout<<"随机生成的第二个整数向量为:"<<endl;
    cout<<m2<<endl;
    // ZZX m2;
    // m2.SetLength(n);
    // SetCoeff(m2,0,1);
    // for (int i = 1; i < n; i++){
    //     SetCoeff(m2,i,0);
    // }
    // cout<<"随机生成的第二个整数向量为:"<<endl;
    // cout<<m2<<endl;

    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;

    ZZ temp,temp2;
    ZZ Sum(0);
    for(int i=0;i<n;i++){
        GetCoeff(temp,m1,i);
        GetCoeff(temp2,m2,i);
        Sum+=temp*temp2;
    }
    cout<<"两个向量的明文内积为:"<<Sum<<endl;
    cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;

    ZZ tempm2;

    ZZX m3;
    m3.SetLength(n);
    int j=1;
    GetCoeff(tempm2,m2,0);
    SetCoeff(m3,0,tempm2);
    for(int i=n-1;i>0;i--){
        GetCoeff(tempm2,m2,i);
        SetCoeff(m3,j,-1*tempm2);
        j++;
    }
    m2=m3;
    cout<<"将第二个向量以特殊的方式填充:"<<endl;
    cout<<m2<<endl;
cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;

    //new
   //2.skgen
    Binomialdistribution Bin_dis;
    srand((unsigned) time(NULL));
    ZZX s[k];
    for(int i=0;i<k;i++){
        s[i].SetLength(n);
        Bin_dis.Binomial_distribution_gen(s[i],n,ita);
      
    }
    //3.pkgen
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
         
        }
    }


    //Multinit
    ZZ_p::init(q);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,n,1);
    ZZ_pXModulus pm(ppm);
    ZZ_pX p_temp1,p_temp2,p_temp3,p_temp4;

    ZZX mm;
    U.Mult(mm,pm,m1,m2,n,q,p_temp1,p_temp2,p_temp3,p_temp4);


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

    //r1,e11,e12
    ZZX r1[k];
    for(int i=0;i<k;i++){
        r1[i].SetLength(n);
        Bin_dis.Binomial_distribution_gen(r1[i],n,ita);
    }
    ZZX e11[k];
    for(int i=0;i<k;i++){
        e11[i].SetLength(n);
        Bin_dis.Binomial_distribution_gen(e11[i],n,ita);
    }
    ZZX e12;
    Bin_dis.Binomial_distribution_gen(e12,n,ita);

    //m and c
    ZZX u1[k];
    ZZX v1;
    clock_t startTime2;
    clock_t endTime2;
    startTime2 = clock();
    //Enc
    Enc enc;
    //q/2d
    ZZ two_d;
    RR two_d_rr;
    power2(two_d,dp);
    two_d_rr=to_RR(two_d);
    RR q_rr;
    q_rr = to_RR(q);
    ZZ q2d;
    q2d=RoundToZZ(q_rr/two_d_rr);
    for(int i=0;i<k;i++){
        U.decompress(t[i],n,q,dt);
    }

    enc.encrypt(m1,u1,v1,n,q,k,q2d,r1,A,e11,e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv);
    endTime2 = clock();
    //cout << "time:" << (double) (endTime2 - startTime2) / CLOCKS_PER_SEC << endl;
    /* //Dec
    ZZX m_dec;
    enc.decrypt(m_dec,u1,v1,n,q,k,s,p_temp1,p_temp2,p_temp3,p_temp4,pm,du,dv,dp);
    cout<<m_dec-m1;*/   

    for(int i = 0;i<k;i++){
        U.decompress(u1[i],n,q,du);
       
    }
    U.decompress(v1,n,q,dv);
   
    //Multinit2
    ZZ_p::init(q2);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,n,1);
    ZZ_pXModulus pm(ppm);
    ZZ_pX pp_temp1,pp_temp2,pp_temp3,pp_temp4;

 
    clock_t startTime;
    clock_t endTime;
    startTime = clock();
   
    // ZZX p1=v1-s0u10-s1u11;
    // ZZX p2;
    // U.Mult(p2,pm,p1,m2,n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    ZZX p1,p2,p3;    
    U.Mult(p1,pm,v1,m2,n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    U.Mult(p2,pm,u1[0],m2,n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    U.Mult(p3,pm,u1[1],m2,n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);


    ZZX s0u10,s1u11;
     U.Mult(s0u10,pm,s[0],p2,n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    U.Mult(s1u11,pm,s[1],p3,n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    ZZX p4=p1-s0u10-s1u11;

     

   
    endTime = clock();
    cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;
    cout << "计算同态内积时间为:" << (double) (endTime - startTime) / CLOCKS_PER_SEC <<"s"<< endl;
    cout<<endl;cout<<endl;
    RR two_dp_rr;
    power2(two_dp_rr,dp);
    RR q2_rr;
    ZZ two_d2;
    ZZ two_d3;
    two_d2=two_d/2;
    two_d3=two_d2*(-1);
    q2_rr = q_rr;
    for(int i=0;i<n;i++){
        RR coeff_rr;
        ZZ coeff_zz;
        GetCoeff(coeff_zz,p4,i);
        coeff_rr=to_RR(coeff_zz);
        coeff_zz=RoundToZZ(coeff_rr*two_dp_rr/q2_rr)%two_d;
        if(coeff_zz>=two_d2){
             coeff_zz= coeff_zz-two_d;
        }else if(coeff_zz<=two_d3){
            coeff_zz=coeff_zz+two_d;
        }
        SetCoeff(p4,i,coeff_zz);
    }
    ZZ inner;
    GetCoeff(inner,p4,0);
    cout<<"解密得到的结果为:"<<inner;
    cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;cout<<endl;

    return 0;
}

