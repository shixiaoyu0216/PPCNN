#include "enc.h"

Enc::Enc()
{

}

void Enc::encrypt(ZZX m, ZZX u[], ZZX &v, int n, ZZ q, int k, ZZ q2d, ZZX r[],ZZX **A,ZZX e1[],ZZX e2,ZZ_pX p_temp1, ZZ_pX p_temp2, ZZ_pX p_temp3, ZZ_pX p_temp4,ZZ_pX pm,ZZX t[],int dt,int du,int dv){
    util U;
    for (int i = 0; i < k; i++) {
        ZZX polytemp1;      
        for (int j = 0; j < k; j++) {
            U.Mult(polytemp1,pm,A[j][i],r[j],n,q,p_temp1,p_temp2,p_temp3,p_temp4);
            u[i] = u[i]+polytemp1;           
        }
        u[i] = u[i]+e1[i];        
        U.Module(u[i],n,q);      
        U.compress(u[i],n,q,du);       
    }
    ZZX polytemp2;
    for (int j = 0; j < k; j++) {
        U.Mult(polytemp2,pm,t[j],r[j],n,q,p_temp1,p_temp2,p_temp3,p_temp4);
        v = v+polytemp2;
    }
    v = v+e2;
    v= v + m * q2d;
    U.Module(v,n,q);
    U.compress(v,n,q,dv);
}

void Enc::decrypt(ZZX& m, ZZX u[], ZZX &v, int n, ZZ q, int k, ZZX s[], ZZ_pX p_temp1, ZZ_pX p_temp2, ZZ_pX p_temp3, ZZ_pX p_temp4, ZZ_pX pm, int du, int dv, int dp){
    util U;
    ZZX stu;
    ZZX polytemp3;
 
    for (int i = 0; i < k; i++) {
        U.decompress(u[i],n,q,du);
        U.Mult(polytemp3,pm,s[i],u[i],n,q,p_temp1,p_temp2,p_temp3,p_temp4);
        stu = stu+polytemp3;
        U.Module(stu,n,q);
    }
    stu=-1*stu;
    U.decompress(v,n,q,dv);
    m = v+stu;
    U.Module(m,n,q);
    U.compress(m,n,q,dp);
    //make plain can be negative
 
}
