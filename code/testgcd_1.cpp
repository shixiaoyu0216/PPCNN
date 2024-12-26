#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/RR.h>
#include <helib/helib.h>
#include <param.h>
#include <ctime>
#include <util.h>
#include <vector>
#include <enc.h>
#include <Binomialdistribution.h>
#include <cloud_util.h>
using namespace std;
using namespace NTL;

int gcdExtended(int a, int b,int *x, int *y){
    if(a==0){
        *x=0;
        *y=1;
        return b;
    }
    int x1,y1;
    int gcd=gcdExtended(b%a, a, &x1, &y1);
    *x=y1-(b/a)*x1;
    *y=x1;
    return gcd;
}

int main()
{ 
    long randnum=2;     
    long m=512, p=7681,r=1;    
    ZZX G;
    G.SetLength(2);
    SetCoeff(G,1,1);
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                            .m(m)
                            .p(p)
                            .r(r)
                            .build();
    // Print the context
    helib::EncryptedArray ea(context,G);
    int slots=ea.size();
    std::vector<long > nr(slots,0);
    std::vector<long > nr_inv(slots,0);
    std::vector<long > m1(slots,0);
  
    ZZX nr_poly;
    ZZX nr_invpoly;
    ZZX m1_poly;
    for(int i=0; i<slots; i++)
    {
        nr[i]=RandomBnd(randnum)%1000+1;
        // cout<<nr[i]<<" ";
        // nr[i]=2; //随机数为1
        int x=0,y=0;
        gcdExtended(nr[i],p,&x,&y);
        //nr*x=1 mod p,  p2和p互素，即p2取小于p的大于0的整数，p2的乘法模逆唯一。
        nr_inv[i]=x;          
        // cout<<endl<<nr_inv[i]<<" ";      
    }
    ea.encode(nr_poly, nr);  
    ea.encode(nr_invpoly, nr_inv); 
    cout<<"nr="<<endl; 
    for(int i=0;i<slots;i++)
        cout<<nr[i]<<" ";    
    cout<<endl;
    // cout<<endl<<nr_invpoly<<endl;
    
   
    // for (int i = 0; i < n; i++){
    //     m1[i]=-RandomBnd(800);
    //     cout<<m1[i]<<" ";
    //     // SetCoeff(m1,i,RandomBnd(1024));
    // }
    // cout<<endl;

    ifstream file("./file/maxpool.csv", ios::in);
    if(!file){
        cout<<"Cannot open weights!\n";
        exit(0);
    }
     cout<<"m1="<<endl;
    if (file.is_open()){    
        for (int j=0; j<16*16; j++){
                file>>m1[j];
            cout<<m1[j]<<" ";
        }
    } 
    file.close(); 
    ea.encode(m1_poly, m1);  
    cout<<"end m1"<<endl;

    Binomialdistribution Bin_dis;    
    ZZX s[k];
    for(int i=0;i<k;i++){
        s[i].SetLength(n);
        Bin_dis.Binomial_distribution_gen(s[i],n,ita);
        //        cout<<s[i]<<endl;
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
            //           cout<<A[i][j]<<endl;
        }
    }


    //Multinit
    ZZ_p::init(q);
    ZZ_pX pm;
    SetCoeff(pm,0,1);
    SetCoeff(pm,n,1);
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
    ZZX u2[k];
    ZZX v2;
  
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
    // cout<<"m1="<<m1<<endl;
    enc.encrypt(m1_poly,u2,v2,n,q,k,q2d,r1,A,e11,e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv);
    // ZZX m_dec1;
    // enc.decrypt(m_dec1,u2,v2,n,q,k,s,p_temp1,p_temp2,p_temp3,p_temp4,pm,du,dv,dp);
    // cout<<"m_dec1="<<m_dec1<<endl;
    // cout<<m_dec1-m1_poly<<endl;


  
    //u,vand c, k=2;
    ZZX u3[k];
    ZZX v3;
   
    // ZZX R_poly;
    // std::vector<long> p1(slots,0);    
    // ea.encode(R_poly, p1);   
    ZZX R;  
    R.SetLength(n); 
    enc_param usp;
    Cloud_util cdu;
    
    usp = U.randnum_gen(); 
    enc.encrypt(R,u3,v3,n,q,k,q2d,usp.r1,A,usp.e11,usp.e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv); 

    //减去maxpool的噪声

    u2[0]=u2[0]+u3[0];
    u2[1]=u2[1]+u3[1];
    v2=v2+v3;
    U.Module(u2[0],slots,q);
    U.Module(u2[1],slots,q);
    U.Module(v2,slots,q);
      
    for(int j=0;j<k;j++){
        U.decompress(u2[j],n,q,du);
    }
    U.decompress(v2,n,q,dv);
 
    
    //乘上用于求relu的噪声
    ZZ_p::init(q2);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,slots,1);
    ZZ_pX pp_temp1,pp_temp2,pp_temp3,pp_temp4; 
    
    ZZX mu0, mu1, mv;    
    U.Mult(mu0,ppm,nr_poly,u2[0],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    U.Mult(mu1,ppm,nr_poly,u2[1],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    U.Mult(mv,ppm,nr_poly,v2,slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);

    ZZX s0u10w1,s1u11w1,h6; 
    RR two_dp_rr;
    power2(two_dp_rr,dp);
    RR q2_rr;
    ZZ two_d2;
    ZZ two_d3;
    two_d2=two_d/2;
    two_d3=two_d2*(-1);
    q2_rr = q_rr;


    U.Mult(s0u10w1,ppm,s[0],mu0,slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    U.Mult(s1u11w1,ppm,s[1],mu1,slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    ZZX h5=mv-s0u10w1-s1u11w1; 
    for(int j=0;j<n;j++){
        RR coeff_rr;
        ZZ coeff_zz;
        GetCoeff(coeff_zz,h5,j);
        coeff_rr=to_RR(coeff_zz);
        coeff_zz=RoundToZZ(coeff_rr*two_dp_rr/q2_rr)%two_d;
        if(coeff_zz>=two_d2){
            coeff_zz= coeff_zz-two_d;
        }
        else if(coeff_zz<=two_d3){
            coeff_zz=coeff_zz+two_d;
        }
        SetCoeff(h5,j,coeff_zz);           
    }  
    U.Module(h5,n,conv<ZZ>(p)); 
   
   //解码再relu
    std::vector<long> pp(slots,0);   
    ea.decode(pp,h5);
    
    // cout<<"h5="<<h5<<endl;
    for(int i=0;i<slots;i++){
       pp[i]=pp[i]%p;
        if(pp[i]>=p/2){
            pp[i]= pp[i]-p;
        }
        else if(pp[i]<=-p/2){
            pp[i]=pp[i]+p;
        }
          cout<<pp[i]<<" ";
    }
    cout<<"end pp"<<endl;

    cout<<"res="<<endl;
    cout<<"p="<<p<<endl;
    std::vector<long> res(slots,0); 
    for(int i=0;i<slots;i++){       
        res[i]=nr[i]*m1[i]-pp[i];
        res[i]=res[i]%p;
        if(res[i]>=p/2){
            res[i]= res[i]-p;
        }
        else if(res[i]<=-p/2){
            res[i]=res[i]+p;
        }
        cout<<res[i]<<" ";
    }
    cout<<endl;
    cout<<q2<<endl;
    cout<<q<<endl;
    cout<<du<<endl;
    cout<<dv<<endl;
    cout<<p<<endl;    
   return 0;
    }