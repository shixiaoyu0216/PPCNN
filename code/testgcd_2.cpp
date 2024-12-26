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
// #include <Binomialdistribution.h>
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
    int kerN=6;  
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
    std::vector<vector<long >> nr(kerN,vector<long>(slots,0));
    std::vector<vector<long >> nr_inv(kerN,vector<long>(slots,0));
    std::vector<vector<long >> m1(kerN,vector<long>(slots,0));
   
  
    ZZX nr_poly[kerN];
    ZZX nr_invpoly[kerN];
    ZZX m1_poly[kerN];
    for(int z=0; z<kerN; z++){
        for(int i=0; i<slots; i++)
        {
            nr[z][i]=RandomBnd(randnum)%1000+1;           
            // nr[z][i]=2; //随机数为1
            int x = 0, y = 0;
            gcdExtended(nr[z][i], p, &x, &y);
            //nr*x=1 mod p,  p2和p互素，即p2取小于p的大于0的整数，p2的乘法模逆唯一。
            nr_inv[z][i] = x;          
        }
        ea.encode(nr_poly[z], nr[z]);  
        ea.encode(nr_invpoly[z], nr_inv[z]); 
    }
    cout<<"nr="<<endl; 
    for(int i=0;i<slots;i++)
        cout<<nr[0][i]<<" ";    
    cout<<endl;
    // cout<<endl<<nr_invpoly<<endl;   
  

    ifstream file("./file/maxpool_2.csv", ios::in);
    if(!file){
        cout<<"Cannot open weights!\n";
        exit(0);
    }
     cout<<"m1="<<endl;
    if (file.is_open()){  
        for(int z=0; z<kerN; z++)  
        for (int j=0; j<16*16; j++){
                file>>m1[z][j];
            cout<<m1[z][j]<<" ";
        }
    } 
    file.close(); 
    for(int i=0; i<kerN; i++)
        ea.encode(m1_poly[i], m1[i]);  
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

    // //r1,e11,e12
    // ZZX r1[k];
    // for(int i=0;i<k;i++){
    //     r1[i].SetLength(n);
    //     Bin_dis.Binomial_distribution_gen(r1[i],n,ita);
    // }
    // ZZX e11[k];
    // for(int i=0;i<k;i++){
    //     e11[i].SetLength(n);
    //     Bin_dis.Binomial_distribution_gen(e11[i],n,ita);
    // }
    // ZZX e12;
    // Bin_dis.Binomial_distribution_gen(e12,n,ita);

     //m and c
    ZZX u2[kerN][k];
    ZZX v2[kerN];
  
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

    enc_param usp1;       
    usp1 = U.randnum_gen(); 
    for(int i=0; i<kerN; i++)
        enc.encrypt(m1_poly[i], u2[i], v2[i], n,q,k,q2d, usp1.r1, A, usp1.e11, usp1.e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv);
    // ZZX m_dec1;
    // enc.decrypt(m_dec1,u2,v2,n,q,k,s,p_temp1,p_temp2,p_temp3,p_temp4,pm,du,dv,dp);
    // cout<<"m_dec1="<<m_dec1<<endl;
    // cout<<m_dec1-m1_poly<<endl;


  
    //u,vand c, k=2;
    ZZX u3[kerN][k];
    ZZX v3[kerN];
   
    // ZZX R_poly;
    // std::vector<long> p1(slots,0);    
    // ea.encode(R_poly, p1);   
    ZZX R[kerN]; 
    cout<<"R[6]="<<endl;
    for(int i=0; i<kerN; i++)
    {
        R[i].SetLength(n);
        SetCoeff(R[i], i, RandomBnd(100));
        cout<<R[i]<<endl;
    }
    enc_param usp2;       
    usp2 = U.randnum_gen(); 
    for(int z=0; z<kerN; z++){
        // R[z].SetLength(n);    
        enc.encrypt(R[z], u3[z], v3[z],n,q,k,q2d,usp2.r1,A,usp2.e11,usp2.e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv); 
    }
    //减去maxpool的噪声
    
    for(int i=0; i<kerN; i++){
        u2[i][0] = u2[i][0] - u3[i][0];
        u2[i][1]= u2[i][1] - u3[i][1];
        v2[i] = v2[i] - v3[i];
        U.Module(u2[i][0], slots, q);
        U.Module(u2[i][1], slots, q);
        U.Module(v2[i], slots, q);
    }
    /*
    //解密解码
    ZZX h66[kerN];
    std::vector <vector<long > >pp2(kerN,vector<long>(slots,0));    
    cout<<"pp2="<<endl;  
    for(int i=0; i<kerN; i++){
        enc.decrypt(h66[i], u2[i], v2[i], slots,q,k,s,p_temp1,p_temp2,p_temp3,p_temp4,pm,du,dv,dp);
        U.Module(h66[i], slots, conv<ZZ>(p));      
        ea.decode(pp2[i], h66[i]);          
    }
    for(int z=0; z<kerN; z++){
        for(int i=0;i<slots;i++){
            pp2[z][i]=pp2[z][i]%p;
                if(pp2[z][i]>=p/2){
                    pp2[z][i]= pp2[z][i]-p;
                }
                else if(pp2[z][i]<=-p/2){
                    pp2[z][i]=pp2[z][i]+p;
                }
                cout<<pp2[z][i]<<" ";
        }
    }
    cout<<endl<<"end pp2"<<endl;
      */  
    for(int i=0; i<kerN; i++){
        for(int j=0;j<k;j++){
            U.decompress(u2[i][j], n, q, du);
        }
        U.decompress(v2[i], n, q, dv);
 
    }
    //乘上用于求relu的噪声
    ZZ_p::init(q2);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,slots,1);
    ZZ_pX pp_temp1,pp_temp2,pp_temp3,pp_temp4; 
    
    ZZX mu0[kerN], mu1[kerN], mv[kerN];   
    for(int i=0; i<kerN; i++){    
        U.Mult(mu0[i],ppm, nr_poly[i], u2[i][0], slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        U.Mult(mu1[i],ppm, nr_poly[i], u2[i][1], slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        U.Mult(mv[i],ppm,nr_poly[i], v2[i], slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    }
    ZZX s0u10w1[kerN], s1u11w1[kerN], h5[kerN]; 
    RR two_dp_rr;
    power2(two_dp_rr,dp);
    RR q2_rr;
    ZZ two_d2;
    ZZ two_d3;
    two_d2=two_d/2;
    two_d3=two_d2*(-1);
    q2_rr = q_rr;
    std::vector<vector<long>> pp(kerN, vector<long>(slots,0));  
    for(int i=0; i<kerN; i++){
        U.Mult(s0u10w1[i], ppm, s[0], mu0[i], slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        U.Mult(s1u11w1[i], ppm,s[1], mu1[i], slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        h5[i] = mv[i] - s0u10w1[i] - s1u11w1[i]; 
        for(int j=0; j<n; j++){
            RR coeff_rr;
            ZZ coeff_zz;
            GetCoeff(coeff_zz, h5[i], j);
            coeff_rr=to_RR(coeff_zz);
            coeff_zz=RoundToZZ(coeff_rr*two_dp_rr/q2_rr)%two_d;
            if(coeff_zz>=two_d2){
                coeff_zz= coeff_zz-two_d;
            }
            else if(coeff_zz<=two_d3){
                coeff_zz=coeff_zz+two_d;
            }
            SetCoeff(h5[i], j, coeff_zz);           
        }  
        U.Module(h5[i], n, conv<ZZ>(p)); 
    
        //解码再relu
         
        ea.decode(pp[i],h5[i]);       
        for(int z=0; z<slots; z++){
        pp[i][z] = pp[i][z] % p;
            if(pp[i][z] >= p/2){
                pp[i][z] = pp[i][z] - p;
            }
            else if(pp[i][z] <= -p/2){
                pp[i][z] = pp[i][z] + p;
            }
            cout<<pp[i][z]<<" ";
        }
    }
    cout<<"end pp"<<endl;

    cout<<"res="<<endl;
    cout<<"p="<<p<<endl;
    std::vector<vector<long>> res(kerN, vector<long>(slots,0)); 
    for(int z=0; z<kerN; z++){
        for(int i=0;i<slots;i++){       
            res[z][i] = nr[z][i] * m1[z][i]-pp[z][i];
            res[z][i] = res[z][i] % p;
            if(res[z][i] >= p/2){
                res[z][i]= res[z][i] - p;
            }
            else if(res[z][i] <= -p/2){
                res[z][i] = res[z][i] + p;
            }
                cout<<res[z][i]<<" ";
        }
    }
    cout<<endl;
    cout<<q2<<endl;
    cout<<q<<endl;
    cout<<du<<endl;
    cout<<dv<<endl;
    cout<<p<<endl;  
   return 0;
    }