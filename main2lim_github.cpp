//
//  main.cpp
//  matValue
//
//  Created by Côme Huré on 14/12/15.
//  Copyright © 2015 Côme Huré. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <random>
#include <math.h>
#include <tuple>
#include <algorithm>
#include <functional>
#include <numeric>
#include <set>
#include <map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <flann/flann.h>



using namespace std;

template <typename T>
vector<T> operator+(const vector<T>& a, const vector<T>& b) {
    assert(a.size()==b.size());
    vector<T> result;
    result.resize(a.size());
    transform(a.begin(), a.end(), b.begin(), result.begin(), plus<T>());
    return result;
}
template <typename T>
vector<T> operator*(const T& a, const vector<T>& b) {
    vector<T> result;
    result.resize(b.size());
    for (int i =0; i<b.size(); i++) {
        result[i]=a*b[i];
    }
    return result;
}
template <typename T>
vector<T> operator-(const vector<T>& a, const vector<T>& b) {
    assert(a.size()==b.size());
    vector<T> result;
    result.reserve(a.size());
    transform(a.begin(), a.end(), b.begin(), back_inserter(result), [](T a, T b) {return a-b;});
    return result;
}

template <typename T>
bool comparetau(const vector<T> &a, const vector<T> &b) {
    assert(a.size() == b.size());
    if (a.back()< b.back()) {
        return 1;
    }
    else return 0;
}

vector<int> e(int j, int K) {
    vector<int> c(K+1,0);
    if (j<=K&&j>=0) {
        c[j]=1;
        return c;
    }
    else {return c;}
}

template <typename T>
ostream& operator<<(ostream& os, vector<T> const& a) {
    for (int i =0; i<a.size(); i++) {
        os << a[i] << " ";
    }
    return os;
}

double intL0(int n) {
    if (n==0) {
        return .25;
    }
    return .75;
}

double intC0(int n) {
    switch (n) {
        case 0:
            return 0;
            break;
        case 1:
            return 0.2;
            break;
        case 2:
            return .4;
            break;
        case 3:
            return .46;
            break;
        case 4:
            return .5;
            break;
        case 5:
            return .55;
            break;
        case 6:
            return .58;
            break;
        case 7:
            return .625;
            break;
        case 8:
            return .65;
            break;
        case 9:
            return .7;
            break;
        case 10:
            return .7;
            break;
        case 11:
            return .75;
            break;
        case 12:
            return .78;
            break;
        default:
            return .8;
            break;
    }
}

double intM(int a, int n) {
    if(a>0) {
    switch (n) {
        case 0:
            return 0;
            break;
        case 1:
            return .25;
            break;
        case 2:
            return .19;
            break;
        case 3:
            return .15;
            break;
        case 4:
            return .125;
            break;
        case 5:
            return .1;
            break;
        case 6:
            return .09;
            break;
        default:
            return .08;
            break;
    }
    }
    else {
        if(n==0) return 0;
        else if (n==1) return .03;
        else if (n==2) return .02;
        else if (n>=3 and n<=13) return .0125;
        else if (n>13 and n <20) return .015;
        else if (n>=20 and n <=31) return .02;
        else return .025;
    }
}

double intL1(int a, int n) {
    if (a>0) {
        if (n==0 or n==1) {
            return 1.1;
        }
        else return .9;
    }
    else {
        if (n==0) {
            return 1.8;
        }
        else if (n==1) return 1.1;
        else if (n==2) return .7;
        else if (n>=3 and n<=8) return .6;
        else if (n>8 and n <=16) return .4;
        else return .035;
    }
}

double intC1(int a, int n) {
    if (a>0) {
        if (n==0) {
            return 0;
        }
        else if (n>=2 and n <=4) return .3;
        else return .35;
    }
    else {
        if (n==0) {
            return 0;
        }
        else if (n==1) return .4;
        else if (n==2) return .6;
        else if (n>=3 and n<=6) return .8;
        else if (n>6 and n<=12) return 1;
        else return 1.15;
    }
}
int A0(vector<int> a,int K, int ainf, int i=0) {
    int n=0;
    int A=abs(a[0]);
    while (A<i+1) {
        n++;
        if (n<K+1) {A+=abs(a[n]);}
        else {A+=abs(ainf);}
    }
    return n;
}
int A1(vector<int> a,int K, int ainf, int i) {
    return A0(a,K,ainf,a[A0(a,K,ainf,0)]+i);
}


int B0(vector<int> b, int K, int binf, int i) {
    int n=0;
    int A=abs(b[0]);
    while (A<i+1) {
        n++;
        if (n<K+1) {A+=abs(b[n]);}
        else {A+=abs(binf);}
    }
    return n;
}

int B1(vector<int> b, int K, int binf, int i) {
    return B0( b, K, binf, abs(b[B0(b, K, binf, 0)])+i);
}

class Carnet {
public:
    //mt19937 eng {random_device{}()};
    static std::mt19937 eng;
    double Tps;
    int ainf; int binf;
    int K;
    int x; int y;
    vector<int> a;
    vector<int> b;
    int dp;
    int pa; int pb;
    int ra; int rb;
    int na; int nb;
    // Intensités:
    vector<double> intL;
    vector<double> intC;vector<double> intLa;vector<double> intCa;vector<double> intLb;vector<double> intCb;
    int sizeLa; int sizeLb; int sizeCa; int sizeCb;
    double intMp;
    double intMm;
    Carnet() : Tps(0), ainf(5), binf(-5), K(2), x(0), y(0), dp(1), pa(101), pb(100), intMp(4.), intMm(2.), na(1), nb(1), ra(0), rb(0) {
        a.resize(K+1);
        a[0]=1; a[1]=1; a[K]=ainf;
        b.resize(K+1);
        b[0]=-1; b[1]=-1; b[K]=binf;
//        intL.resize(K+1);
//        intC.resize(K+1);
//        for (int i =0; i<K; i++) {
//            intL[i]= (K-i);
//            intC[i] = .5*double(i+1)/K;
//        }
        
        intL.resize(K+1);
        intL[0]=3.5; intL[1]=2.; intL[2]=0.;
        intC.resize(K+1);
        intC[0]=0.1; intC[1]=0.3; intC[2]=0.;
        
        
        //intL[0]=intL0(a[0]); intL[1]=intL1(a[0], a[1]);
        
        //intC[0]=intC0(a[0]); intC[1]=intC1(a[0], a[1]); intC[2]=0;
    }
    double test() {
        uniform_real_distribution<>diss(1,2);
        return diss(eng);
    }
    
    void loadData(float* dataset, int i ,int cols) {
        x=(int)dataset[i*cols];
        y=(int)dataset[i*cols+1];
        for (int k=0; k<K; k++) {
            a[k]=(int)dataset[i*cols+2+k];
        }
        for (int k=0; k<K; k++) {
            b[k]=(int)dataset[i*cols+2+K+k];
        }
        na=(int)dataset[i*cols+2+K+K];
        nb=(int)dataset[i*cols+2+K+K+1];
        pa=(int)dataset[i*cols+2+K+K+2];
        pb=(int)dataset[i*cols+2+K+K+3];
        ra=(int)dataset[i*cols+2+K+K+4];
        rb=(int)dataset[i*cols+2+K+K+5];
    }
    void loadVector(vector<int> v) {
        x=(int)v[0];
        y=(int)v[1];
        for (int k=0; k<K; k++) {
            a[k]=(int)v[2+k];
        }
        for (int k=0; k<K; k++) {
            b[k]=(int)v[2+K+k];
        }
        na=(int)v[2+K+K];
        nb=(int)v[2+K+K+1];
        pa=(int)v[2+K+K+2];
        pb=(int)v[2+K+K+3];
        ra=(int)v[2+K+K+4];
        rb=(int)v[2+K+K+5];
    }
    
    Carnet& operator=(const Carnet& c ) {
        Tps=c.Tps;
        ainf=c.ainf; binf=c.binf; K=c.K; x=c.x; y=c.y; a=c.a; b=c.b; dp=c.dp; pa=c.pa; pb=c.pb; ra=c.ra; rb=c.rb; na=c.na; nb=c.nb; intL=c.intL; intC=c.intC; intLa=c.intLa; intLb=c.intLb; intCa=c.intCa;intCb=c.intCb; sizeLa=c.sizeLa; sizeLb=c.sizeLb; sizeCa=c.sizeCa; sizeCb=c.sizeCb; Tps=c.Tps;
        return *this;
    }
    Carnet(const Carnet& c) {
        Tps=c.Tps;
        ainf=c.ainf; binf=c.binf; K=c.K; x=c.x; y=c.y; a=c.a; b=c.b; dp=c.dp; pa=c.pa; pb=c.pb; ra=c.ra; rb=c.rb; na=c.na; nb=c.nb; intL=c.intL; intC=c.intC;
    }
    int A0(int i=0) {
        int n=0;
        int A=abs(a[0]);
        while (A<i+1) {
            n++;
            if (n<K+1) {A+=abs(a[n]);}
            else {A+=abs(ainf);}
        }
        return n;
    }
    int A1(int i) {
        return A0(a[A0(0)]+i);
    }
    
    int A2(int i) {
        return A0(a[A0(0)]+a[A1(0)]+i);
    }
    
    
    int B0(int i) {
        int n=0;
        int A=abs(b[0]);
        while (A<i+1) {
            n++;
            if (n<K+1) {A+=abs(b[n]);}
            else {A+=abs(binf);}
        }
        return n;
    }
    
    int B1(int i) {
        return B0(abs(b[B0(0)])+i);
    }
    
    int B2(int i) {
        return B0(abs(b[A0(0)])+abs(b[B1(0)])+i);
    }
    
    /*
    void intensites(int la, int lb) {
//        intL[0]=intL0(a[0]); intL[1]=intL1(a[0], a[1]);
//        intC[0]=intC0(a[0]); intC[1]=intC1(a[0], a[1]);
        int a0=A0(0);
        if ((la!=1)&&(lb!=1)) {
            intLa=intL;
            intCa.resize(intC.size());
            for (int i=0; i<intC.size(); i++) {
                intCa[i]=intC[i]*a[i];
            }
//            intCa[0]=intC0(a[0]);
//            intCa[1]=intC1(a[0], a[1]);
            intLb=intL;
            intCb.resize(intC.size());
            for (int i=0; i<intC.size(); i++) {
                intCb[i]=intC[i]*abs(b[i]);
            }
//            intCb[0]=intC0(abs(b[0]));
//            intCb[1]=intC1(abs(b[0]), abs(b[1]));
            sizeLa=K+1;
            sizeCa=K+1;
            sizeLb=K+1;
            sizeCb=K+1;
        }
        if (la!=1&&lb==1) {
            intCa.resize(K-a0+rb+1);
            for (int i =0; i<K-a0+rb; i++) {
                intCa[i]=intC[i]*a[a0-rb+i];
            }
            intCa[K-a0+rb]=0;
            //intCa[0]=intC0(a[a0-rb]);
            //intCa[1]=intC1(a[a0-rb], a[a0-rb+1]);
            intLa.resize(K-a0+rb+1);
            for (int i =0; i<K-a0+rb; i++) {
                intLa[i]=intL[i];
            }
            intLa[K-a0+rb]=0;
            intLb=intL;
            intCb.resize(intC.size());
            for (int i=0; i<intC.size(); i++) {
                intCb[i]=intC[i]*abs(b[i]);
            }
//            intCb[0]=intC0(abs(b[0]));
//            intCb[1]=intC1(abs(b[0]), abs(b[1]));
            sizeLa=K-a0+rb+1;
            sizeCa=K-a0+rb+1;
            sizeLb=K+1;
            sizeCb=K+1;
        }
        
        if (la==1&&lb!=1) {
            intCa.resize(intC.size());
            for (int i=0; i<intC.size(); i++) {
                intCa[i]=intC[i]*a[i];
            }
//            intCa[0]=intC0(a[0]);
//            intCa[1]=intC1(a[0], a[1]);
            intLa=intL;
            intCb.resize(K-a0+ra+1);
            for (int i =0; i<K-a0+ra; i++) {
                intCb[i]=intC[i]*abs(b[a0-ra+i]);
            }
            intCb[K-a0+ra]=0;
//            intCb[0]=intC0(abs(b[a0-ra]));
//            intCb[1]=intC1(abs(b[a0-ra]), abs(b[a0-ra+1]));
            intLb.resize(K-a0+ra+1);
            for (int i =0; i<K-a0+ra; i++) {
                intLb[i]=intL[i];
            }
            intLb[K-a0+ra]=0;
            sizeLa=K+1;
            sizeCa=K+1;
            sizeLb=K-a0+ra+1;
            sizeCb=K-a0+ra+1;
        }
        
        if (la==1&&lb==1) {
            intCa.resize(K-a0+rb+1);
            for (int i =0; i<K-a0+rb; i++) {
                intCa[i]=intC[i]*a[a0-rb+i];
            }
            intCa[K-a0+rb]=0;
//            intCa[0]=intC0(a[a0-rb]);
//            intCa[1]=intC1(a[a0-rb], a[a0-rb+1]);
            intLa.resize(K-a0+rb+1);
            for (int i =0; i<K-a0+rb; i++) {
                intLa[i]=intL[i];
            }
            intLa[K-a0+rb]=0;
            intCb.resize(K-a0+ra+1);
            for (int i =0; i<K-a0+ra; i++) {
                intCb[i]=intC[i]*abs(b[a0-ra+i]);
            }
            intCb[K-a0+ra]=0;
//            intCb[0]=intC0(abs(b[a0-ra]));
//            intCb[1]=intC1(abs(b[a0-ra]), abs(b[a0-ra+1]));
            intLb.resize(K-a0+ra+1);
            for (int i =0; i<K-a0+ra; i++) {
                intLb[i]=intL[i];
            }
            intLb[K-a0+ra]=0;
            sizeLa=K-a0+rb+1;
            sizeCa=K-a0+rb+1;
            sizeLb=K-a0+ra+1;
            sizeCb=K-a0+ra+1;
        }
    }
    */
    
    void intensites(int la, int lb) {
                intL[0]=intL0(a[0]); intL[1]=intL1(a[0], a[1]);
                intC[0]=intC0(a[0]); intC[1]=intC1(a[0], a[1]);
        int a0=A0(0);
        if ((la!=1)&&(lb!=1)) {
            intLa=intL;
            intCa.resize(intC.size());
            for (int i=0; i<intC.size(); i++) {
                intCa[i]=intC[i]*a[i];
            }
//            intCa[0]=intC0(a[0]);
//            intCa[1]=intC1(a[0], a[1]);
            intLb=intL;
            intCb.resize(intC.size());
            for (int i=0; i<intC.size(); i++) {
                intCb[i]=intC[i]*abs(b[i]);
            }
//            intCb[0]=intC0(abs(b[0]));
//            intCb[1]=intC1(abs(b[0]), abs(b[1]));
            sizeLa=K+1;
            sizeCa=K+1;
            sizeLb=K+1;
            sizeCb=K+1;
        }
        if (la!=1&&lb==1) {
            intCa.resize(K-a0+rb+1);
            for (int i =0; i<K-a0+rb; i++) {
                intCa[i]=intC[i]*a[a0-rb+i];
            }
            intCa[K-a0+rb]=0;
//            intCa[0]=intC0(a[a0-rb]);
//            intCa[1]=intC1(a[a0-rb], a[a0-rb+1]);
//            intCa[K-a0+rb]=0;
            intLa.resize(K-a0+rb+1);
            for (int i =0; i<K-a0+rb; i++) {
                intLa[i]=intL[i];
            }
            intLa[K-a0+rb]=0;
            intLb=intL;
            intCb.resize(intC.size());
            for (int i=0; i<intC.size(); i++) {
                intCb[i]=intC[i]*abs(b[i]);
            }
//            intCb[0]=intC0(abs(b[0]));
//            intCb[1]=intC1(abs(b[0]), abs(b[1]));
            sizeLa=K-a0+rb+1;
            sizeCa=K-a0+rb+1;
            sizeLb=K+1;
            sizeCb=K+1;
        }
        
        if (la==1&&lb!=1) {
            intCa.resize(intC.size());
            for (int i=0; i<intC.size(); i++) {
                intCa[i]=intC[i]*a[i];
            }
            intCa[0]=intC0(a[0]);
            intCa[1]=intC1(a[0], a[1]);
            intLa=intL;
            intCb.resize(K-a0+ra+1);
            for (int i =0; i<K-a0+ra; i++) {
                intCb[i]=intC[i]*abs(b[a0-ra+i]);
            }
            intCb[K-a0+ra]=0;
//            intCb[0]=intC0(abs(b[a0-ra]));
//            intCb[1]=intC1(abs(b[a0-ra]), abs(b[a0-ra+1]));
//            intCb[K-a0+ra]=0;
            intLb.resize(K-a0+ra+1);
            for (int i =0; i<K-a0+ra; i++) {
                intLb[i]=intL[i];
            }
            intLb[K-a0+ra]=0;
            sizeLa=K+1;
            sizeCa=K+1;
            sizeLb=K-a0+ra+1;
            sizeCb=K-a0+ra+1;
        }
        
        if (la==1&&lb==1) {
            intCa.resize(K-a0+rb+1);
            for (int i =0; i<K-a0+rb; i++) {
                intCa[i]=intC[i]*a[a0-rb+i];
            }
            intCa[K-a0+rb]=0;
//            intCa[0]=intC0(a[a0-rb]);
//            intCa[1]=intC1(a[a0-rb], a[a0-rb+1]);
//            intCa[K-a0+rb]=0;
            intLa.resize(K-a0+rb+1);
            for (int i =0; i<K-a0+rb; i++) {
                intLa[i]=intL[i];
            }
            intLa[K-a0+rb]=0;
            intCb.resize(K-a0+ra+1);
            for (int i =0; i<K-a0+ra; i++) {
                intCb[i]=intC[i]*abs(b[a0-ra+i]);
            }
            intCb[K-a0+ra]=0;
//            intCb[0]=intC0(abs(b[a0-ra]));
//            intCb[1]=intC1(abs(b[a0-ra]), abs(b[a0-ra+1]));
//            intCb[K-a0+ra]=0;
            intLb.resize(K-a0+ra+1);
            for (int i =0; i<K-a0+ra; i++) {
                intLb[i]=intL[i];
            }
            intLb[K-a0+ra]=0;
            sizeLa=K-a0+rb+1;
            sizeCa=K-a0+rb+1;
            sizeLb=K-a0+ra+1;
            sizeCb=K-a0+ra+1;
        }
    }
    
    void Jmp() {
        vector<int> c(K+1);
        int a1=A0(1);
        int a0=A0(0);
        for (int i = 0; i<a1-a0; i++) {
            c[i]=0;
        }
        for (int i=a1-a0; i<K; i++) {
            c[i]=b[i-(a1-a0)];
        }
        c[K]=binf;
        b=c;
    }
    
    void Jcp() {
        vector<int> c(K+1);
        int a1=A0(1);
        int a0=A0(0);
        for (int i = 0; i<a1-a0; i++) {
            c[i]=0;
        }
        for (int i=a1-a0; i<K; i++) {
            c[i]=b[i-(a1-a0)];
        }
        c[K]=binf;
        b=c;
    }
    
    void J00Lp(int j, int a0) {
        vector<int> c(K+1);
        if (j<=a0) {
            for (int i = 0; i<K+1+j-a0; i++) {
                c[i]=b[i+a0-j];
            }
            for (int i=K+1+j-a0; i<K+1; i++) {
                c[i]=binf;
            }
            b=c;
        }
    }
    
    void J10Lp( int j, int a0) {
        if (j<=a0) {
            vector<int> c(K+1);
            for (int i = 0; i<K+1+j-a0; i++) {
                c[i]=b[i+a0-j];
            }
            for (int i=K+1+j-a0; i<K+1; i++) {
                c[i]=binf;
            }
            b=c;
        }
    }
    
    void J01Lp(int j) {
        if (j<=rb) {
            vector<int> c(K+1);
            for (int i = 0; i<K+1+j-rb; i++) {
                c[i]=b[i+rb-j];
            }
            for (int i=K+1+j-rb; i<K+1; i++) {
                c[i]=binf;
            }
            b=c;
        }
    }
    
    void J11Lp( int j) {
        if (j<=rb) {
            vector<int> c(K+1);
            for (int i = 0; i<K+1+j-rb; i++) {
                c[i]=b[i+rb-j];
            }
            for (int i=K+1+j-rb; i<K+1; i++) {
                c[i]=binf;
            }
            b=c;
        }
    }
    
    void Jmm() {
        vector<int> c(K+1);
        int b1=B0(1);
        int b0=A0(0);
        for (int i = 0; i<b1-b0; i++) {
            c[i]=0;
        }
        for (int i=b1-b0; i<K; i++) {
            c[i]=a[i-(b1-b0)];
        }
        c[K]=ainf;
        a=c;
    }
    
    void Jcm() {
        vector<int> c(K+1);
        int b1=B0(1);
        int b0=A0(0);
        for (int i = 0; i<b1-b0; i++) {
            c[i]=0;
        }
        for (int i=b1-b0; i<K; i++) {
            c[i]=a[i-(b1-b0)];
        }
        c[K]=ainf;
        a=c;
    }
    
    void J00Lm(int j) {
        int a0=A0(0);
        if (j<=a0){
            vector<int> c(K+1);
            for (int i = 0; i<K+1+j-a0; i++) {
                c[i]=a[i+a0-j];
            }
            for (int i=K+1+j-a0; i<K+1; i++) {
                c[i]=ainf;
            }
            a=c;
        }
    }
    
    void J01Lm(int j) {
        int a0=A0(0);
        if (j<=a0) {
            vector<int> c(K+1);
            for (int i = 0; i<K+1+j-a0; i++) {
                c[i]=a[i+a0-j];
            }
            for (int i=K+1+j-a0; i<K+1; i++) {
                c[i]=ainf;
            }
            a=c;
        }
    }
    
    void J10Lm(int j) {
        if (j<=ra) {
            vector<int> c(K+1);
            for (int i = 0; i<K+1+j-ra; i++) {
                c[i]=a[i+ra-j];
            }
            for (int i=K+1+j-ra; i<K+1; i++) {
                c[i]=ainf;
            }
            a=c;
        }
    }
    
    void J11Lm(int j) {
        if (j<=ra) {
            vector<int> c(K+1);
            for (int i = 0; i<K+1+j-ra; i++) {
                c[i]=a[i+ra-j];
            }
            for (int i=K+1+j-ra; i<K+1; i++) {
                c[i]=ainf;
            }
            a=c;
        }
    }
    
    bool operator<(const Carnet& c) const {
        if (Tps <c.Tps) return true;
        else if (c.Tps<Tps) return false;
        return 0;
    }
    
    double liquideI() {
        double V=0;
        if(y<0) {
            int temp=a[0];
            int k=0;
            while (-y>temp) {
                k++;
                if (k<K+1) {
                    temp+=a[k];
                }
                else
                    temp+=ainf;
            }
            if (k==0) {
                return y*pa;
            }
            int a0=A0(0);
            for (int i=a0; i<min(k, K+1); i++) {
                V+=(pa+(i-a0)*dp)*a[i];
            }
            for (int i  =K+1; i<k; i++) {
                V+=ainf*(pa+(i-a0)*dp);
            }
            if (k<K+1) {
                V+=(pa+(k-a0)*dp)*(-y-temp+a[k]);
            }
            else
                V+=(pa+(k-a0)*dp)*(-y-(temp-ainf));
            return -V;
        }
        if(y>0) {
            int temp=-b[0];
            int k=0;
            while (y>temp) {
                k++;
                if (k<K+1) {
                    temp-=b[k];
                }
                else
                    temp+=-binf;
            }
            if (k==0) {
                return y*pb;
            }
            int b0=A0(0);
            for (int i=b0; i<min(k, K+1); i++) {
                V+=(pb-(i-b0)*dp)*-b[i];
            }
            for (int i  =K+1; i<k; i++) {
                V+=-binf*(pb-(i-b0)*dp);
            }
            if (k<K+1) {
                V+=(pb-(k-b0)*dp)*(y-temp-b[k]);
            }
            else
                V+=(pb-(k-b0)*dp)*(y-(temp+binf));
            return V;
        }
        return 0;
    }
    
    double PL() {
        double V=liquideI();
        return x+V;
    }
    
    bool operator==(const Carnet& c) const {
        return x==c.x && y==c.y && a==c.a && b==c.b && na==c.na && nb==c.nb && pa==c.pa && pb==c.pb;
    }
    
    //    double dist(const Carnet& c) {
    //        double res=pow(x-c.x,2)+ pow(y-c.y,2);
    //        for (int i=0; i<K; i++) {
    //            res+=pow(a[i]-c.a[i], 2);
    //        }
    //        for (int i=0; i<K; i++) {
    //            res+=pow(b[i]-c.b[i], 2);
    //        }
    //        res+=pow(na-c.na, 2)+pow(nb-c.nb, 2)+pow(ra-c.ra, 2)+pow(rb-c.rb, 2)+pow(pa-c.pa, 2)+pow(pb-c.pb, 2);
    //        return res;
    //    }
    
    void Jump(int la, int lb) {
        na = na + (la==0)*(-1-na) +(la==1)*((ra>A0(0))+(ra==-1))*(a[A0(0)]+1 -na) +(la==2)*(ra <= A0(0))*(a[A1(0)]+1 - na);
        nb = nb
        +(lb==0)*(-1-nb)
        +(lb==1)*((rb>A0(0))+(rb==-1))*(abs(b[A0(0)])+1 -nb)
        +(lb==2)*(rb <= A0(0))*(abs(b[B1(0)])+1 - nb);
        pa= pa + (la!=1)* ((0<=ra)*(ra <= A0(0))*(A0(0)-ra));
        pb= pb - (lb!=1)* ( (0<=rb)*(rb <= A0(0))*(A0(0)-rb));
        ra= ra + (la==0)*(-1-ra) +(la==1)*((ra>A0(0))+(ra==-1)) *(A0(0)-ra) +(la==2)*((ra>A1(0))+(ra <= A0(0))) *(A1(0)-ra);
        rb= rb + (lb==0)*(-1-rb) +(lb==1)*((rb>A0(0))+(rb==-1)) *(A0(0)-rb) +(lb==2)*((rb>B1(0))+(rb <= B0(0))) *(B1(0)-rb);
    }
    
    void dMp(int la, int lb) {
        int a0=A0(0);
        int a1=A0(1);
        x=x+(la==1)*pa*(na==1);
        y=y-(la==1)*(na==1);
        //int na_=na+la*(-1*(na>1)+(a[a0]+1-na)*(na==1))+(1-la)*(a[a0]*(a[a0]>1)+(a[a1]+1)*(a[a0]==1)-na);
        bool change = 0;
        if(lb!=0 && rb+a1-a0 >= K && !(ra<=a0&& na==1)) {change =1;}
        pa+=dp*(la==2)*(ra-a0)*(a[a0]==1)
        +dp*(la==1)*(a0-ra)
        +dp*(la==0)*(a1-a0);
        ra=ra+(la==1)*(na==1)*(a0-ra);
        rb=rb+(lb==1)*((la!=1)+(la==1)*(na>1))*(a1-a0)+(lb==2)*(((la!=1)+(la==1)*(na>1))*(a1-a0));
        if (rb>K) rb=K;
        if ((la!=1) or ((la==1)&&(na>1))) {
            Jmp();
        }
        int bb0=B1(0);
//        if(lb==2 && bb0==K) nb=binf -1;
        if(change) nb=-binf +1;
        int na_=na+(la==1)*(-1*(na>1)+(a[a0]+1-na)*(na==1));
        a=a-int((la!=1))*e(a0,K)-(la==1)*(na>1)*e(a0,K);
        na=na_;
    }
    
    void dLp(int la, int lb, int j) {
        int a0=A0(0);
        int aa0= A1(0);
        int b0=a0;
        na=na
        +(la==1)*(2-na)*((lb!=1)*(j<ra)+(lb==1)*(j<ra-(a0-rb)))
        +(la==2)*((2-na)*((lb==1)*(j>rb)*(j<rb-a0+ra)+(lb!=1)*(j>a0)*(j<ra))
                  +(a[a0]+1-na)*((lb==1)*(j<rb)+(lb!=1)*(j<a0))
                  );
        pa=pa
        +dp*(la==2)*((lb==1)*(-(rb-j)*(j<rb))+(lb!=1)*(-(a0-j)*(j<a0)))
        +dp*(la==1)*(-(lb==1)*(ra-(j+a0-rb))*(rb-(a0-ra)>j)-(lb!=1)*(ra-j)*(j<ra))
        +dp*(la==0)*((lb==1)*(-(rb-j)*(j<rb))+(lb!=1)*(-(a0-j)*(j<a0)))
        ;
        ra=ra
        + (la==1)*((lb!=1)*(j-ra)*(j<ra)+(lb==1)*(j+b0-rb-ra)*(j<rb-(b0-ra)))
        +(la==2)*((lb==1)*((a0-ra)*(j<rb)+((j+a0-rb)-ra)*(j>rb)*(j<ra-(a0-rb)))+(lb!=1)*((j-ra)*(j>a0)*(j<ra)+(j<a0)*(a0-ra)));
        int rb_=rb+ (lb==1)*(j-rb)*(j<rb)+(lb==2)*(j+rb-b0-rb)*(j<a0);
        a=a+int((lb!=1))*e(j,K)+int((lb==1))*e(j+b0-rb,K);
        if (lb!=1) {
            if ((la!=1)&&(j<=b0-1)) {
                J00Lp(j,a0);
            }
            if((la==1)&& (b0>j)) {
                J10Lp(j,a0);
            }
        }
        else {
            if ((la==1) && (rb>j)) {
                J11Lp(j);
            }
            if((la!=1)&&(rb-1>=j)) {
                J01Lp(j);
            }
        }
        rb=rb_;
//        if (all_of(a.begin(), --a.end(), [](int i) {return i==0;})) {
//            if (ra==K) {
//                na=ainf+1;
//            }
//            if (rb==K) {
//                nb=abs(binf)+1;
//            }
//            a[K]=ainf;
//            b[K]=binf;
//        }
    }
    
    void dCp(int la, int lb, const int& j) {
        int a0=A0(0);
        int aa0=A1(0);
        int aaa0=A2(0);
        int b0=a0;
        int a1=A0(1);
        bernoulli_distribution Xa0(float((na-1))/a[a0]);
        int pro0= int(Xa0(eng));
        bernoulli_distribution Xa1(float((na-1))/a[aa0]);
        int pro1= int(Xa1(eng));
        na = na
        +(la==1)*((-pro0)*(lb!=1)*(j==a0)-pro0*(lb==1)*(j==rb))
        +(la==2)*((lb==1)*-pro1*(j==rb+aa0-a0)+(lb!=1)*-pro1*(j==aa0)
        //+ (a[a0]==1)*((lb==1)*(j==rb)+(lb!=1)*(j==a0))*(a[aaa0]+1-na))
                  );
        pa = pa
        + dp*(la==2)*((lb==1)*((ra-a0)*(a[a0]==1)*(j==rb)) + (lb!=1)*((ra-a0)*(a[a0]==1)*(j==a0)))
        +dp*(la==0)*((lb==1)*((a1-a0)*(j==rb)) + (lb!=1)*((a1-a0)*(j==a0)));
        int rb_=rb + (lb==1)*(a1-a0)*(j==rb)+(lb==2)*(a1-a0)*(j==a0);
        if (rb_>K) rb_=K;
        if ((lb!=1)&& j==a0) {
            Jcp();
        }
        if ((lb==1) && j==rb) {
            Jcp();
        }
        a=a-int((lb!=1))*e(j,K)-int((lb==1))*e(j+b0-rb,K);
        rb=rb_;
        if (rb==K) nb=abs(binf)+1;
    }
    
    void dMm(const int la, const int lb) {
        int a0=A0(0);
        int b0=a0;
        int b1=B0(1);
        x+=-(lb==1)*pb*(nb==1);
        y+=(lb==1)*(nb==1);
        bool change=0;
        if(la!=0 && ra+b1-b0>=K && !(rb<=a0 && nb==1)) {change=1;};
        int nb_=nb+(lb==1)*(-1*(nb>1)+(abs(b[b0])+1-nb)*(nb==1));
        pb=pb-dp*(lb==1)*(b0-rb)-dp*(lb==2)*(rb-a0)*(abs(b[b0])==1)-(lb==0)*(b1-b0);
        ra=ra+(la==1)*((lb!=1)+(lb==1)*(nb>1))*(b1-b0)+(la==2)*(((lb!=1)+(lb==1)*(nb>1))*(b1-b0));
        if(ra >K) ra=K;
        rb=rb+(lb==1)*(nb==1)*(b0-rb);
        if ((lb!=1) or ((lb==1)&&(nb>1))) {
            Jmm();
        }
        int aa0=A1(0);
        if (change) na=ainf+1;
//        if(la==2&& aa0==K) na=ainf+1;
        vector<int> ee=e(b0,K);
        //cout << ee << " ; " << b << endl;
        b=b+int((lb!=1))*ee+(lb==1)*(nb>1)*ee;
        nb=nb_;
    }
    
    void dLm(int la, int lb, const int& j) {
        int a0=A0(0);
        int b0=a0;
        int bb0= B1(0);
        nb=nb
        +(lb==1)*(2-nb)*((la!=1)*(j<rb)+(la==1)*(j<rb-(b0-ra)))
        +(lb==2)*((2-nb)*((la==1)*(ra<j)*(j<ra+rb-b0)+(la!=1)*(j>a0)*(j<rb))
                  +(abs(b[b0])+1-nb)*((la==1)*(j<ra)+(la!=1)*(j<b0)));
        pb=pb
        -dp*(lb==2)*((la==1)*(-(ra-j)*(j<ra))+(la!=1)*(-(b0-j)*(j<b0)))
        -dp*(lb==1)*(-(la==1)*(rb-(j+b0-ra))*(ra-(b0-rb)>j)-(la!=1)*(rb-j)*(j<rb))
        -dp*(lb==0)*((la==1)*(-(ra-j)*(j<ra))+(la!=1)*(-(b0-j)*(j<b0)))
        ;
        int ra_=ra+ (la==1)*(j-ra)*(j<ra)+(la==2)*((j-a0)*(j<a0));
        rb=rb
        +(lb==1)*((la!=1)*(j-rb)*(j<rb)+(la==1)*(j+a0-ra-rb)*(j<ra-(a0-rb)))
        +(lb==2)*((la==1)*((b0-rb)*(j<ra)+(j>ra)*(j<rb-(b0-ra))*((j+b0-ra)-rb))
                  +(la!=1)*((j-rb)*(j<rb)*(j>b0)+(b0-rb)*(j<a0)));
        b=b-int((la!=1))*e(j,K)-int((la==1))*e(j+a0-ra,K);
        if (la!=1) {
            if ((lb!=1)&&(j<=a0-1)) {
                J00Lm(j);
            }
            if ((lb==1)&&(a0>j)) {
                J01Lm(j);
            }
        }
        if (la==1) {
            if (lb==1 && (ra>j)) {
                J11Lm(j);
            }
            if (lb!=1 && (ra-1>=j)) {
                J10Lm(j);
            }
        }
        ra=ra_;
    }
    
    void dCm(int la, int lb, const int& j) {
        int a0=A0(0);
        int b0=a0;
        int b1=B0(1);
        int bb0 = B1(0);
        int bbb0= B2(0);
        bernoulli_distribution Xb0(float((nb-1))/abs(b[b0]));
        bernoulli_distribution Xb1(float((nb-1))/abs(b[bb0]));
        int prob0 =Xb0(eng);
        int prob1 =Xb1(eng);
        nb=nb
        +(lb==1)*(-prob0*(la==1)*(j==ra)-prob0*(la!=1)*(j==a0))
        +(lb==2)*((la==1)*-prob1*(j==ra+bb0-b0) +(la!=1)*-prob1*(j==bb0)
                  //+(abs(b[b0])==1)*((la==1)*(j==ra)+(lb!=1)*(j==b0))*(abs(bbb0)+1-na)
                  );
        pb=pb
        - dp*(lb==2)*((la==1)*((rb-b0)*(abs(b[a0])==1)*(j==ra)) + (la!=1)*((rb-b0)*(abs(b[a0])==1)*(j==a0)))
        - dp*(lb==0)*((la==1)*(b1-b0)*(j==ra) + (la!=1)*((b1-a0)*(j==b0)));
        int ra_ = ra + (la==1)*(b1-b0)*(j==ra)+ (la==2)*(j==a0)*(b1-b0);
        if (ra_ > K) {ra_=K;}
        if (la!=1 &&(j==b0)) {
            Jcm();
        }
        if (la==1 &&(j==ra)) {
            Jcm();
        }
        b=b + int((la!=1))*e(j,K) + int((la==1))*e(j+a0-ra,K);
        ra=ra_;
        if(ra==K) {na=ainf+1;}
//        if (all_of(a.begin(), --a.end(), [](int i) {return i==0;})) {
//            if (ra==K) {
//                na=ainf+1;
//            }
//            if (rb==K) {
//                nb=abs(binf)+1;
//            }
//            a[K]=ainf;
//            b[K]=binf;
//        }
    }
    
    void toArray(float* tab) {
        tab[0]=x;
        tab[1]=y;
        for (int k=0; k<K; k++) {
            tab[2+k]= a[k];
        }
        for (int k=0; k<K; k++) {
            tab[2+K+k]=b[k];
        }
        tab[2+K+K]=na;
        tab[2+K+K+1]=nb;
        tab[2+K+K+2]=pa;
        tab[2+K+K+3]=pb;
        tab[2+K+K+4]=ra;
        tab[2+K+K+5]=rb;
    }
    vector<double> toVector() {
        vector<double> v(2+K+K+5+1+1);
        v[0]=x;
        v[1]=y;
        for (int k=0; k<K; k++) {
            v[2+k]= a[k];
        }
        for (int k=0; k<K; k++) {
            v[2+K+k]=b[k];
        }
        v[2+K+K]=na;
        v[2+K+K+1]=nb;
        v[2+K+K+2]=pa;
        v[2+K+K+3]=pb;
        v[2+K+K+4]=ra;
        v[2+K+K+5]=rb;
        v[2+K+K+6]=Tps;
        return v;
    }
};

static bool seeded = false;
ostream& operator<<(ostream& os, Carnet const& c) {
    os << c.x << " "<< c.y << " ";
    for (int i=0; i<c.K; i++) {
        os << c.a[i] << " ";
    }
    for (int i=0; i<c.K; i++) {
        os << c.b[i] << " ";
    }
    os << c.na << " " << c.nb << " " << c.pa << " " << c.pb << " " << c.ra << " " << c.rb ;
    return os;
}

class CarnetControle : public Carnet {
public:
    int la;
    int lb;
    CarnetControle(int la=1, int lb=1 ) : la(la), lb(lb) {}
    Carnet& operator=(const CarnetControle& c ) {
        Tps=c.Tps;
        ainf=c.ainf; binf=c.binf; K=c.K; x=c.x; y=c.y; a=c.a; b=c.b; dp=c.dp; pa=c.pa; pb=c.pb; ra=c.ra; rb=c.rb; na=c.na; nb=c.nb; intL=c.intL; intC=c.intC; intLa=c.intLa; intLb=c.intLb; intCa=c.intCa;intCb=c.intCb; sizeLa=c.sizeLa; sizeLb=c.sizeLb; sizeCa=c.sizeCa; sizeCb=c.sizeCb; Tps=c.Tps;
        la=c.la; lb=c.lb;
        return *this;
    }
    
    void dMp() {
        Carnet::dMp(la,lb);
    }
    
    void dLp(int j) {
        Carnet::dLp(la, lb, j);
    }
    
    void dCp(const int& j) {
        Carnet::dCp(la, lb, j);
    }
    
    void dMm() {
        Carnet::dMm(la, lb);
    }
    
    void dLm(const int& j) {
        Carnet::dLm(la, lb, j);
    }
    
    void dCm(const int& j) {
        Carnet::dCm(la, lb, j);
    }
    
    tuple<int,double> genEv() {
        int a0=A0(0);
        intensites(la, lb);
        double lambda=intMp+intMm+accumulate(intLa.begin(), intLa.end(), 0.)+accumulate(intLb.begin(), intLb.end(), 0.)+accumulate(intCa.begin(), intCa.end(), 0.)+accumulate(intCb.begin(), intCb.end(), 0.);
        exponential_distribution<double> distExp(lambda);
        double tau=distExp(eng);
        vector<double> prob= {intMp,intMm};
        prob.insert(prob.end(), intLa.begin(), intLa.end());
        prob.insert(prob.end(), intLb.begin(), intLb.end());
        prob.insert(prob.end(), intCa.begin(), intCa.end());
        prob.insert(prob.end(), intCb.begin(), intCb.end());
        prob=1./lambda*prob;
        discrete_distribution<int> dist_disc(prob.begin(),prob.end());
        int event=dist_disc(eng);
        return make_tuple(event,tau);
    }
    
    void move(int event, ostream & out) {
        //cout<< " event:" << event << endl;
        if (event==0) {
            //cout << "dMp :" << endl;
            out <<"dMp :" << flush;
            dMp();
        }
        if (event==1) {
            //cout << "dMm :" << endl;
            out <<"dMm :" << flush ;
            dMm();
        }
        for (int i=2; i<sizeLa+2; i++) {
            if (event ==i) {
                //cout << "dLp :" << i-2 << endl;
                out << "dLp :" << i-2 << flush;
                dLp(i-2);
            }
        }
        for (int i=sizeLa+2; i<sizeLa+2+sizeLb; i++) {
            if (event ==i) {
                //cout << "dLm :" << i-(2+sizeLa) << endl;
                out <<"dLm :" << i-(2+sizeLa) << flush ;
                dLm(i-(2+sizeLa));
            }
        }
        for (int i= sizeLa+2+sizeLb; i<sizeLa+2+sizeLb+sizeCa;i++) {
            if (event==i) {
                //cout << "dCp :" << i-(sizeLa+2+sizeLb) << endl;
                out << "dCp :" << i-(sizeLa+2+sizeLb) << flush ;
                dCp(i-(sizeLa+2+sizeLb));
            }
        }
        for (int i=sizeLa+2+sizeLb+sizeCa; i<sizeLa+2+sizeLb+sizeCa+sizeCb; i++) {
            if (event==i) {
                //cout << "dCm :" << i-(sizeLa+2+sizeLb+sizeCa) << endl;
                out << "dCm :" << i-(sizeLa+2+sizeLb+sizeCa) << flush ;
                dCm(i-(sizeLa+2+sizeLb+sizeCa));
            }
        }
    }
    
    void move(int event, vector<double>& record) {
        //cout<< " event:" << event << endl;
        if (event==0) {
            record[0]+=1;
            //cout << "dMp :" << endl;
            dMp();
        }
        if (event==1) {
            record[1]+=1;
            //cout << "dMm :" << endl;
            dMm();
        }
        for (int i=2; i<sizeLa+2; i++) {
            if (event ==i) {
                record[i]+=1;
                //cout << "dLp :" << i-2 << endl;
                dLp(i-2);
            }
        }
        for (int i=sizeLa+2; i<sizeLa+2+sizeLb; i++) {
            if (event ==i) {
                record[2+K +i-sizeLa -2]+=1;
                //cout << "dLm :" << i-(2+sizeLa) << endl;
                dLm(i-(2+sizeLa));
            }
        }
        for (int i= sizeLa+2+sizeLb; i<sizeLa+2+sizeLb+sizeCa;i++) {
            if (event==i) {
                //cout << "dCp :" << i-(sizeLa+2+sizeLb) << endl;
                dCp(i-(sizeLa+2+sizeLb));
            }
        }
        for (int i=sizeLa+2+sizeLb+sizeCa; i<sizeLa+2+sizeLb+sizeCa+sizeCb; i++) {
            if (event==i) {
                //cout << "dCm :" << i-(sizeLa+2+sizeLb+sizeCa) << endl;
                dCm(i-(sizeLa+2+sizeLb+sizeCa));
            }
        }
    }
    
    
    void move(int event, double& record) {
        //cout<< " event:" << event << endl;
        if (event==0) {
            record+=1;
            //cout << "dMp :" << endl;
            dMp();
        }
        if (event==1) {
            
            //cout << "dMm :" << endl;
            dMm();
        }
        for (int i=2; i<sizeLa+2; i++) {
            if (event ==i) {
                //cout << "dLp :" << i-2 << endl;
                dLp(i-2);
        }
    }
    for (int i=sizeLa+2; i<sizeLa+2+sizeLb; i++) {
        if (event ==i) {
            //cout << "dLm :" << i-(2+sizeLa) << endl;
            dLm(i-(2+sizeLa));
        }
    }
    for (int i= sizeLa+2+sizeLb; i<sizeLa+2+sizeLb+sizeCa;i++) {
        if (event==i) {
            //cout << "dCp :" << i-(sizeLa+2+sizeLb) << endl;
            dCp(i-(sizeLa+2+sizeLb));
        }
    }
    for (int i=sizeLa+2+sizeLb+sizeCa; i<sizeLa+2+sizeLb+sizeCa+sizeCb; i++) {
        if (event==i) {
            //cout << "dCm :" << i-(sizeLa+2+sizeLb+sizeCa) << endl;
            dCm(i-(sizeLa+2+sizeLb+sizeCa));
        }
    }
}


    void move(int event) {
        //cout<< " event:" << event << endl;
        if (event==0) {
            //cout << "dMp :" << endl;
            dMp();
        }
        if (event==1) {
            //cout << "dMm :" << endl;
            dMm();
        }
        for (int i=2; i<sizeLa+2; i++) {
            if (event ==i) {
                //cout << "dLp :" << i-2 << endl;
                dLp(i-2);
            }
        }
        for (int i=sizeLa+2; i<sizeLa+2+sizeLb; i++) {
            if (event ==i) {
                //cout << "dLm :" << i-(2+sizeLa) << endl;
                dLm(i-(2+sizeLa));
            }
        }
        for (int i= sizeLa+2+sizeLb; i<sizeLa+2+sizeLb+sizeCa;i++) {
            if (event==i) {
                //cout << "dCp :" << i-(sizeLa+2+sizeLb) << endl;
                dCp(i-(sizeLa+2+sizeLb));
            }
        }
        for (int i=sizeLa+2+sizeLb+sizeCa; i<sizeLa+2+sizeLb+sizeCa+sizeCb; i++) {
            if (event==i) {
                //cout << "dCm :" << i-(sizeLa+2+sizeLb+sizeCa) << endl;
                dCm(i-(sizeLa+2+sizeLb+sizeCa));
            }
        }
    }
};

ostream& operator<<(ostream& os, CarnetControle const& c) {
    os << c.x << " "<< c.y << " ";
    for (int i=0; i<c.K+1; i++) {
        os << c.a[i] << " ";
    }
    for (int i=0; i<c.K+1; i++) {
        os << c.b[i] << " ";
    }
    os << c.na << " " << c.nb << " " << c.pa << " " << c.pb << " " << c.ra << " " << c.rb << " "<< c.la << " " << c.lb << " " << c.Tps << endl;
    return os;
}

ostream& operator<<(ostream& os, vector<double> const& a) {
    for (int i =0; i<a.size(); i++) {
        os << a[i] << " ";
    }
    return os;
}

class Proc {
public:
    int N;
    double T;
    vector<CarnetControle> proc;
    
    Proc(int n=500, double t=10): N(n),T(t) {
        proc.resize(N+1);
        CarnetControle carnet;
        proc[0]=carnet;
    }
    
    void simule() {
        ofstream screen("screen", ios::out);
        CarnetControle carnet;
        uniform_int_distribution<int> Unif1(0,1);
        uniform_int_distribution<int> Unif2(0,2);
        for (int i =1; i<N+1; i++) {
            cout << "tour : " << i << endl;
            if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                carnet.la=Unif1(carnet.eng);
                carnet.lb=Unif1(carnet.eng);
            }
            else {
                
                carnet.la=Unif2(carnet.eng);
                carnet.lb=Unif2(carnet.eng);
            }
            cout << "avant: " << carnet << endl;
            screen <<"avant: " << carnet << flush;
            carnet.Jump(carnet.la,carnet.lb);
            if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                if (carnet.ra !=-1) {
                    carnet.pa+=carnet.K-carnet.ra;
                }
                carnet.ra=carnet.K;
                carnet.na=carnet.ainf+1;
                if (carnet.rb!=-1) {
                    carnet.pb-=carnet.K-carnet.rb;
                }
                carnet.rb=carnet.K;
                carnet.nb=abs(carnet.binf)+1;
                carnet.a[carnet.K]=carnet.ainf;
                carnet.b[carnet.K]=carnet.binf;
            }
            carnet.intensites(carnet.la, carnet.lb);
            tuple<int,double> pair=carnet.genEv();
            int event = get<0>(pair);
            double tau= get<1>(pair);
            if (carnet.Tps+tau <T) {
                carnet.Tps+=tau;
                cout << "apres jump : " << carnet << endl;
                screen <<"apres jump : " << carnet << flush;
                carnet.move(event, screen);
                int a0=carnet.A0(0);
                if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                    //cout << "modifie" << endl;
                    if (carnet.ra !=-1) {
                        carnet.pa+=carnet.K-carnet.ra;
                    }
                    carnet.ra=carnet.K;
                    carnet.na=carnet.ainf +1;
                    if (carnet.rb!=-1) {
                        carnet.pb-=carnet.K-carnet.rb;
                    }
                    carnet.rb=carnet.K;
                    carnet.nb=abs(carnet.binf)+1;
                    carnet.a[carnet.K]=carnet.ainf;
                    carnet.b[carnet.K]=carnet.binf;
                }
                int aa0=carnet.A1(0);
                int bb0=carnet.B1(0);
                a0=carnet.A0(0);
                int b0=carnet.B0(0);
                cout << "apres move " << carnet << endl;
                screen <<"apres move " << carnet << endl;
                //assert(carnet.pa-carnet.pb==carnet.ra+carnet.rb-a0+1);
                assert(carnet.ra<=carnet.K && carnet.rb <=carnet.K);
                assert(carnet.b[carnet.K]==-5 && carnet.a[carnet.K]==5);
                assert(all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i>=0;}));
                assert(all_of(carnet.b.begin(), --carnet.b.end(), [](int i) {return i<=0;}));
                assert(carnet.pa>carnet.pb);
                assert(carnet.pa-carnet.pb<=carnet.A0(0)+1);
                if(carnet.la!=0 && carnet.lb !=0) assert(carnet.pa-carnet.pb==min(carnet.A0(0),carnet.rb)+min(carnet.A0(0),carnet.ra)-carnet.A0(0)+1);
                assert(carnet.pa-carnet.pb<=2*carnet.K+1);
                assert(carnet.A0(0)==carnet.B0(0));
                if (carnet.la==1) assert(carnet.na <= carnet.a[a0]+1);
                if (carnet.la==2) assert(carnet.na <= max(carnet.a[aa0], carnet.a[a0])+1);
                if (carnet.lb==1) assert(carnet.nb<=abs(carnet.b[a0])+1);
                if (carnet.lb==2)assert(carnet.nb<=max(abs(carnet.b[bb0]),abs(carnet.b[b0]))+1);
                proc[i]=carnet;
            }
            else {
                carnet.Tps=T;
                proc[i]=carnet;
            }
        }
    }
};

class DataPrc {
public:
    vector<vector<vector<double>>> points;
    int N;
    double T;
    int NbSimul;
    DataPrc(int nbSimul=10000, int n=50, double tps=1) : NbSimul(nbSimul), N(n), T(tps) { // nbSimul = NbPoints ; N = nb de subdivisions
        CarnetControle carnet;
        points.resize(N+1);
        for (int i=0; i<n; i++) {
            points[i].reserve(nbSimul);
        }
        for (int i=0; i<nbSimul; i++) {
            points[0].push_back(carnet.toVector());
        }
    }
    
    void simule() {
        uniform_int_distribution<int> Unif1(0,1);
        uniform_int_distribution<int> Unif2(0,2);
        CarnetControle carnet;
        CarnetControle c_init;
        for (int j =0 ; j<NbSimul; j++) {
            carnet=c_init;
            for (int i =1; i<N+1; i++) {
                cout << "tour : " << i << endl;
                if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                    carnet.la=Unif1(carnet.eng);
                    carnet.lb=Unif1(carnet.eng);
                }
                else {
                    carnet.la=Unif2(carnet.eng);
                    carnet.lb=Unif2(carnet.eng);
                }
                cout << "avant: " << carnet << endl;
                //carnet.intensites(carnet.la, carnet.lb);
                carnet.Jump(carnet.la,carnet.lb);
                if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                    //cout << "modifie" << endl;
                    if (carnet.ra !=-1) {
                        carnet.pa+=carnet.K-carnet.ra;
                    }
                    carnet.ra=carnet.K;
                    carnet.na=carnet.ainf+1;
                    if (carnet.rb!=-1) {
                        carnet.pb-=carnet.K-carnet.rb;
                    }
                    carnet.rb=carnet.K;
                    carnet.nb=abs(carnet.binf)+1;
                    carnet.a[carnet.K]=carnet.ainf;
                    carnet.b[carnet.K]=carnet.binf;
                }
                carnet.intensites(carnet.la, carnet.lb);
                //cout << "tour :" << i << " temps: "<< i*float(T)/N << endl;
                tuple<int,double> pair=carnet.genEv();
                int event = get<0>(pair);
                double tau= get<1>(pair);
                if (carnet.Tps+tau <T) {
                    carnet.Tps+=tau;
                    cout << "apres jump : " << carnet << endl;
                    carnet.move(event);
                    int a0=carnet.A0(0);
                    if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                        //cout << "modifie" << endl;
                        if (carnet.ra !=-1) {
                            carnet.pa+=carnet.K-carnet.ra;
                        }
                        carnet.ra=carnet.K;
                        carnet.na=carnet.ainf +1;
                        if (carnet.rb!=-1) {
                            carnet.pb-=carnet.K-carnet.rb;
                        }
                        carnet.rb=carnet.K;
                        carnet.nb=abs(carnet.binf)+1;
                        carnet.a[carnet.K]=carnet.ainf;
                        carnet.b[carnet.K]=carnet.binf;
                    }
                    int aa0=carnet.A1(0);
                    int bb0=carnet.B1(0);
                    a0=carnet.A0(0);
                    int b0=carnet.B0(0);
                    cout << "apres move " << carnet << endl;
                    //assert(carnet.pa-carnet.pb==carnet.ra+carnet.rb-a0+1);
                    assert(carnet.b[carnet.K]==-5 && carnet.a[carnet.K]==5);
                    assert(all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i>=0;}));
                    assert(all_of(carnet.b.begin(), --carnet.b.end(), [](int i) {return i<=0;}));
                    assert(carnet.pa>carnet.pb);
                    assert(carnet.pa-carnet.pb<=carnet.A0(0)+1);
                    if(carnet.la!=0 && carnet.lb !=0) assert(carnet.pa-carnet.pb==min(carnet.A0(0),carnet.rb)+min(carnet.A0(0),carnet.ra)-carnet.A0(0)+1);
                    assert(carnet.pa-carnet.pb<=2*carnet.K+1);
                    assert(carnet.A0(0)==carnet.B0(0));
                    if (carnet.la==1) assert(carnet.na <= carnet.a[a0]+1);
                    if (carnet.la==2) assert(carnet.na <= max(carnet.a[aa0], carnet.a[a0])+1);
                    if (carnet.lb==1) assert(carnet.nb<=abs(carnet.b[a0])+1);
                    if (carnet.lb==2)assert(carnet.nb<=max(abs(carnet.b[bb0]),abs(carnet.b[b0]))+1);
                    points[i].push_back(carnet.toVector());
                }
                else {
                    carnet.Tps=T;
                    points[i].push_back(carnet.toVector());
                }
            }
            
            
            
//            for (int i =1; i<N+1; i++) {
//                bernoulli_distribution B(0.5);
//                carnet.la=B(carnet.eng);
//                carnet.lb=B(carnet.eng);
//                //cout << "tour :" << i << " temps: "<< i*float(T)/N << endl;
//                tuple<int,double> pair=carnet.genEv();
//                int event = get<0>(pair);
//                double tau= get<1>(pair);
//                if (carnet.Tps+tau <T) {
//                    carnet.Tps+=tau;
//                    carnet.move(event);
//                    int a0=carnet.A0(0);
//                    if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
//                        //cout << "modifie" << endl;
//                        if (carnet.ra==carnet.K) {
//                            carnet.na=carnet.ainf+1;
//                        }
//                        if (carnet.rb==carnet.K) {
//                            carnet.nb=abs(carnet.binf)+1;
//                        }
//                        carnet.a[carnet.K]=carnet.ainf;
//                        carnet.b[carnet.K]=carnet.binf;
//                    }
//                    assert(carnet.pa-carnet.pb==carnet.ra+carnet.rb-a0+1);
//                    assert(carnet.na <= carnet.a[a0]+1);
//                    assert(carnet.nb<=abs(carnet.b[a0])+1);
//                    cout << carnet.toVector();
//                    points[i].push_back(carnet.toVector());
//                }
//                else {
//                    carnet.Tps=T;
//                    points[i].push_back(carnet.toVector());
//                }
//            }
        }
    }
};

/*


class DataPrc11 {
public:
    vector<vector<vector<double>>> points;
    
    int N;
    double T;
    int NbSimul;
    double h;
    DataPrc11(int nbSimul=10000, int n=50, double tps=1) : NbSimul(nbSimul), N(n), T(tps) { // nbSimul = NbPoints ; N = nb de subdivisions
        CarnetControle carnet;
        points.resize(N+1);
        for (int i=0; i<n; i++) {
            points[i].reserve(nbSimul);
        }
        for (int i=0; i<nbSimul; i++) {
            points[0].push_back(carnet.toVector());
        }
        h=T/N;
    }
    
    void simule() {
        CarnetControle carnet;
        carnet.la=1; carnet.lb=1;
        CarnetControle c_init;
        c_init.la=1; c_init.lb=1;
        for (int j =0 ; j<NbSimul; j++) {
            carnet=c_init;
            for (int i =1; i<N+1; i++) {
                //cout << "tour :" << i << " temps: "<< i*float(T)/N << endl;
                tuple<int,double> pair=carnet.genEv();
                int event = get<0>(pair);
                double tau= get<1>(pair);
                carnet.Tps+=tau;
                while (carnet.Tps < i*float(T)/N) {
                    carnet.move(event);
                    int a0=carnet.A0();
                    if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                        if (carnet.ra==carnet.K) {
                            carnet.na=carnet.ainf+1;
                        }
                        if (carnet.rb==carnet.K) {
                            carnet.nb=abs(carnet.binf)+1;
                        }
                        carnet.a[carnet.K]=carnet.ainf;
                        carnet.b[carnet.K]=carnet.binf;
                    }
                    assert(carnet.pa-carnet.pb==carnet.ra+carnet.rb-a0+1);
                    assert(carnet.na <= carnet.a[a0]+1);
                    assert(carnet.nb<=abs(carnet.b[a0])+1);
                    tuple<int,double> pair=carnet.genEv();
                    event = get<0>(pair);
                    tau= get<1>(pair);
                    carnet.Tps+=tau;
                }
                carnet.Tps=i*float(T)/N;
                //                bernoulli_distribution B(0.5);
                //                carnet.la=B(carnet.eng);
                //                carnet.lb=B(carnet.eng);
                //cout<< carnet;
                points[i].push_back(carnet.toVector());
            }
        }
    }
};
 
*/


pair<int, int> argmax(vector<vector<double>>& temp, bool vide=0) {
    pair<int, int> res;
    res=make_pair(0, 0);
    double max=temp[0][0];
    vector<int> v;
    if(vide) {v={0,1};}
    else {v={0,1,2};
        //v={0,1};
    }
    for(int la:v) {
        for(int lb: v) {
            if (temp[la][lb]>max) {
                res.first=la;
                res.second=lb;
                max=temp[la][lb];
            }
        }
    }
    //cout << "max en " << res.first << " " << res.second << endl;
    return res;
}

void writeset() {
    
}

class Data {
public:
    vector<vector<vector<double>>> points;
    
    int N;
    double T;
    int NbSimul;
    double h;
    Data(int nbSimul=10000, int n=50, double tps=1) : NbSimul(nbSimul), N(n), T(tps) { // nbSimul = NbPoints ; N = nb de subdivisions
        CarnetControle carnet;
        points.resize(N+1);
        points[0].push_back(carnet.toVector());
        for (int t=1; t<N+1; t++) {
            points[t].reserve(NbSimul);
        }
    }
    
    void simule() {
        vector<double> record;
        CarnetControle carnet;
        CarnetControle c_init;
        uniform_int_distribution<int> Unif1(0,1);
        uniform_int_distribution<int> Unif2(0,2);
        uniform_int_distribution<int> UnifRemplace(0,1);
        bernoulli_distribution B(.5);
        
        record.resize(2+carnet.K+carnet.K);
        
        for (int j =0 ; j<NbSimul; j++) {
            carnet=c_init;
            for (int i =1; i<N+1; i++) {
                //cout << "tour : " << i << endl;
            
                //cout << "avant: " << carnet << endl;
                //carnet.intensites(carnet.la, carnet.lb);
            
                carnet.Jump(carnet.la,carnet.lb);
                if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                    //cout << "modifie" << endl;
                    if (carnet.ra !=-1) {
                        carnet.pa+=carnet.K-carnet.ra;
                    }
                    carnet.ra=carnet.K;
                    carnet.na=carnet.ainf+1;
                    if (carnet.rb!=-1) {
                        carnet.pb-=carnet.K-carnet.rb;
                    }
                    carnet.rb=carnet.K;
                    carnet.nb=abs(carnet.binf)+1;
                    carnet.a[carnet.K]=carnet.ainf;
                    carnet.b[carnet.K]=carnet.binf;
                }
                carnet.intensites(carnet.la, carnet.lb);
                //cout << "tour :" << i << " temps: "<< i*float(T)/N << endl;
                tuple<int,double> pair=carnet.genEv();
                int event = get<0>(pair);
                double tau= get<1>(pair);
                if (carnet.Tps+tau <T) {
                    carnet.Tps+=tau;
                    //cout << "apres jump : " << carnet << endl;
                    carnet.move(event, record);
                    int a0=carnet.A0(0);
                    if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                        //cout << "modifie" << endl;
                        if (carnet.ra !=-1) {
                            carnet.pa+=carnet.K-carnet.ra;
                        }
                        carnet.ra=carnet.K;
                        carnet.na=carnet.ainf +1;
                        if (carnet.rb!=-1) {
                            carnet.pb-=carnet.K-carnet.rb;
                        }
                        carnet.rb=carnet.K;
                        carnet.nb=abs(carnet.binf)+1;
                        carnet.a[carnet.K]=carnet.ainf;
                        carnet.b[carnet.K]=carnet.binf;
                    }
                    int aa0=carnet.A1(0);
                    int bb0=carnet.B1(0);
                    a0=carnet.A0(0);
                    int b0=carnet.B0(0);
                    //cout << "apres move " << carnet << endl;
                    //assert(carnet.pa-carnet.pb==carnet.ra+carnet.rb-a0+1);
                    assert(carnet.b[carnet.K]==-5 && carnet.a[carnet.K]==5);
                    assert(all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i>=0;}));
                    assert(all_of(carnet.b.begin(), --carnet.b.end(), [](int i) {return i<=0;}));
                    assert(carnet.pa>carnet.pb);
                    assert(carnet.pa-carnet.pb<=carnet.A0(0)+1);
                    if(carnet.la!=0 && carnet.lb !=0) assert(carnet.pa-carnet.pb==min(carnet.A0(0),carnet.rb)+min(carnet.A0(0),carnet.ra)-carnet.A0(0)+1);
                    assert(carnet.pa-carnet.pb<=2*carnet.K+1);
                    assert(carnet.A0(0)==carnet.B0(0));
                    if (carnet.la==1) assert(carnet.na <= carnet.a[a0]+1);
                    if (carnet.la==2) assert(carnet.na <= max(carnet.a[aa0], carnet.a[a0])+1);
                    if (carnet.lb==1) assert(carnet.nb<=abs(carnet.b[a0])+1);
                    if (carnet.lb==2)assert(carnet.nb<=max(abs(carnet.b[bb0]),abs(carnet.b[b0]))+1);
                    points[i].push_back(carnet.toVector());
                }
                else {
                    carnet.Tps=T;
                    points[i].push_back(carnet.toVector());
                }
                
                if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                    carnet.la=Unif1(carnet.eng);
                    carnet.lb=Unif1(carnet.eng);
                }
                else {
                    carnet.la=Unif2(carnet.eng);
                    carnet.lb=Unif2(carnet.eng);
                }
//                carnet.la=B(carnet.eng);
//                carnet.lb=B(carnet.eng);
            }
        }
        for (int i=0; i<N+1; i++) {
            //cout << points[i].size() << endl;
            sort(points[i].begin(),points[i].end(),comparetau<double>);
        }
        cout << "dMp moyen " << 1./NbSimul*record;
     }
};


void newP(vector<double>& c, Carnet& cc ) {
    vector<int> a;
    vector<int> b;
    a.reserve(cc.K);
    b.reserve(cc.K);
    for (int i =2; i<2+cc.K; i++) {
        a.push_back(c[i]);
        b.push_back(c[i+cc.K]);
    }
    //cout << a << endl;
    //cout << b << endl;
    int na=c[2+cc.K+cc.K];
    int nb= c[2+2*cc.K+1];
    int ra=c[2+2*cc.K+2+2];
    int rb=c[2+2*cc.K+2+3];
    
    if(na==-1 and nb>=0) {
        int a0=A0(a,cc.K,cc.ainf,0);
        c[2+cc.K+cc.K]=a[a0]+1;
        c[2+2*cc.K+2+2]=a0;
    }
    else if (na>=0 and nb ==-1) {
        int b0=B0(b,cc.K,cc.binf,0);
        c[2+2*cc.K+1]=abs(b[b0])+1;
        c[2+2*cc.K+2+3]=b0;
    }
    else if (na==-1 and nb==-1) {
        int a0=A0(a,cc.K,cc.ainf,0);
        c[2+cc.K+cc.K]=a[a0]+1;
        c[2+2*cc.K+2+2]=a0;
        int b0=B0(b,cc.K,cc.binf,0);
        c[2+2*cc.K+1]=abs(b[b0])+1;
        c[2+2*cc.K+2+3]=b0;
    }
}

void newP(Carnet& c ) {
    if(c.na==-1 and c.nb>=0) {
        int a0=c.A0(0);
        c.na=c.a[a0]+1;
        c.ra=a0;
    }
    else if (c.na>=0 and c.nb ==-1) {
        int b0=c.B0(0);
        c.nb=abs(c.b[b0])+1;
        c.rb=b0;
    }
    else if (c.na==-1 and c.nb==-1) {
        int a0=c.A0(0);
        c.na=c.a[a0]+1;
        c.ra=a0;
        int b0=c.B0(0);
        c.nb=abs(c.b[b0])+1;
        c.rb=b0;
    }
}


void newP(float* c, Carnet& cc ) {
    vector<int> a;
    vector<int> b;
    a.reserve(cc.K);
    b.reserve(cc.K);
    for (int i =2; i<2+cc.K; i++) {
        a.push_back(c[i]);
        b.push_back(c[i+cc.K]);
    }
    //cout << a << endl;
    //cout << b << endl;
    int na=c[2+cc.K+cc.K];
    int nb= c[2+2*cc.K+1];
    int ra=c[2+2*cc.K+2+2];
    int rb=c[2+2*cc.K+2+3];
    
    if(na==-1 and nb>=0) {
        int a0=A0(a,cc.K,cc.ainf,0);
        c[2+cc.K+cc.K]=a[a0]+1;
        c[2+2*cc.K+2+2]=a0;
    }
    else if (na>=0 and nb ==-1) {
        int b0=B0(b,cc.K,cc.binf,0);
        c[2+2*cc.K+1]=abs(b[b0])+1;
        c[2+2*cc.K+2+3]=b0;
    }
    else if (na==-1 and nb==-1) {
        int a0=A0(a,cc.K,cc.ainf,0);
        c[2+cc.K+cc.K]=a[a0]+1;
        c[2+2*cc.K+2+2]=a0;
        int b0=B0(b,cc.K,cc.binf,0);
        c[2+2*cc.K+1]=abs(b[b0])+1;
        c[2+2*cc.K+2+3]=b0;
    }
}

void backP(float* c, int la, int lb) {
    Carnet cc;
    if(la==0 and lb==0) {
        c[2+cc.K+cc.K]=-1;
        c[2+2*cc.K+1]=-1;
        c[2+2*cc.K+2+2]=-1;
        c[2+2*cc.K+2+3]=-1;
    }
    else if (la!=0 and lb==0) {
        c[2+2*cc.K+1]=-1;
        c[2+2*cc.K+2+3]=-1;
    }
    else if (la==0 and lb!=0) {
        c[2+cc.K+cc.K]=-1;
        c[2+2*cc.K+2+2]=-1;
    }
}

void backP(Carnet& c, int la, int lb) {
    Carnet cc;
    if(la==0 and lb==0) {
        c.na=-1;
        c.nb=-1;
        c.ra=-1;
        c.rb=-1;
    }
    else if (la!=0 and lb==0) {
        c.nb=-1;
        c.rb=-1;
    }
    else if (la==0 and lb!=0) {
        c.na=-1;
        c.ra=-1;
    }
}


void genereData(int nbSim=10000, int subdiv=500, double tps=10, int parts=1) {
    Carnet cc;
    auto t1=chrono::high_resolution_clock::now();
    Data data(nbSim, subdiv,tps);
    data.simule();
    auto t2=chrono::high_resolution_clock::now();
    cout << " time : "<< std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " milliseconds\n";
    ofstream fout("set_tps"+to_string(0)+"nbSim"+to_string(nbSim)+"subd"+to_string(subdiv)+".dat");
    if (fout.is_open()) {
        fout << 1 << endl;
        data.points[0][0].pop_back();// Enleve la composante temps
        fout << data.points[0][0];
    }
    ofstream f2out("QEset_tps"+to_string(0)+"nbSim"+to_string(nbSim)+"subd"+to_string(subdiv)+".dat");
    f2out << 0;
    for (int t =1; t<subdiv+1; t++) {
        ofstream fout("QEset_tps"+to_string(t)+"nbSim"+to_string(nbSim)+"subd"+to_string(subdiv)+".dat");
        for (int i =0; i<parts; i++) {
            if (fout.is_open()) {
                //fout << (data.points[t][i*nbSim/parts].back() + data.points[t][(i+1)*nbSim/parts-1].back())/2. << endl;
                fout << data.points[t][(i+1)*nbSim/parts-1].back() << endl;
            }
        }
        ofstream f2out("set_tps"+to_string(t)+"nbSim"+to_string(nbSim)+"subd"+to_string(subdiv)+".dat");
        //ofstream f3out("setMod_tps"+to_string(t)+"nbSim"+to_string(nbSim)+"subd"+to_string(subdiv)+".dat");
        set<vector<double>> temp;
        set<vector<double>> temp2;
        for(int p =0; p < data.points[t].size();p++) {
            data.points[t][p].pop_back();// Enleve la composante temps
            newP(data.points[t][p],cc);
            temp.insert(data.points[t][p]);
//            newP(data.points[t][p],cc);
//            temp2.insert(data.points[t][p]);
        }
        f2out << temp.size() << endl;
        //f3out << temp.size() << endl;
        for(auto& p : temp) {
            f2out << p << endl;
            //f3out << p << endl;
        }
    }
}
//
//void testProc() {
//    Proc prc(500,10);
//    vector<vector<CarnetControle>> data;
//    data.reserve(2000);
//    auto t1=chrono::high_resolution_clock::now();
//    for (int i=0; i<2000; i++) {
//        prc.simule();
//        data.push_back(prc.proc);
//        prc = Proc();
//    }
//    auto t2=chrono::high_resolution_clock::now();
//    cout << " time : "<< std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " milliseconds\n";
//}

float* read_points(const string filename, int cols, int& taille)
{
    float* data;
    float* p;
    int i,j;
    
    ifstream fin(filename);
    if (!fin) {
        cout << "Cannot open input file : "+filename+".\n";
        exit(1);
    }
    
    int rows;
    fin >> rows;
    taille = rows;
    
    data = new float[rows*cols];
    if (!data) {
        cout << "Cannot allocate memory." << endl;
        exit(1);
    }
    p = data;
    
    for (i=0;i<rows;++i) {
        for (j=0;j<cols;++j) {
            fin >> *p;
            p++;
        }
    }
    return data;
}

void read_temps(const string filename, vector<double>& tps, int parts)
{
    float p;
    int i;
    
    ifstream fin(filename);
    if (!fin) {
        cout << "Cannot open input file : "+filename+".\n";
        exit(1);
    }
    for (i=0;i<parts;++i) {
        fin >> p;
        tps[i]=p;
    }
}

void inventaire(int t) {
    float* data;
    int rows;
    double tm;
    data=read_points("11set_tps"+to_string(t)+".dat", 12, rows);
    ofstream fout("inventaire_tps"+to_string(t)+".dat");
    if (fout.is_open()) {
        for (int i=0; i<rows; i++) {
            fout << data[i*12+1] << endl;
        }
    }
}

void PL(int t) {
    float* data;
    int rows;
    data=read_points("11set_tps"+to_string(t)+".dat", 12, rows);
    ofstream fout("PL_tps"+to_string(t)+".dat");
    if (fout.is_open()) {
        for (int i=0; i<rows; i++) {
            Carnet c;
            c.loadData(data,i,12);
            fout << c.PL() << endl;
        }
    }
}

void write_results(const string filename, int *data, int rows, int cols)
{
    int* p;
    int i,j;
    
    ofstream fout(filename);
    if (!fout) {
        cout << "Cannot open output file." << endl;
        exit(1);
    }
    
    p = data;
    for (i=0;i<rows;++i) {
        for (j=0;j<cols;++j) {
            fout << *p;
            p++;
        }
        fout << endl;
    }
}

void write_vector(const string filename, vector<float> res)
{
    ofstream fout(filename);
    if (!fout) {
        cout << "Cannot open output file." << endl;
        exit(1);
    }
    for (int i=0; i<res.size(); i++) {
        fout << res[i] << "\n";
    }
}
void write_array(const string filename, float* res, int rows)
{
    ofstream fout(filename);
    if (!fout) {
        cout << "Cannot open output file." << endl;
        exit(1);
    }
    for (int i=0; i<rows; i++) {
        fout << res[i] << "\n";
    }
}

int projeteT(double& t, vector<double>& t2) {
    int n=0;
    int taille =t2.size();
    while (t>t2[n] and n < taille-1) {
        n++;
    }
    return n;
}



using namespace flann;

class Quantif3 {
public:
    vector<vector<vector<vector<map<int,double>>>>> matProb;
    vector<vector<vector<double>>> Valeurs;
    vector<vector<vector<pair<int, int>>>> StratOpt;
    
    int Parts;
    int N;
    double T;
    int NbSimul;
    Quantif3(int nbSimul=10000, int n=50, double tps=1, int parts=1) : Parts(parts), NbSimul(nbSimul), N(n), T(tps) { // nbSimul = NbPoints ; N = nb de subdivisions
        Valeurs.resize(n+1);
        Valeurs[0].resize(1);
        Valeurs[0][0].resize(1);
        StratOpt.resize(n);
        matProb.resize(3);
        for (int i=0; i<3; i++) {
            matProb[i].resize(3);
            for (int j=0; j<3; j++) {
                matProb[i][j].resize(N);
                matProb[i][j][0].resize(1);
            }
        }
    }

    
    void fillmat(int t) {
        ofstream value("Value"+to_string(t)+".dat");
        int parts;
        if (t==0) {
            parts=1;
        }
        else {parts =Parts;};
        Carnet c;
        Carnet c_;
        float* testset = new float[12];
        int tcount = 1;
        int nn=1;
        int* result = new int[tcount*nn];
        float* dists =  new float[tcount*nn];
        struct FLANNParameters p;
        float speedup;
        flann_index_t index_id;
        p = DEFAULT_FLANN_PARAMETERS;
        p.cores=1;
        p.algorithm = FLANN_INDEX_KDTREE;
        p.trees = 8;
        p.log_level = FLANN_LOG_INFO;
        p.checks = 64;
        //        p.algorithm = FLANN_INDEX_AUTOTUNED;
        //        p.target_precision=0.95;
        float* set1;
        //float* set1Mod;
        float* set2;
        int rows1;
        int row1Mod;
        vector<double> t1(parts);
        int rows2;
        vector<double> t2(Parts);
        
        int cols=12;
        int ind=0;
        set1=read_points("set_tps"+to_string(t)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat", cols,rows1);
        //set1Mod=read_points("setMod_tps"+to_string(t)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat", cols,rows1Mod);
        read_temps("QEset_tps"+to_string(t)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat", t1, parts);
        //cout << t1 << endl;
        //        for (int i=0; i<parts; i++) {
        //            Carnet cc;
        //            cc.loadData(set1[i],0,cols);
        //            cout << cc << " ; " ;
        //        }
        //        cout << endl;
        set2=read_points("set_tps"+to_string(t+1)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat", cols,rows2);
        read_temps("QEset_tps"+to_string(t+1)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat", t2, Parts);
        //        cout << t2 << endl;
        index_id = flann_build_index(set2, rows2, cols, &speedup, &p);
        //
        //double moyenne=0;
        //if(t==0) cout << "aaaaaaa"<< endl;
        for(int count=0; count<rows1; count++) {
            //            cout << "aaaa" << endl;
            c.loadData(set1,count,cols); // crée le carnet
            Carnet cbis=c;
            //moyenne+=c.a[0];
            
            bool vide =0;
            vector<int> v;
            if (!all_of(c.a.begin(), --c.a.end(), [](int i) {return i==0;})) {
                v={0,1,2};
                //v={0,1};
            }
            else {
                v={0,1};
                vide=1;
            }
            for( int la : v) {
                for (int lb : v) {
                    c=cbis;
                    //value << "blabla" << endl;
                    matProb[la][lb][t].resize(rows1);
                    //cout << "avant : " << c << endl;
                    //cout << "la " << la <<" lb " << lb << endl;
                    int avant=c.PL();
                    Carnet cc =c;
                    //cout << "avant jump : " << c << endl;
                    //cout << "PL avant " << c.PL() << endl;
                    backP(c,la,lb);
                    //cout << "avant jump backP : " << c << endl;
                    c.Jump(la, lb);
                    //cout << "apres jump : "<< c << endl;
                    int apres=c.PL();
//                    if(avant!=apres) { cout << cc << endl ;
//                        cout << c << endl;
//                        cout << "avant " << avant << " apres : " << apres << endl;
//                    }
                    //cout << "PL apres " << c.PL() << endl;
                    if (all_of(c.a.begin(), --c.a.end(), [](int i) {return i==0;})) {
                        if (c.ra !=-1) {
                            c.pa+=c.K-c.ra;
                        }
                        c.ra=c.K;
                        c.na=c.ainf+1;
                        if (c.rb!=-1) {
                            c.pb-=c.K-c.rb;
                        }
                        c.rb=c.K;
                        c.nb=abs(c.binf)+1;
                        c.a[c.K]=c.ainf;
                        c.b[c.K]=c.binf;
                    }
                    //cout << "apres jump : " << c << endl;
                    c.intensites(la, lb);
                    
                    //cout << "apres jump newP " << c << endl;
                    double Q=c.intMp+c.intMm+accumulate(c.intLa.begin(), c.intLa.end(), 0.)+accumulate(c.intLb.begin(), c.intLb.end(), 0.)+accumulate(c.intCa.begin(), c.intCa.end(), 0.)+accumulate(c.intCb.begin(), c.intCb.end(), 0.);
                    c_=c;
                    //                    cout << c_ << endl;
                    c_.dMp(la, lb);
                    //cout << " a " << c.a << endl;
                    if (all_of(c_.a.begin(), --c_.a.end(), [](int i) {return i==0;})) {
                        if (c_.ra !=-1) {
                            c_.pa+=c_.K-c_.ra;
                        }
                        c_.ra=c_.K;
                        c_.na=c_.ainf+1;
                        if (c_.rb!=-1) {
                            c_.pb-=c_.K-c_.rb;
                        }
                        c_.rb=c_.K;
                        c_.nb=abs(c_.binf)+1;
                        c_.a[c.K]=c_.ainf;
                        c_.b[c.K]=c_.binf;
                    }
                    //cout << "apres dMp : " <<  c_ << endl;
                    c_.toArray(testset);
                    newP(testset, c);
//                    cout << "à projeter : " ;
//                    for (int i =0; i <12; i++) {
//                        cout << testset[i] << " " ;
//                    }
                    //cout << endl;
                    
                    flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
                    ind=*result;
                    //c_.loadData(set2,ind,cols);
                    //cout << "projeté sur : " << c_ << endl;
                    //
                    matProb[la][lb][t][count][ind]+=float(c.intMp)/Q;
                    c_=c;
                    //cout << "apres jump : " << c_ << endl;
                    c_.dMm(la, lb);
                    if (all_of(c_.a.begin(), --c_.a.end(), [](int i) {return i==0;})) {
                        if (c_.ra !=-1) {
                            c_.pa+=c_.K-c_.ra;
                        }
                        c_.ra=c_.K;
                        c_.na=c_.ainf+1;
                        if (c_.rb!=-1) {
                            c_.pb-=c_.K-c_.rb;
                        }
                        c_.rb=c_.K;
                        c_.nb=abs(c_.binf)+1;
                        c_.a[c.K]=c_.ainf;
                        c_.b[c.K]=c_.binf;
                    }
                    //cout << "apres dMm : " <<  c_ << endl;
                    c_.toArray(testset);
                    newP(testset, c);
                    flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
                    ind=*result;
                    //cout << "projeté sur : " << c_ << endl;
                    matProb[la][lb][t][count][ind]+=float(c.intMm)/Q;
                    for (int j =0; j<c.sizeLa; j++) {
                        c_=c;
                        c_.dLp(la,lb, j);
                        if (all_of(c_.a.begin(), --c_.a.end(), [](int i) {return i==0;})) {
                            if (c_.ra !=-1) {
                                c_.pa+=c_.K-c_.ra;
                            }
                            c_.ra=c_.K;
                            c_.na=c_.ainf+1;
                            if (c_.rb!=-1) {
                                c_.pb-=c_.K-c_.rb;
                            }
                            c_.rb=c_.K;
                            c_.nb=abs(c_.binf)+1;
                            c_.a[c.K]=c_.ainf;
                            c_.b[c.K]=c_.binf;
                        }
                        c_.toArray(testset);
                        newP(testset, c);
                        flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
                        ind=*result;
                        matProb[la][lb][t][count][ind]+=float(c.intLa[j])/Q;
                    }
                    for (int j =0; j<c.sizeLb; j++) {
                        c_=c;
                        c_.dLm(la, lb, j);
                        if (all_of(c_.a.begin(), --c_.a.end(), [](int i) {return i==0;})) {
                            if (c_.ra !=-1) {
                                c_.pa+=c_.K-c_.ra;
                            }
                            c_.ra=c_.K;
                            c_.na=c_.ainf+1;
                            if (c_.rb!=-1) {
                                c_.pb-=c_.K-c_.rb;
                            }
                            c_.rb=c_.K;
                            c_.nb=abs(c_.binf)+1;
                            c_.a[c.K]=c_.ainf;
                            c_.b[c.K]=c_.binf;
                        }
                        c_.toArray(testset);
                        newP(testset, c);
                        flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
                        ind=*result;
                        matProb[la][lb][t][count][ind]+=float(c.intLb[j])/Q;
                        //else {matProb[la][lb][t][p1][p2][count][ind]+=float(c.intLb[j])/Q;}
                    }
                    for (int j =0; j<c.sizeCa; j++) {
                        c_=c;
                        c_.dCp(la, lb, j);
                        if (all_of(c_.a.begin(), --c_.a.end(), [](int i) {return i==0;})) {
                            if (c_.ra !=-1) {
                                c_.pa+=c_.K-c_.ra;
                            }
                            c_.ra=c_.K;
                            c_.na=c_.ainf+1;
                            if (c_.rb!=-1) {
                                c_.pb-=c_.K-c_.rb;
                            }
                            c_.rb=c_.K;
                            c_.nb=abs(c_.binf)+1;
                            c_.a[c.K]=c_.ainf;
                            c_.b[c.K]=c_.binf;
                        }
                        c_.toArray(testset);
                        newP(testset, c);
                        flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
                        ind=*result;
                        matProb[la][lb][t][count][ind]+=float(c.intCa[j])/Q;
                        
                    }
                    for (int j =0; j<c.sizeCb; j++) {
                        c_=c;
                        c_.dCm(la, lb, j);
                        if (all_of(c_.a.begin(), --c_.a.end(), [](int i) {return i==0;})) {
                            if (c_.ra !=-1) {
                                c_.pa+=c_.K-c_.ra;
                            }
                            c_.ra=c_.K;
                            c_.na=c_.ainf+1;
                            if (c_.rb!=-1) {
                                c_.pb-=c_.K-c_.rb;
                            }
                            c_.rb=c_.K;
                            c_.nb=abs(c_.binf)+1;
                            c_.a[c.K]=c_.ainf;
                            c_.b[c.K]=c_.binf;
                        }
                        c_.toArray(testset);
                        newP(testset, c);
                        flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
                        ind=*result;
                        matProb[la][lb][t][count][ind]+=float(c.intCb[j])/Q;
                    }
                }
            }
            
        }
        //cout << "building Value ... 1 " << endl;
        //cout << "t : " << t << " : " <<  moyenne/rows1;
        flann_free_index(index_id, &p);
        //
        //            cout << "baba" << endl;
        //            //cout << "Valeurs "<< t << endl;
        Valeurs[t].resize(parts);
        StratOpt[t].resize(parts);
        //cout << "building Value ... 2 " << endl;
        for (int p1=0; p1<parts; p1 ++) {
            ////                if(t==0) cout << "bb" << endl;
            ////                //cout << "t1[p1] " << t1[p1] << endl;
            Valeurs[t][p1].resize(rows1);
            StratOpt[t][p1].resize(rows1);
            //cout << "building Value ... 2a " << endl;
            ////                //cout << "p1 "<< p1 << endl;
            for(int i = 0; i < rows1; i++ ) {
                Carnet c;
                c.loadData(set1,i,cols);
                //cout << c << endl;
                if( t1[p1]<T){
//                    vector<vector<double>> temp;
//                    temp.resize(2);
//                    for (int j=0; j<2; j++) {
//                        temp[j].resize(2);
//                    }
                    vector<vector<double>> temp;
                    temp.resize(3);
                    for (int j=0; j<3; j++) {
                        temp[j].resize(3);
                    }
                    int ind=projeteT(t1[p1], t2);
                    //cout <<"t1[p1], t2 " << t1[p1] << " ; " << t2 << endl;
                    //cout << ind << endl;
                    ////                        //if (ind ==0) cout << ind << "aaaaazzzz " << endl ;
                    ////                        cout << "baba" << endl;
                    bool vide =0;
                    vector<int> v;
                    if (!all_of(c.a.begin(), --c.a.end(), [](int i) {return i==0;})) {
                        v={0,1,2};
                        //v={0,1};
                    }
                    else {
                        v={0,1};
                        vide=1;
                    }
                    for( int la : v) {
                        for (int lb : v) {
                            //value << "blab" << endl;
                            Carnet cc=c;
                            backP(c, la, lb);
                            c.Jump(la, lb);
                            if (all_of(c.a.begin(), --c.a.end(), [](int i) {return i==0;})) {
                                if (c.ra !=-1) {
                                    c.pa+=c.K-c.ra;
                                }
                                c.ra=c.K;
                                c.na=c.ainf+1;
                                if (c.rb!=-1) {
                                    c.pb-=c.K-c.rb;
                                }
                                c.rb=c.K;
                                c.nb=abs(c.binf)+1;
                                c.a[c.K]=c.ainf;
                                c.b[c.K]=c.binf;
                            }
                            //cout << "compute intensity " << la << " " << lb << endl;
                            c.intensites(la, lb);
                            newP(c);
                            //cout << "intensity computed" << endl;
                            double Q=c.intMp+c.intMm+accumulate(c.intLa.begin(), c.intLa.end(), 0.)+accumulate(c.intLb.begin(), c.intLb.end(), 0.)+accumulate(c.intCa.begin(), c.intCa.end(), 0.)+accumulate(c.intCb.begin(), c.intCb.end(), 0.);
                            double s = 0;
                            if(ind<Parts-1) {
                                for (int p2=ind+1; p2<Parts-1 ; p2++) {
                                    for(auto const& it : matProb[la][lb][t][i]) {
                                        s+= (exp(-Q*(t2[p2-1] - t1[p1]))-exp(-Q*(t2[p2] - t1[p1])))*it.second * Valeurs[t+1][p2][it.first];
                                    }
                                }
                                for(auto const& it : matProb[la][lb][t][i]) {
                                    s+= (1-exp(-Q*(t2[ind] - t1[p1])))*it.second * Valeurs[t+1][ind][it.first];
                                }
                                
                                for(auto const& it : matProb[la][lb][t][i]) {
                                    assert(t2[Parts-2] - t1[p1]>0);
                                    s+= (exp(-Q*(t2[Parts-2] - t1[p1]))-exp(-Q*(T-t1[p1])))*it.second * Valeurs[t+1][Parts-1][it.first];
                                }
                            }
                            else {
                                for(auto const& it : matProb[la][lb][t][i]) {
                                    //assert(t2[Parts-2] - t1[p1]>0);
                                    s+= (1-exp(-Q*(T-t1[p1])))*it.second * Valeurs[t+1][Parts-1][it.first];
                                }
                            }
                            //cout << s << endl;
                            //cout << " s " << s << endl;
                            Carnet ccc;
                            ccc=c;
                            ccc.Jump(0, 0);
                            temp[la][lb]=exp(-Q*(T-t1[p1]))*ccc.PL()+s;
                            value << "temp " << la << " " << lb << " : " << temp[la][lb] << " ; s : " << s  ;
                            value << endl;
                            assert(T-t1[p1]>0);
                        }
                    }
                    StratOpt[t][p1][i]=argmax(temp,vide);
                    Valeurs[t][p1][i]= temp[StratOpt[t][p1][i].first][StratOpt[t][p1][i].second];
                    if(Valeurs[t][p1][i]>8) { cout << "i : " << i << " ; " <<Valeurs[t][p1][i];}
                    ////                        //cout << Valeurs[t][p1][i] << " " ;
                }
                else {
                    Valeurs[t][p1][i]=0;
                }
                ////                    //if(t==N-1) cout <<Valeurs[t][p1][i] << " ";
            }
            //cout << "building Value ... 2b " << endl;
        }
        //cout << "building Value ... 3" << endl;
        
        /*
        if(t==0) {
            for(int la :{0,1}) {
                for (int lb : {0,1}) {
                    cout << "la , lb " << la << " " << lb << endl;
                    for(auto& it : matProb[la][lb][0][0] ){
                        Carnet c;
                        int ind=it.first;
                        c.loadData(set2, ind, 12);
                        cout <<  c << " " << it.second << endl;
                    }
                }
            }
            cout << "--------------------------" << endl;
            Carnet c;
            c.loadData(set1, 0, 12);
            Carnet cbis;
            cbis=c;
            for(int la :{0,1}) {
                for (int lb : {0,1}) {
                    c=cbis;
                    cout << "la lb " << la << " " << lb << endl;
                    cout << " c : " << c << endl;
            backP(c,la,lb);
            cout << "avant jump backP : " << c << endl;
            c.Jump(la, lb);
            cout << "apres jump : "<< c << endl;
                    //                    if(avant!=apres) { cout << cc << endl ;
            //                        cout << c << endl;
            //                        cout << "avant " << avant << " apres : " << apres << endl;
            //                    }
            //cout << "PL apres " << c.PL() << endl;
            if (all_of(c.a.begin(), --c.a.end(), [](int i) {return i==0;})) {
                if (c.ra !=-1) {
                    c.pa+=c.K-c.ra;
                }
                c.ra=c.K;
                c.na=c.ainf+1;
                if (c.rb!=-1) {
                    c.pb-=c.K-c.rb;
                }
                c.rb=c.K;
                c.nb=abs(c.binf)+1;
                c.a[c.K]=c.ainf;
                c.b[c.K]=c.binf;
            }
            cout << "apres jump : " << c << endl;
            c.intensites(la, lb);
            
            //cout << "apres jump newP " << c << endl;
            double Q=2*c.intM+accumulate(c.intLa.begin(), c.intLa.end(), 0.)+accumulate(c.intLb.begin(), c.intLb.end(), 0.)+accumulate(c.intCa.begin(), c.intCa.end(), 0.)+accumulate(c.intCb.begin(), c.intCb.end(), 0.);
            c_=c;
            //                    cout << c_ << endl;
            c_.dMp(la, lb);
            //cout << " a " << c.a << endl;
            if (all_of(c_.a.begin(), --c_.a.end(), [](int i) {return i==0;})) {
                if (c_.ra !=-1) {
                    c_.pa+=c_.K-c_.ra;
                }
                c_.ra=c_.K;
                c_.na=c_.ainf+1;
                if (c_.rb!=-1) {
                    c_.pb-=c_.K-c_.rb;
                }
                c_.rb=c_.K;
                c_.nb=abs(c_.binf)+1;
                c_.a[c.K]=c_.ainf;
                c_.b[c.K]=c_.binf;
            }
            cout << "apres dMp : " <<  c_ << endl;
            c_.toArray(testset);
            newP(testset, c);
            index_id = flann_build_index(set2, rows2, cols, &speedup, &p);
            flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
                    ind=*result;
            c_.loadData(set2,ind,cols);
            cout << "projeté sur : " << c_ << endl;
                }
            }
        }
        */
        
        
        delete[] set1;
        delete[] set2;
        delete[] testset;
        delete[] result;
        delete[] dists;
    }
    
    void fillmat() {
        Valeurs[N].resize(Parts);
        Carnet c;
        int cols=8+2*c.K;
        float* set2;
        int rows2;
        vector<double> t2;
        t2.resize(Parts);
        for (int p =0; p<Parts; p++) {
            //cout << " ahah" << endl;
            set2=read_points("set_tps"+to_string(N)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat", cols,rows2);
            read_temps("set_tps"+to_string(N)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat", t2, Parts);
            Valeurs[N][p].resize(rows2);
            double moyenne =0;
            for (int i=0; i<rows2; i++) {
                
                Carnet cc;
                cc.loadData(set2, i, 12);
                //cout << cc << endl;
                //cout << cc.PL() << endl;
                moyenne+=cc.PL();
                //cout << cc.PL() << " " ;
                //c.loadData(set2[p],i,cols);
                //Valeurs[N][p][i]=c.PL();
                Valeurs[N][p][i]=0;
                //cout << Valeurs[N][p][i] << " " ;
            }
            //cout << moyenne/rows2 << endl;
        }
        delete []  set2;
        for (int t=N-1; t>-1; t--) {
            cout << t << " " << flush;
            fillmat(t);
//            if (t==1) {
//                for(auto k : Valeurs[1][0]) {
//                    cout << k << " " ;
//                }
//            }
            //cout << "start cleaning ..." << endl;
            Valeurs[t+1].clear();
            if(0<t and t<N-1) {
                for(int la :{0,1,2}){
                    for( int lb : {0,1,2}) {
                        matProb[la][lb][t+1].clear();
                    }
                }
            }
            //cout << Valeurs[t][0] << endl;
        }

        cout << "Valeur0 : " << Valeurs[0][0][0] << endl;
    }
    
    void writeStrat() {
        ofstream fout("Strat_tps"+to_string(0)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat");
        fout<<StratOpt[0][0][0].first << " " << StratOpt[0][0][0].second<< endl;
        fout.close();
        for (int t=1; t<N; t++) {
            ofstream fout("Strat_tps"+to_string(t)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat");
            for (int i=0; i<Parts; i++) {
                for (int k=0; k<StratOpt[t][i].size(); k++) {
                    fout <<StratOpt[t][i][k].first << " " << StratOpt[t][i][k].second << " " ;
                }
                fout<< endl;
            }
        }
    }
    
};


class Test {
public:
    int NbSimul;
    int N;
    double T;
    int Parts;
    vector<vector<pair<int, int>>> StratOpt;
    int cols;
    
    Test(int nbSimul=10000,int n=500, double t=10, int Cols = 12, int parts=1) : NbSimul(nbSimul), N(n), T(t), cols(Cols) , Parts(parts) {
        Carnet c;
        cols= 8+2*c.K;
    }
    
    
    float* readSet(int t, int& rows) {
        float* p;
        float* data;
        ifstream setin("set_tps"+to_string(t)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat");
        if (!setin) {
            cout << "Cannot open input Set." << endl;
            exit(1);
        }
        setin >> rows;
        
        data = new float[rows*cols];
        if (!data) {
            cout << "Cannot allocate memory." << endl;
            exit(1);
        }
        p = data;
        
        for (int i=0;i<rows;++i) {
            for (int j=0;j<cols;++j) {
                setin >> *p;
                p++;
            }
        }
        return data;
    }
    
    vector<double> simuleStratOpt(int nbSimul=1000) {
        float* set;
        int taille;
        float* test = new float[cols];
        int nn=1;
        int tcount=1;
        int* result = new int[tcount*nn];
        float* dists =  new float[tcount*nn];
        struct FLANNParameters p;
        float speedup;
        flann_index_t index_id;
        p = DEFAULT_FLANN_PARAMETERS;
        p.cores=1;
        p.algorithm = FLANN_INDEX_KDTREE;
        p.trees = 8;
        p.log_level = FLANN_LOG_INFO;
        p.checks = 64;
        
        vector<CarnetControle> carnets(nbSimul);
        Carnet c;
        
        pair<int, int> temp;
        
        for (int t=0; t<N; t++) {
            int parts=Parts;
            if(t==0) parts=1;
            StratOpt.resize(parts);
            
            // Lecture de set_t
            set=readSet(t, taille);
            //cout << "taille=" << taille << endl;
            
            //            for (int l=0; l<taille*12; l++) {
            //                cout << set[l] << endl;
            //            }
            
            ifstream stratin("Strat_tps"+to_string(t)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat");
            if (!stratin) {
                cout << "Cannot open strat file. ahah" << endl;
                exit(1);
            }
            for (int i =0; i<parts; i++) {
                StratOpt[i].resize(taille);
                for (int k=0; k<taille; k++) {
                    stratin >> StratOpt[i][k].first;
                    stratin >> StratOpt[i][k].second;
                }
            }
            stratin.close();
            
            
            index_id = flann_build_index(set, taille, cols, &speedup, &p);
            
            //if (t==1) cout << StratOpt[1][1].first << " " << StratOpt[1][1].second << endl;
            
            vector<double> t1;
            t1.resize(parts);
            read_temps("QEset_tps"+to_string(t)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat", t1, parts);
            //cout << t1 <<endl;
            //cout << t1 << endl;
            //cout << "taille : " << t1.size() << endl;
            //cout << "t : " << t << endl;
            for (int j=0; j<nbSimul; j++) {
                if(carnets[j].Tps<T) {
                    carnets[j].toArray(test);
                    //cout << "carnet " << carnets[j];
                    newP(test,c);
                    //cout << "nouveau carnet " ;
//                    for(int i =0 ; i<12; i++) {
//                        cout << test[i] << " " ;
//                    }
//                    cout << endl;
                    //cout << "avant : " << c << endl;
                    flann_find_nearest_neighbors_index(index_id, test, tcount, result, dists, nn, &p);
                    //cout << "proj :"<< *result << endl;
                    int proj= projeteT(carnets[j].Tps, t1);
                    Carnet cc;
                    cc.loadData(set, *result, 12);
                    //cout << "tour :" << t << endl;
                    //cout << "projeté sur : " << cc<< endl;
                    //cout << " carnets[j].Tps " << carnets[j].Tps << endl;
                    //cout << "Projeté sur " <<t1[proj] << endl;
                    carnets[j].la=StratOpt[proj][*result].first;
                    carnets[j].lb=StratOpt[proj][*result].second;
                    //cout << StratOpt[proj][*result].first << " " << StratOpt[proj][*result].second << endl;
                    //cout <<"strat " << StratOpt[proj][*result].first << " " << StratOpt[proj][*result].second << endl ;
                    carnets[j].Jump(carnets[j].la,carnets[j].lb);
                    if (all_of(carnets[j].a.begin(), --carnets[j].a.end(), [](int i) {return i==0;})) {
                        //cout << "modifie" << endl;
                        if (carnets[j].ra !=-1) {
                            carnets[j].pa+=carnets[j].K-carnets[j].ra;
                        }
                        carnets[j].ra=carnets[j].K;
                        carnets[j].na=carnets[j].ainf+1;
                        if (carnets[j].rb!=-1) {
                            carnets[j].pb-=carnets[j].K-carnets[j].rb;
                        }
                        carnets[j].rb=carnets[j].K;
                        carnets[j].nb=abs(carnets[j].binf)+1;
                        carnets[j].a[carnets[j].K]=carnets[j].ainf;
                        carnets[j].b[carnets[j].K]=carnets[j].binf;
                    }
                    carnets[j].intensites(carnets[j].la, carnets[j].lb);
                    
                    tuple<int,double> pair=carnets[j].genEv();
                    int event = get<0>(pair);
                    double tau= get<1>(pair);
                    //cout << tau ;
                    if (carnets[j].Tps+tau < T) {
                        //cout <<tau;
                        carnets[j].Tps+=tau;
                        carnets[j].move(event);
                        int a0=carnets[j].A0(0);
                        if (all_of(carnets[j].a.begin(), --carnets[j].a.end(), [](int i) {return i==0;})) {
                            //cout << "modifie" << endl;
                            if (carnets[j].ra !=-1) {
                                carnets[j].pa+=carnets[j].K-carnets[j].ra;
                            }
                            carnets[j].ra=carnets[j].K;
                            carnets[j].na=carnets[j].ainf +1;
                            if (carnets[j].rb!=-1) {
                                carnets[j].pb-=carnets[j].K-carnets[j].rb;
                            }
                            carnets[j].rb=carnets[j].K;
                            carnets[j].nb=abs(carnets[j].binf)+1;
                            carnets[j].a[carnets[j].K]=carnets[j].ainf;
                            carnets[j].b[carnets[j].K]=carnets[j].binf;
                        }
                        int aa0=carnets[j].A1(0);
                        int bb0=carnets[j].B1(0);
                        a0=carnets[j].A0(0);
                        int b0=carnets[j].B0(0);
                        //cout << "apres move " << carnet << endl;
                        //assert(carnet.pa-carnet.pb==carnet.ra+carnet.rb-a0+1);
                        assert(carnets[j].b[carnets[j].K]==-5 && carnets[j].a[carnets[j].K]==5);
                        assert(all_of(carnets[j].a.begin(), --carnets[j].a.end(), [](int i) {return i>=0;}));
                        assert(all_of(carnets[j].b.begin(), --carnets[j].b.end(), [](int i) {return i<=0;}));
                        assert(carnets[j].pa>carnets[j].pb);
                        assert(carnets[j].pa-carnets[j].pb<=carnets[j].A0(0)+1);
                        if(carnets[j].la!=0 && carnets[j].lb !=0) assert(carnets[j].pa-carnets[j].pb==min(carnets[j].A0(0),carnets[j].rb)+min(carnets[j].A0(0),carnets[j].ra)-carnets[j].A0(0)+1);
                        assert(carnets[j].pa-carnets[j].pb<=2*carnets[j].K+1);
                        assert(carnets[j].A0(0)==carnets[j].B0(0));
                        if (carnets[j].la==1) assert(carnets[j].na <= carnets[j].a[a0]+1);
                        if (carnets[j].la==2) assert(carnets[j].na <= max(carnets[j].a[aa0], carnets[j].a[a0])+1);
                        if (carnets[j].lb==1) assert(carnets[j].nb<=abs(carnets[j].b[a0])+1);
                        if (carnets[j].lb==2)assert(carnets[j].nb<=max(abs(carnets[j].b[bb0]),abs(carnets[j].b[b0]))+1);
                    }
                }
                else {
                    carnets[j].Tps =T;
                }
            }
            StratOpt.clear();
            t1.clear();
            flann_free_index(index_id, &p);
            delete [] set;
        }
        
        vector<double> res;
        res.reserve(nbSimul);
        for(auto& el : carnets) {
            res.push_back(el.PL());
        }
        delete [] test;
        delete [] result;
        delete [] dists;
        return res;
    }
    
 
 
 
 
 
 
 
 
 
 

    map<double,map<pair<bool,bool>, int>> plotstratImb(int t) {
        float* set;
        int taille;
        set=readSet(t, taille);

        ifstream stratin("Strat_tps"+to_string(t)+"nbSim"+to_string(NbSimul)+"subd"+to_string(N)+".dat");
        if (!stratin) {
            cout << "Cannot open strat file. ahah" << endl;
            exit(1);
        }
        StratOpt.resize(Parts);
        for (int i =0; i<Parts; i++) {
            StratOpt[i].resize(taille);
            for (int k=0; k<taille; k++) {
                stratin >> StratOpt[i][k].first;
                stratin >> StratOpt[i][k].second;
            }
        }
        stratin.close();
        //cout << "ahah" << endl;
        map<double,map<pair<bool,bool>, int>> stat;
        
        for (int k=0; k<taille; k++) {
            Carnet c;
            c.loadData(set, k, 12);
            //cout << c << endl;
            int a0= c.A0();
            //cout << a0 << endl;
            double bal=(double)c.a[a0]/(abs(c.b[a0])+c.a[a0]);
            //cout << bal << endl;
            for (int i =0; i<Parts; i++) {
                stat[bal][StratOpt[i][k]]+=1;
            }
        }
//        for( auto& i : stat) {
//            cout << "balance : " << i.first << endl;
//            for(auto& j : i.second) {
//                cout << "strat " << j.first.first << " " << j.first.second  << " : " << j.second << endl;
//            }
//        }
        delete[] set;
        return stat;
    }
    
    void plotstratImb() {
        map<double,map<pair<bool,bool>, int>> stat;
        for (int t=1; t<N; t++) {
            map<double,map<pair<bool,bool>, int>> temp=plotstratImb(t);
            for(auto &i : temp) {
                for(auto & j : i.second) {
                    stat[i.first][j.first]+=j.second;
                }
            }
        }
        for( auto& i : stat) {
            cout << "balance : " << i.first << endl;
            for(auto& j : i.second) {
                cout << "strat " << j.first.first << " " << j.first.second  << " : " << j.second << endl;
            }
        }
        pair<bool,bool> temp;
        for (bool i : {0,1}) {
            for(bool j : {0,1}) {
                ofstream outstat("statStrat"+to_string(i)+to_string(j)+".dat");
                temp.first=i;
                temp.second=j;
                for(auto & t : stat) {
                    outstat << t.first << " " << t.second[temp] << endl;
                }
            }
        }
    }
    
        //
        ////    void InventaireOpt(int nbSimul=10) {
        ////        float* set;
        ////        int taille;
        ////        float* test = new float[cols];
        ////        int nn=1;
        ////        int tcount=1;
        ////        int* result = new int[tcount*nn];
        ////        float* dists =  new float[tcount*nn];
        ////        struct FLANNParameters p;
        ////        float speedup;
        ////        flann_index_t index_id;
        ////        p = DEFAULT_FLANN_PARAMETERS;
        ////        p.cores=1;
        ////        p.algorithm = FLANN_INDEX_KDTREE;
        ////        p.trees = 8;
        ////        p.log_level = FLANN_LOG_INFO;
        ////        p.checks = 64;
        ////
        ////        vector<CarnetControle> carnets(nbSimul);
        ////        vector<vector<int>> inventaires(nbSimul);
        ////        ifstream stratin("StratN"+to_string(N)+"NbSimuls"+to_string(NbSimul)+"T"+to_string(T)+".dat");
        ////        if (!stratin) {
        ////            cout << "Cannot open strat file." << endl;
        ////            exit(1);
        ////        }
        ////        pair<bool, bool> temp;
        ////
        ////        for (int t=0; t<N; t++) {
        ////            // Lecture de set_t
        ////            set=readSet(t, taille);
        ////            //cout << "taille=" << taille << endl;
        ////            index_id = flann_build_index(set, taille, cols, &speedup, &p);
        ////
        ////            // Lecture de strat_t
        ////            StratOpt.reserve(taille);
        ////            for(int k=0 ; k<taille ; k++) {
        ////                stratin >> temp.first;
        ////                stratin >> temp.second;
        ////                //cout << temp.first << temp.second;
        ////                StratOpt.push_back(temp);
        ////            }
        ////            //cout << endl;
        ////
        ////
        ////            for (int j=0; j<nbSimul; j++) {
        ////                carnets[j].toArray(test);
        ////                flann_find_nearest_neighbors_index(index_id, test, tcount, result, dists, nn, &p);
        ////                //cout << "proj :"<< *result << endl;
        ////                carnets[j].la=StratOpt[*result].first;
        ////                carnets[j].lb=StratOpt[*result].second;
        ////                tuple<int,double> pair=carnets[j].genEv();
        ////                int event = get<0>(pair);
        ////                double tau= get<1>(pair);
        ////                carnets[j].Tps+=tau;
        ////                while (carnets[j].Tps < t*float(T)/N) {
        ////                    carnets[j].move(event);
        ////                    int a0=carnets[j].A0();
        ////                    if (all_of(carnets[j].a.begin(), --carnets[j].a.end(), [](int i) {return i==0;})) {
        ////                        if (carnets[j].ra==carnets[j].K) {
        ////                            carnets[j].na=carnets[j].ainf+1;
        ////                        }
        ////                        if (carnets[j].rb==carnets[j].K) {
        ////                            carnets[j].nb=abs(carnets[j].binf)+1;
        ////                        }
        ////                        carnets[j].a[carnets[j].K]=carnets[j].ainf;
        ////                        carnets[j].b[carnets[j].K]=carnets[j].binf;
        ////                    }
        ////                    assert(carnets[j].pa-carnets[j].pb==carnets[j].ra+carnets[j].rb-a0+1);
        ////                    assert(carnets[j].na <= carnets[j].a[a0]+1);
        ////                    assert(carnets[j].nb<=abs(carnets[j].b[a0])+1);
        ////                    tuple<int,double> pair=carnets[j].genEv();
        ////                    event = get<0>(pair);
        ////                    tau= get<1>(pair);
        ////                    carnets[j].Tps+=tau;
        ////                }
        ////                carnets[j].Tps=t*float(T)/N;
        ////                inventaires[j].push_back(carnets[j].y);
        ////            }
        ////            StratOpt.clear();
        ////            delete [] set;
        ////            flann_free_index(index_id, &p);
        ////        }
        ////
        ////        for (int j=0; j<nbSimul; j++) {
        ////            tuple<int,double> pair=carnets[j].genEv();
        ////            int event = get<0>(pair);
        ////            double tau= get<1>(pair);
        ////            carnets[j].Tps+=tau;
        ////            carnets[j].move(event);
        ////            if (all_of(carnets[j].a.begin(), --carnets[j].a.end(), [](int i) {return i==0;})) {
        ////                if (carnets[j].ra==carnets[j].K) {
        ////                    carnets[j].na=carnets[j].ainf+1;
        ////                }
        ////                if (carnets[j].rb==carnets[j].K) {
        ////                    carnets[j].nb=abs(carnets[j].binf)+1;
        ////                }
        ////                carnets[j].a[carnets[j].K]=carnets[j].ainf;
        ////                carnets[j].b[carnets[j].K]=carnets[j].binf;
        ////            }
        ////            inventaires[j].push_back(carnets[j].y);
        ////        }
        ////        vector<double> res;
        ////        res.reserve(nbSimul);
        ////        for(auto& el : carnets) {
        ////            res.push_back(el.PL());
        ////        }
        ////
        ////        delete [] test;
        ////        delete [] result;
        ////        delete [] dists;
        ////
        ////        ofstream fout("inventaire.dat");
        ////        for (int j=0; j<nbSimul; j++) {
        ////            for (auto& it : inventaires[j]) {
        ////                fout << it << " ";
        ////            }
        ////            fout << endl;
        ////        }
        ////    }
    };
//
//    
//    
//    
//    */
//    
//    //class Test00 { // Calcul d'esperances de fonction du carnet d'ordre à maturité
//    //public:
//    //    int N;
//    //    double T;
//    //    int cols;
//    //    double h;
//    //    vector<vector<map<int,double>>> matProb;
//    //    vector<float*> Valeurs;
//    //    double Q;
//    //
//    //    Test00(int n=500, double t=10, int Cols = 12) : N(n), T(t), cols(Cols) {
//    //        Valeurs.resize(N+1);
//    //        matProb.resize(N);
//    //        h=T/N;
//    //        cout << "h=" << h << endl;
//    //    }
//    //
//    //    void fillmatT() {
//    //        //double bench=0;
//    //        Carnet c;
//    //        Carnet c_;
//    //        float* testset = new float[12];
//    //        int tcount = 1;
//    //        int nn=1;
//    //        int* result = new int[tcount*nn];
//    //        float* dists =  new float[tcount*nn];
//    //        struct FLANNParameters p;
//    //        float speedup;
//    //        flann_index_t index_id;
//    //        p = DEFAULT_FLANN_PARAMETERS;
//    //        p.algorithm = FLANN_INDEX_KDTREE;
//    //        p.trees = 8;
//    //        p.log_level = FLANN_LOG_INFO;
//    //        p.checks = 64;
//    //        float* set2;
//    //        float* set1;
//    //        int rows2;
//    //        int rows1;
//    //        int cols=12;
//    //        set2=read_points("00set_tps"+to_string(N)+".dat", cols,rows2);
//    //
//    //        Valeurs[N] = new float[rows2];
//    //        for (int i=0; i<rows2; i++) {
//    //            c.loadData(set2,i,cols);
//    //            Valeurs[N][i]=c.a[0]; // Esperance de a[0]
//    //        }
//    //        //cout << "esperance : " << bench/rows2 << endl;
//    //        for (int t = N-1; t>-1; --t) {
//    //            cout << t <<" " << flush;
//    //            //cout << "tour : " << t << endl;
//    //            index_id = flann_build_index(set2, rows2, cols, &speedup, &p);
//    //
//    //            set1=read_points("00set_tps"+to_string(t)+".dat", cols,rows1);
//    //
//    //            matProb[t].resize(rows1);
//    //
//    //            for(int count=0; count<rows1; count++) {
//    //                c.loadData(set1,count,cols); // crée le carnet
//    //                bool la =0;
//    //                bool lb=0;
//    //                c.intensites(la, lb);
//    //                Q=2.*c.intM+accumulate(c.intLa.begin(), c.intLa.end(), 0.)+accumulate(c.intLb.begin(), c.intLb.end(), 0.)+accumulate(c.intCa.begin(), c.intCa.end(), 0.)+accumulate(c.intCb.begin(), c.intCb.end(), 0.);
//    //                //cout << "Q=" << Q << endl;
//    //                c_=c;
//    //                //cout << c_ << endl;
//    //                //cout << c_ << " subit dMp donne ";
//    //                c_.dMp(la, lb);
//    //                //cout << c_ << endl;
//    //                //cout << c_ << " et se projete sur ";
//    //                c_.toArray(testset);
//    ////                for (int kk=0; kk<12; kk++) {
//    ////                    cout << testset[kk] << " ";
//    ////                } cout << endl;
//    //                flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                int ind=*result;
//    //                //c_.loadData(set2,ind,cols);
//    //                //cout << c_ << " avec proba " << min(h*Q,1.)*double(c.intM)/Q<<  endl;
//    //                //cout << "projeté sur : " << c_ << endl;
//    //                matProb[t][count][ind]+=min(h*Q,1.)*double(c.intM)/Q;
//    //                c_=c;
//    //                //cout << c_ << endl;
//    //                //cout << c_ << " subit dMm donne ";
//    //                c_.dMm(la, lb);
//    //                //cout << c_ << " et se projete sur ";
//    //                //cout << c_ << endl;
//    //                c_.toArray(testset);
//    //                flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                ind=*result;
//    //                //c_.loadData(set2,ind,cols);
//    //                //cout << c_ << " avec proba " << min(h*Q,1.)*double(c.intM)/Q << endl;
//    //                matProb[t][count][ind]+=min(h*Q,1.)*double(c.intM)/Q;
//    //                for (int j =0; j<c.sizeLa; j++) {
//    //                    c_=c;
//    //                    //cout << c_ << endl;
//    //                    //cout << c_ << " subit dLa " << j << " donne ";
//    //                    c_.dLp(la, lb, j);
//    //                    //cout << c_ << " et se projete sur ";
//    //                    c_.toArray(testset);
//    //                    flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                    ind=*result;
//    //                    //c_.loadData(set2,ind,cols);
//    //                    //cout << c_ << " avec proba " << min(h*Q,1.)*float(c.intLa[j])/Q <<  endl;
//    //                    matProb[t][count][ind]+=min(h*Q,1.)*double(c.intLa[j])/Q;
//    //                }
//    //                for (int j =0; j<c.sizeLb; j++) {
//    //                    c_=c;
//    //                    //cout << c_ << endl;
//    //                    //cout << c_ << " subit dLb " << j << "donne " ;
//    //                    c_.dLm(la, lb, j);
//    //                    //cout << c_ << " et se projete sur ";
//    //                    c_.toArray(testset);
//    //                    flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                    ind=*result;
//    //                    //c_.loadData(set2,ind,cols);
//    //                    //cout << c_ << " avec proba : " << min(h*Q,1.)*float(c.intLb[j])/Q << endl;
//    //                    matProb[t][count][ind]+=min(h*Q,1.)*double(c.intLb[j])/Q;
//    //                }
//    //                for (int j =0; j<c.sizeCa; j++) {
//    //                    c_=c;
//    //                    //cout << c_ << endl;
//    //                    //cout << c_ << " subit dCa" << j <<" donne ";
//    //                    c_.dCp(la, lb, j);
//    //                    //cout << c_ << " et se projete sur ";
//    //                    c_.toArray(testset);
//    //                    flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                    ind=*result;
//    //                    //c_.loadData(set2,ind,cols);
//    //                    //cout << c_ <<  " avec proba " << min(h*Q,1.)*float(c.intCa[j])/Q << endl;
//    //                    matProb[t][count][ind]+=min(h*Q,1.)*double(c.intCa[j])/Q;
//    //                }
//    //                for (int j =0; j<c.sizeCb; j++) {
//    //                    c_=c;
//    //                    //cout << c_ << endl;
//    //                    //cout << c_ << " subit dCb" << j << " donne ";
//    //                    c_.dCm(la, lb, j);
//    //                    //cout << c_ << " et se projete sur ";
//    //                    //cout << "point diffusé :" << c_ << endl;
//    //                    c_.toArray(testset);
//    //                    flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                    ind=*result;
//    //                    //c_.loadData(set2,ind,cols);
//    //                    //cout << c_ << " avec proba " <<min(h*Q,1.)*float(c.intCb[j])/Q << endl;
//    //                    //c_.loadData(set2,ind,cols);
//    //                    //cout << "projeté sur : " << c_ << endl;
//    //                    matProb[t][count][ind]+=min(h*Q,1.)*double(c.intCb[j])/Q;
//    //                }
//    //                c_=c;
//    //                c_.toArray(testset);
//    //                //cout << c << " ne bouge pas avec proba : "<<1-min(h*Q,1.) << endl;
//    //                flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                ind=*result;
//    //                matProb[t][count][ind]+=1.-min(h*Q,1.);
//    //            } // calcule MatProb[t]
//    //            flann_free_index(index_id, &p);
//    //
//    //            Valeurs[t]= new float[rows1];
//    //            for(int i = 0; i < rows1; i++ ) {
//    //                double s = 0;
//    //                for(auto const& it : matProb[t][i]) {
//    //                    //cout << " tps :" << t << "point : " << i << " vers " << it.first << " avec proba = " << it.second << endl;
//    //                    s+=it.second*Valeurs[t+1][it.first];
//    //                }
//    //                //cout << "somme : "<< summ << endl;
//    //                Valeurs[t][i]=s;
//    //            }
//    //            delete [] Valeurs[t+1];
//    //            delete [] set2;
//    //            set2=set1;
//    //            rows2=rows1;
//    //            matProb[t].clear();
//    //        }
//    //        delete [] set1;
//    //        delete[] testset;
//    //        delete[] result;
//    //        delete[] dists;
//    //    }
//    //
//    //    void testFLANN() {
//    //        //double bench=0;
//    //        Carnet c;
//    //        Carnet c_;
//    //        float* testset = new float[12];
//    //        int tcount = 1;
//    //        int nn=1;
//    //        int* result = new int[tcount*nn];
//    //        float* dists =  new float[tcount*nn];
//    //        struct FLANNParameters p;
//    //        float speedup;
//    //        flann_index_t index_id;
//    //        p = DEFAULT_FLANN_PARAMETERS;
//    //        p.algorithm = FLANN_INDEX_KDTREE;
//    //        p.trees = 8;
//    //        p.log_level = FLANN_LOG_INFO;
//    //        p.checks = 64;
//    //        float* set2;
//    //        float* set1;
//    //        int rows2;
//    //        int rows1;
//    //        int cols=12;
//    //        set2=read_points("00set_tps"+to_string(N)+".dat", cols,rows2);
//    //
//    //        Valeurs[N] = new float[rows2];
//    //        for (int i=0; i<rows2; i++) {
//    //            c.loadData(set2,i,cols);
//    //            Valeurs[N][i]=c.a[0]; // Esperance de a[0]
//    //        }
//    //        //cout << "esperance : " << bench/rows2 << endl;
//    //        for (int t = N-1; t>-1; --t) {
//    //            //cout << "tour : " << t << endl;
//    //            index_id = flann_build_index(set2, rows2, cols, &speedup, &p);
//    //
//    //            set1=read_points("00set_tps"+to_string(t)+".dat", cols,rows1);
//    //
//    //            matProb[t].resize(rows1);
//    //
//    //            for(int count=0; count<rows1; count++) {
//    //                c.loadData(set1,count,cols); // crée le carnet
//    //                bool la =0;
//    //                bool lb=0;
//    //                c.intensites(la, lb);
//    //                Q=2*c.intM+accumulate(c.intLa.begin(), c.intLa.end(), 0.)+accumulate(c.intLb.begin(), c.intLb.end(), 0.)+accumulate(c.intCa.begin(), c.intCa.end(), 0.)+accumulate(c.intCb.begin(), c.intCb.end(), 0.);
//    //                //cout << Q << endl;
//    //                c_=c;
//    //                //cout << c_ << endl;
//    //                //cout << c_ << " subit dMp donne ";
//    //                c_.dMp(la, lb);
//    //                //cout << c_ << " et se projete sur ";
//    //                c_.toArray(testset);
//    //                flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                int ind=*result;
//    //                c_.loadData(set2,ind,cols);
//    //                //cout << c_  << " avec proba " << min(h*Q,1.)*double(c.intM)/Q<<  endl;
//    //                //cout << "projeté sur : " << c_ << endl;
//    //                matProb[t][count][ind]+=min(h*Q,1.)*double(c.intM)/Q;
//    //                c_=c;
//    //                //cout << c_ << endl;
//    //                //cout << c_ << " subit dMm donne ";
//    //                c_.dMm(la, lb);
//    //                //cout << c_ << " et se projete sur ";
//    //                //cout << c_ << endl;
//    //                c_.toArray(testset);
//    //                flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                ind=*result;
//    //                //c_.loadData(set2,ind,cols);
//    //                //cout << c_ << " avec proba " << min(h*Q,1.)*double(c.intM)/Q << endl;
//    //                matProb[t][count][ind]+=min(h*Q,1.)*double(c.intM)/Q;
//    //                for (int j =0; j<c.sizeLa; j++) {
//    //                    c_=c;
//    //                    //cout << c_ << endl;
//    //                    //cout << c_ << " subit dLa " << j << " donne ";
//    //                    c_.dLp(la, lb, j);
//    //                    //cout << c_ << " et se projete sur ";
//    //                    c_.toArray(testset);
//    //                    flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                    ind=*result;
//    //                    //c_.loadData(set2,ind,cols);
//    //                    //cout << c_ << " avec proba " << min(h*Q,1.)*float(c.intLa[j])/Q <<  endl;
//    //                    matProb[t][count][ind]+=min(h*Q,1.)*double(c.intLa[j])/Q;
//    //                }
//    //                for (int j =0; j<c.sizeLb; j++) {
//    //                    c_=c;
//    //                    //cout << c_ << endl;
//    //                    //cout << c_ << " subit dLb " << j << "donne " ;
//    //                    c_.dLm(la, lb, j);
//    //                    //cout << c_ << " et se projete sur ";
//    //                    c_.toArray(testset);
//    //                    flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                    ind=*result;
//    //                    //c_.loadData(set2,ind,cols);
//    //                    //cout << c_ << " avec proba : " << min(h*Q,1.)*float(c.intLb[j])/Q << endl;
//    //                    matProb[t][count][ind]+=min(h*Q,1.)*double(c.intLb[j])/Q;
//    //                }
//    //                for (int j =0; j<c.sizeCa; j++) {
//    //                    c_=c;
//    //                    //cout << c_ << endl;
//    //                    //cout << c_ << " subit dCa" << j <<" donne ";
//    //                    c_.dCp(la, lb, j);
//    //                    //cout << c_ << " et se projete sur ";
//    //                    c_.toArray(testset);
//    //                    flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                    ind=*result;
//    //                    //c_.loadData(set2,ind,cols);
//    //                    //cout << c_ <<  " avec proba " << min(h*Q,1.)*float(c.intCa[j])/Q << endl;
//    //                    matProb[t][count][ind]+=min(h*Q,1.)*double(c.intCa[j])/Q;
//    //                }
//    //                for (int j =0; j<c.sizeCb; j++) {
//    //                    c_=c;
//    //                    //cout << c_ << endl;
//    //                    //cout << c_ << " subit dCb" << j << " donne ";
//    //                    c_.dCm(la, lb, j);
//    //                    //cout << c_ << " et se projete sur ";
//    //                    //cout << "point diffusé :" << c_ << endl;
//    //                    c_.toArray(testset);
//    //                    flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                    ind=*result;
//    //                    //c_.loadData(set2,ind,cols);
//    //                    //cout << c_ << " avec proba " <<min(h*Q,1.)*float(c.intCb[j])/Q << endl;
//    //                    //c_.loadData(set2,ind,cols);
//    //                    //cout << "projeté sur : " << c_ << endl;
//    //                    matProb[t][count][ind]+=min(h*Q,1.)*double(c.intCb[j])/Q;
//    //                }
//    //                c_=c;
//    //                c_.toArray(testset);
//    //                //cout << c << " ne bouge pas avec proba : "<<1-min(h*Q,1.) << endl;
//    //                flann_find_nearest_neighbors_index(index_id, testset, tcount, result, dists, nn, &p);
//    //                ind=*result;
//    //                matProb[t][count][ind]+=1.-min(h*Q,1.);
//    //            } // calcule MatProb[t]
//    //            flann_free_index(index_id, &p);
//    //
//    //            Valeurs[t]= new float[rows1];
//    //            for(int i = 0; i < rows1; i++ ) {
//    //                double s = 0;
//    //                for(auto const& it : matProb[t][i]) {
//    //                    //cout << " tps :" << t << "point : " << i << " vers " << it.first << " avec proba = " << it.second << endl;
//    //                    s+=it.second*Valeurs[t+1][it.first];
//    //                }
//    //                //cout << "somme : "<< summ << endl;
//    //                Valeurs[t][i]=s;
//    //            }
//    //            delete [] Valeurs[t+1];
//    //            delete [] set2;
//    //            set2=set1;
//    //            rows2=rows1;
//    //            matProb[t].clear();
//    //        }
//    //        delete [] set1;
//    //        delete[] testset;
//    //        delete[] result;
//    //        delete[] dists;
//    //    }
//    //};
//    
//
//    
    class Test11 {  //Calcul d'esperances de fonction du carnet d'ordre à maturité
        
    public:
        int NbSimul;
        int N;
        double T;
        int Parts;
        int cols;
        vector<Carnet> points;
        
        Test11(int nbSimul=10000,int n=500, double t=10, int Cols = 12, int parts=1) : NbSimul(nbSimul), N(n), T(t), cols(Cols) , Parts(parts) {
            points.reserve(NbSimul);
        }
        
        vector<double> simule11() {
            double record;
            CarnetControle carnet;
            CarnetControle c_init;
            
            for (int j =0 ; j<NbSimul; j++) {
                carnet=c_init;
                for (int i =1; i<N+1; i++) {
                    //cout << "tour : " << i << endl;
                    //if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                    carnet.la=1;
                    carnet.lb=1;
                    //                }
                    //                else {
                    //                    carnet.la=Unif2(carnet.eng);
                    //                    carnet.lb=Unif2(carnet.eng);
                    //                }
                    
                    
                    
                    //cout << "avant: " << carnet << endl;
                    //carnet.intensites(carnet.la, carnet.lb);
                    carnet.Jump(carnet.la,carnet.lb);
                    if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                        //cout << "modifie" << endl;
                        if (carnet.ra !=-1) {
                            carnet.pa+=carnet.K-carnet.ra;
                        }
                        carnet.ra=carnet.K;
                        carnet.na=carnet.ainf+1;
                        if (carnet.rb!=-1) {
                            carnet.pb-=carnet.K-carnet.rb;
                        }
                        carnet.rb=carnet.K;
                        carnet.nb=abs(carnet.binf)+1;
                        carnet.a[carnet.K]=carnet.ainf;
                        carnet.b[carnet.K]=carnet.binf;
                    }
                    carnet.intensites(carnet.la, carnet.lb);
                    //cout << "tour :" << i << " temps: "<< i*float(T)/N << endl;
                    tuple<int,double> pair=carnet.genEv();
                    int event = get<0>(pair);
                    double tau= get<1>(pair);
                    if (carnet.Tps+tau <T) {
                        carnet.Tps+=tau;
                        //cout << "apres jump : " << carnet << endl;
                        carnet.move(event, record);
                        int a0=carnet.A0(0);
                        if (all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i==0;})) {
                            //cout << "modifie" << endl;
                            if (carnet.ra !=-1) {
                                carnet.pa+=carnet.K-carnet.ra;
                            }
                            carnet.ra=carnet.K;
                            carnet.na=carnet.ainf +1;
                            if (carnet.rb!=-1) {
                                carnet.pb-=carnet.K-carnet.rb;
                            }
                            carnet.rb=carnet.K;
                            carnet.nb=abs(carnet.binf)+1;
                            carnet.a[carnet.K]=carnet.ainf;
                            carnet.b[carnet.K]=carnet.binf;
                        }
                        int aa0=carnet.A1(0);
                        int bb0=carnet.B1(0);
                        a0=carnet.A0(0);
                        int b0=carnet.B0(0);
                        //cout << "apres move " << carnet << endl;
                        //assert(carnet.pa-carnet.pb==carnet.ra+carnet.rb-a0+1);
                        assert(carnet.b[carnet.K]==-5 && carnet.a[carnet.K]==5);
                        assert(all_of(carnet.a.begin(), --carnet.a.end(), [](int i) {return i>=0;}));
                        assert(all_of(carnet.b.begin(), --carnet.b.end(), [](int i) {return i<=0;}));
                        assert(carnet.pa>carnet.pb);
                        assert(carnet.pa-carnet.pb<=carnet.A0(0)+1);
                        if(carnet.la!=0 && carnet.lb !=0) assert(carnet.pa-carnet.pb==min(carnet.A0(0),carnet.rb)+min(carnet.A0(0),carnet.ra)-carnet.A0(0)+1);
                        assert(carnet.pa-carnet.pb<=2*carnet.K+1);
                        assert(carnet.A0(0)==carnet.B0(0));
                        if (carnet.la==1) assert(carnet.na <= carnet.a[a0]+1);
                        if (carnet.la==2) assert(carnet.na <= max(carnet.a[aa0], carnet.a[a0])+1);
                        if (carnet.lb==1) assert(carnet.nb<=abs(carnet.b[a0])+1);
                        if (carnet.lb==2)assert(carnet.nb<=max(abs(carnet.b[bb0]),abs(carnet.b[b0]))+1);
                    }
                    else {
                        carnet.Tps=T;
                        
                    }
                }
                points.push_back(carnet);
            }
            cout << "dMp moyen " << record/NbSimul;
            vector<double> res;
            res.reserve(NbSimul);
            for(auto& el : points) {
                double pl=el.PL();
                res.push_back(el.PL());
            }
            return res;
        }
        
    };
//
//
//
//    
//    /*
//    
//    void convergenceQuantif(int subd, double T,int parts, int NbSimultest) {
//        ofstream fout;
//        fout.open("convergenceQuantif.dat", std::ios_base::app);
//        #pragma omp parallel for num_threads(1) schedule (dynamic)
//        for (int NbSimul=1000; NbSimul<20001; NbSimul+=1000) {
//            genereData(NbSimul,subd,T,parts);
//            Quantif3 quant(NbSimul,subd,T,parts);
//            quant.fillmat();
//            quant.writeStrat();
//            quant.StratOpt.clear();
//            Test test(NbSimul,subd,T,12,parts);
//            vector<double> res =test.simuleStratOpt(NbSimultest);
//            ofstream outres("PLoptT"+to_string(T)+"NbSimul"+to_string(NbSimul)+"parts"+to_string(parts)+"subd"+to_string(subd)+".dat");
//            for(auto i :res) {
//                outres << i << endl;
//            }
//        }
//    }
//    //}
//    //void valuefct(double T) {
//    //    ofstream fout;
//    //    int NbSimul=50000;
//    //    fout.open("valuefct.dat", std::ios_base::app);
//    //#pragma omp parallel for num_threads(10) schedule (dynamic)
//    //    for (double t=.01; t<T; t+=.01) {
//    //        int subd=t*1000;
//    //        genereData(NbSimul,subd,t);
//    //        Quantif3 quant(NbSimul,subd,t);
//    //        quant.fillmatT();
//    //        fout << t << " " << quant.Valeurs[0][0] << endl;
//    //    }
//    //}
//    //
//    //void writeResults(vector<int>& pl, int N, int nbtests, double T) {
//    //    ofstream fout("PL_N"+to_string(N)+"subd.dat");
//    //    if (!fout) {
//    //        cout << "Cannot open output file." << endl;
//    //        exit(1);
//    //    }
//    //        for (int i=0; i<nbtests; i++) {
//    //            fout<< pl[i] << "\n";
//    //        }
//    //}
//    //void writeResults(vector<double>& pl, int N, int nbtests, double T, string k) {
//    //    ofstream fout(k+to_string(N)+"subd.dat");
//    //    if (!fout) {
//    //        cout << "Cannot open output file." << endl;
//    //        exit(1);
//    //    }
//    //    for (int i=0; i<nbtests; i++) {
//    //        fout<< pl[i] << "\n";
//    //    }
//*/
//    


    std::mt19937 Carnet::eng{std::random_device{}()};
    
    int main(int argc, const char * argv[]) {
        int NbSimul=2000000;
        int parts = 20;
        int subd = 160;
        double T=8.;
        int NbSimulTest = 100000;
        
        CarnetControle c;
        
        
//        cout << !(1==0)*4 << endl;
//
//        Carnet c;
//        c.b[0]=0;
//        cout << c << endl;
//        cout << c.B1(0) << endl;
//        c.b[0]=0;
//        c.a[1]=1;
//        c.a[2]=1;
//        c.a[3]=4;
//        //c.a[0]=3;
//        cout << c <<endl;
//        cout << c.b[c.B1(0)] << endl;
//        cout << c.A0(0) << endl;
//        cout << c.B1(0) << endl;
//        c.Jump(1,2);
//        cout << c;
        
//        double moyenne = 0;
//        auto t1=chrono::high_resolution_clock::now();
//        for (int k=0; k<10000; k++) {
//            Proc proc(200,10.);
//            proc.simule();
//            cout <<proc.proc[25] << endl;
//            moyenne+=proc.proc[25].a[0];
//        }
//        cout << moyenne/10000;
//        auto t2=chrono::high_resolution_clock::now();
//        cout << "time : "<< std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " milliseconds\n" << endl;

        
//        t1=chrono::high_resolution_clock::now();
//        Data data(NbSimul, subd,T);
//        data.simule();
//        t2=chrono::high_resolution_clock::now();
//        cout << "time : "<< std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " milliseconds\n" << endl;
//        
        
        
        
        
//        for(auto & i : proc.proc) {
//            cout << i ;
//        }
     
        
        
        //Génère une DataBase
        genereData(NbSimul,subd,T,parts);

        //Fait une quantif
        Quantif3 quant(NbSimul,subd,T, parts);
        quant.fillmat();

        cout << quant.Valeurs[0][0];
        quant.writeStrat();

        // Test Quantif sur le proc controlé
        Test test(NbSimul,subd,T,12,parts);
        auto t1=chrono::high_resolution_clock::now();
        vector<double> res = test.simuleStratOpt(NbSimulTest);
//        writeResults(res, subd, 10000, T,"PL");
        auto t2=chrono::high_resolution_clock::now();
        cout << "time : "<< std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " milliseconds\n";
        ofstream outres("PLoptT"+to_string(T)+"NbSimul"+to_string(NbSimul)+"parts"+to_string(parts)+"subd"+to_string(subd)+".dat");
        for(auto i :res) {
        outres << i << endl;
        }
        cout << "PL moyen " << (double)accumulate(res.begin(), res.end(), 0.)/NbSimulTest << endl;

//
        
        
//        // PlotstratImb
//        Test testImb(NbSimul,subd,T,12,parts);
//        testImb.plotstratImb();
        
        
        //        test.InventaireOpt();
        //        cout << "nbSimuls : " << NbSimul << endl;
        //        cout << "subd :" << subd << endl;
        //        cout <<"Temps :" << T << endl;
        
        
        
        //    ///////////////// Convergence quantif
        //convergenceQuantif(subd, T, parts, NbSimulTest); // PL à maturité T, en suivant les strategies optimales issues de N simulations de trajectoires du carnet d'ordres.  (1<N<1000000)

        //debug
//        Test11 debug11(NbSimul,subd,T,12,parts);
//        vector<double> resbug(1);
//        for (int i=0; i<1000000; i++) {
//            resbug=debug11.simuleStrat11(1);
//        }
        
        
        
            //////////// Test PRocs11
            
                    // Test PL
        Test11 test11(NbSimulTest,subd,T,12,parts);
        double bench =0;
        ofstream fout("PL11"+to_string(NbSimul)+"nbSimulT"+to_string(T)+".dat");
        vector<double> res11(NbSimulTest);
        res11=test11.simule11();
        for(auto i : res11) {
        fout << i << endl;
        }
        bench=1./NbSimulTest*accumulate(res11.begin(), res11.end(), 0.);
        cout << "bench :" << bench << endl;
        
        //        // Test Quantif vs MC
        //        double bench=0;
        //            //genereData11(NbSimul,subd,T);
        //            Test11 test(subd,T);
        //            test.fillmatT(bench);
        //            cout << "nbSimuls : " << NbSimul << endl;
        //            cout << "subd :" << subd << endl;
        //            cout <<"Temps :" << T << endl;
        //            cout << "EsperanceQuantif :" << test.Valeurs[0][0] << endl;
        //            cout << "EsperanceMC :" << bench << endl;
        //    
        //    //////////// Test Random Generator
        //    
        //            Carnet A;
        //            Carnet B;
        //        A=B;
        //        cout << A.test() << endl;
        //            A=B;
        //        cout << A.test() << endl;
        //        assert(A.test()==B.test());
        //            cout << B.test() << " " << A.test() << endl;
        //            cout << B.test() << " " << A.test() << endl;
        //            cout << B.test() << " " << A.test() << endl;
        //            cout << B.test() << " " << A.test() << endl;
        //    
        return 0;
    }
