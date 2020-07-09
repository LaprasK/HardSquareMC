#include "myrandom.h"

Rand::Rand(){
    length = 100;
    i0 = time(NULL);
    for(int i = 0; i < (length+1); ++i)
        list.push_back(Rand::rdnm1());
}


double Rand::rdnm1(){
    int a = 16807;
    int q = INT32_MAX / a;
    int p = INT32_MAX % a;
    int x = i0/q, y = i0%q;
    i0 = a * y - p * x;
    if(i0 <0) i0 += INT32_MAX;
    return double(i0)/INT32_MAX;
}


double Rand::rdnm(){
    double ret = list[int(list[length]*length)];
    list[length] = ret;
    list[int(list[length]*length)]=rdnm1();
    return ret;
}


double Rand::normal_rdnm(double mu, double& sigma){
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    double u1, u2;
    do{
        u1 = rdnm();
        u2 = rdnm();
    }
    while(u1 <= epsilon);

    double z0 = sqrt(-2.0*log(u1))*cos(M_PI*2*u2);
    return z0*sigma + mu;
}

void Rand::test_average(){
    int R = 1e7;
    for(int n = 1; n < 5; ++n){
        double shu = 0;
        for(int i = 0; i < R; ++i)
            shu += pow(rdnm(), n);
        std::cout << shu/R << std::endl;
    }
    return;
}


void Rand::test_correct(){
    int R = 1e7;
    for(int k=1; k < 5;++k){
        double shu = 0;
        std::deque<double> lstshu;
        for(int j = 0; j < k+1; ++j)
            lstshu.push_back(rdnm());
        for(int i = 0; i < R; ++i){
            shu += (lstshu[0] - 0.5)*(lstshu[k] - 0.5);
            lstshu.push_back(rdnm());
            lstshu.pop_front();
        }
        std::cout << shu/R << std::endl;
    }
    return;
}