#ifndef MyRandom_Included
#define MyRandom_Included


#include <time.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <deque>
#include <limits>


class Rand{
public:
    Rand();
    double rdnm1();
    double rdnm();
    double normal_rdnm(double, double&);
    void test_correct();
    void test_average();

private: 
    int i0, length;
    std::vector<double> list;
};

#endif