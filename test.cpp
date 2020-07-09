#include <iostream>
#include <vector>
#include "VectorMath.h"
#include "geometry.h"
#include "myrandom.h"
#include <stack>
#include <algorithm>
#include <time.h>
#include <iomanip>
#include <numeric>


using namespace std;

const unsigned int INVALID_NODE = 0xffffffff;

struct Point{
    Point():x(1), y(2){}
    double x;
    double y;
};

int main(){
    vec2<double> a(3,4);
    vec2<double> b(3,4);
    cout << dot(a,b) << endl;
    
    /*
    double sidelength = sqrt(10/0.25);
    Square sq1(vector<double>({0.40571, 0.892707}), 1.86641);
    vector<vec2<double> > verts = sq1.get_vertices();
    for(auto vert: verts){
        std::cout << vert.x << " " << vert.y << endl;
    }
    Square sq2(vector<double>({1.02252, 1.70252}),1.96455);
    vec2<double> box(sidelength, sidelength);
    verts = sq2.get_vertices();
    for(auto vert: verts){
        std::cout << vert.x << " " << vert.y << endl;
    }
    std::cout << test_polygon_overlap(sq1, sq2, box) << endl;
    vector<double> pos(2, 0);
    Square temp;
    temp = Square(pos, M_PI_4);
    
    cout << temp.get_vertices()[0].x << endl;
    cout << temp.orientation << endl;
    

    Rand m;
    vector<double> list;
    double mu = 1;
    double std = 3.0;
    for(int i = 0; i < 1e7; ++i)
        list.push_back(m.normal_rdnm(mu, std));
    double sum = accumulate(list.begin(), list.end(), 0.0);
    double mean = sum / 1e7;
    double sq_sum = inner_product(list.begin(), list.end(), list.begin(), 0.0);
    double stdev = sqrt(sq_sum/1e7 - mean*mean);
    cout << stdev<< endl;
    vec2<double> pos1(9.5,9.5);
    vec2<double> pos2(9.5,0.5);
    Square sq1(pos1, 0);
    Square sq2(pos2, 0);
    vec2<double> box(10,10);
    cout << test_polygon_overlap(sq1, sq2, box) << endl;*/
}
