#ifndef Geometry_Included
#define Geometry_Included

#include "VectorMath.h"
#include <vector>
#include <string>
#include <cmath>


struct Square{
    Square(): position(vec2<double>(0,0)), orientation(0){
        vertices.push_back(vec2<double>(0, sqrt(0.5)));
        vertices.push_back(vec2<double>(sqrt(0.5), 0));
        vertices.push_back(vec2<double>(0, -sqrt(0.5)));
        vertices.push_back(vec2<double>(-sqrt(0.5), 0));       
    }

    Square(const std::vector<double>& _pos, const double _ori): position(vec2<double>(_pos[0],_pos[1])), orientation(_ori){
        vertices.push_back(vec2<double>(0, sqrt(0.5)));
        vertices.push_back(vec2<double>(sqrt(0.5), 0));
        vertices.push_back(vec2<double>(0, -sqrt(0.5)));
        vertices.push_back(vec2<double>(-sqrt(0.5), 0));
    }

    Square(const std::vector<double>& _pos, const double _ori, std::string pt_type): position(vec2<double>(_pos[0], _pos[1])), orientation(_ori){
        if(pt_type=="diag"){
            vertices.push_back(vec2<double>(0, sqrt(0.5)));
            vertices.push_back(vec2<double>(sqrt(0.5), 0));
            vertices.push_back(vec2<double>(0, -sqrt(0.5)));
            vertices.push_back(vec2<double>(-sqrt(0.5), 0));
        }
        else if(pt_type == "edge"){
            vertices.push_back(vec2<double>(-0.5, 0.5));
            vertices.push_back(vec2<double>( 0.5, 0.5));
            vertices.push_back(vec2<double>(0.5, -0.5));
            vertices.push_back(vec2<double>(-0.5,-0.5));
        }
    }

    std::vector<vec2<double> > get_vertices(){
        std::vector<vec2<double> > res;
        for(auto vert: vertices){
            res.push_back(rotation_matrix(orientation, vert) + position);
        }
        return res;
    }

    std::vector<vec2<double> > vertices; 
    vec2<double> position;
    double orientation;
};


inline bool is_inside(const vec2<double>& point, const std::vector<vec2<double> >& vertices){
    int N = vertices.size();
    int i, j = N-1;
    bool oddNodes = false;
    for(i = 0; i < N; ++i){
        if((vertices[i].y < point.y && vertices[j].y >= point.y) || (vertices[i].y >= point.y && vertices[j].y < point.y)){
            if(vertices[i].x + (point.y - vertices[i].y)/(vertices[j].y - vertices[i].y)* (vertices[j].x - vertices[i].x) < point.x){
                oddNodes = !oddNodes;
            }
        }
        j = i;
    }
    return oddNodes;
}


inline unsigned int tri_orientation(const vec2<double>& a, const vec2<double>& b, const vec2<double>& c){
    const double precision = 1e-6;
    double v = ((c.y - a.y)*(b.x - a.x) - (b.y - a.y)*(c.x - a.x));

    if(fabs(v) < precision) return 0;
    else if(v > 0) return 1;
    else return 2;
}


inline bool segment_intersect(const vec2<double>& a, const vec2<double>& b, const vec2<double>& c, const vec2<double>& d){
    unsigned int r1 = tri_orientation(a, c, d);
    unsigned int r2 = tri_orientation(b, c, d);
    unsigned int r3 = tri_orientation(a, b, c);
    unsigned int r4 = tri_orientation(a, b, d);

    if(r1 != r2 && r3 != r4) return true;

    if(r1 == 0 && r2 == 0 && r3==0 && r4==0){
        vec2<double> v = b-a;
        double p1 = dot(a,v);
        double p2 = dot(b,v);
        double min_1 = std::min(p1,p2);
        double max_1 = std::max(p1,p2);

        double p3 = dot(c,v);
        double p4 = dot(d,v);

        if((p3 > min_1 &&  p3 < max_1)||(p4 > min_1 && p4 < max_1))
            return true;
    }

    return false;
}

inline bool test_polygon_overlap(const Square& sq1, const Square& sq2, const vec2<double>& boxSize){
    vec2<double> dr = sq2.position - sq1.position;
    if(dr.x < -0.5*boxSize.x) dr.x += boxSize.x;
    else if(dr.x >= 0.5*boxSize.x) dr.x -= boxSize.x;
    if(dr.y < -0.5*boxSize.y) dr.y += boxSize.y;
    else if(dr.y >= 0.5*boxSize.y) dr.y -= boxSize.y;

    double ab_rotation = -remainder(sq2.orientation, 2*M_PI) + remainder(sq1.orientation, 2*M_PI);

    vec2<double> ab_t = rotation_matrix<double>(-sq2.orientation, dr);

    int j = sq1.vertices.size() - 1;

    vec2<double> prev_a = rotation_matrix<double>(ab_rotation, sq1.vertices[j]) - ab_t;

    for(int i = 0; i < sq1.vertices.size(); ++i){
        vec2<double> cur_a = rotation_matrix(ab_rotation, sq1.vertices[i]) - ab_t;

        int k = sq2.vertices.size() - 1;
        vec2<double> prev_b = sq2.vertices[k];

        for(int l=0; l < sq2.vertices.size(); ++l){
            vec2<double> cur_b = sq2.vertices[l];
            if(segment_intersect(prev_a, cur_a, prev_b, cur_b))
                return true;
            k = l;
            prev_b = cur_b;
        }

        if(is_inside(cur_a, sq2.vertices))
            return true;

        j = i;
        prev_a = cur_a;
    }

    ab_rotation = -ab_rotation;

    ab_t = rotation_matrix(-sq1.orientation, -dr);

    for(int i=0; i < sq2.vertices.size(); ++i){
        vec2<double> cur = rotation_matrix(ab_rotation, sq2.vertices[i]) - ab_t;
        if(is_inside(cur, sq1.vertices))
            return true;
    }

    return false;
}


#endif