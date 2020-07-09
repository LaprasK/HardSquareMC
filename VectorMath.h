#ifndef VecMath_Included
#define VecMath_Included

#include <math.h>


template <typename T>
struct vec2{
    vec2(): x(0), y(0) {}
    
    vec2(const T& _x, const T& _y): x(_x), y(_y){}

    vec2(const vec2<T>& other): x(other.x), y(other.y) {}

    T x;
    T y;
};

template<typename T>
inline vec2<T> operator+ (const vec2<T>& a, const vec2<T>&b){
    return vec2<T>(a.x + b.x, a.y+b.y);
}

template<typename T>
inline vec2<T> operator- (const vec2<T>& a, const vec2<T>&b){
    return vec2<T>(a.x - b.x, a.y - b.y);
}

template<typename T>
inline vec2<T> operator* (const vec2<T>& a, const vec2<T>& b){
    return vec2<T>(a.x*b.x, a.y*b.y);
}

template<typename T>
inline vec2<T> operator/ (const vec2<T>& a, const vec2<T>& b){
    return vec2<T>(a.x/b.x, a.y/b.y);
}

template<typename T>
inline vec2<T> operator- (const vec2<T>& a){
    return vec2<T>(-a.x, -a.y);
}

template <typename T>
inline vec2<T>& operator+= (vec2<T>& a, const vec2<T>&b){
    a.x += b.x;
    a.y += b.y;
    return a;
}

template<typename T>
inline vec2<T>& operator-= (vec2<T>& a, const vec2<T>& b){
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

template<typename T>
inline vec2<T>& operator*= (vec2<T>& a, const vec2<T>& b){
    a.x*=b.x;
    a.y*=b.y;
    return a;
}

template<typename T>
inline vec2<T>& operator/= (vec2<T>& a, const vec2<T>& b){
    a.x/=b.x;
    a.y/=b.y;
    return a;
}

template<typename T>
inline vec2<T> operator*(const vec2<T>& a, const T b){
    return vec2<T>(a.x*b, a.y*b);
}

template<typename T>
inline vec2<T> operator*(const T b, const vec2<T>& a){
    return vec2<T>(a.x*b, a.y*b);
}

template<typename T>
inline vec2<T> operator/(const vec2<T>& a, const T b){
    T q = T(1.0)/b;
    return a*q;
}

template<typename T>
inline vec2<T>& operator*= (vec2<T>& a, const T b){
    a.x*=b;
    a.y*=b;
    return a;
}

template<typename T>
inline vec2<T>& operator/= (vec2<T>& a, const T b){
    a.x/=b;
    a.y/=b;
    return a;
}

template<typename T>  
inline bool operator==(const vec2<T>& a, const vec2<T>& b){
    return (a.x==b.x)&&(a.y==b.y);
}

template<typename T>
inline bool operator!=(const vec2<T>& a, const vec2<T>& b){
    return (a.x!=b.x)||(a.y!=b.y);
}

template<typename T>
inline T dot(const vec2<T>& a, const vec2<T>& b){
    return (a.x*b.x + a.y*b.y);
}

template<typename T>
inline vec2<T> perp(const vec2<T>& a){
    return vec2<T>(-a.y, a.x);
}

template<typename T>
inline T perpdot(const vec2<T>& a, const vec2<T>& b){
    return dot(perp(a,b));
}

template<typename T>
inline vec2<T> rotation_matrix(double alpha, const vec2<T>& a){
    double c = cos(alpha), s = sin(alpha);
    return vec2<T>(c*a.x - s*a.y, s*a.x + c*a.y);
}
#endif