#pragma once
#ifndef LIGHT_H
#define LIGHT_H
#include"Vector.hpp"

class Light
{
public:
    Light(Vec3& p, color& i) : position(p), intensity(i) {}
    Light( color& i) : position(Vec3(0,0,0)), intensity(i) {}

    virtual ~Light() = default;
    Vec3 position;
    color intensity;
};
////√Êπ‚‘¥
class AreaLight : public Light
{
public:
    AreaLight(Vec3& p, color& i, Vec3& v0, Vec3 v1, Vec3 v2) : Light(p, i)
    {
        V0 = v0;
        V1 = v1;
        V2 = v2;
        normal = Cal_Normal_3D(v0, v2, v1);
    }
    AreaLight(color& i, Vec3& v0, Vec3 v1, Vec3 v2) : Light(i)
    {
        V0 = v0;
        V1 = v1;
        V2 = v2;
        normal = Cal_Normal_3D(v0, v2, v1);
    }
    double GetArea()
    {
        Vec3 u = V0 - V1;
        Vec3 v = V2 - V0;
        Vec3 n = u.cross_multiply(v);
        return 0.5 * n.len;
    }
    float length;
    Vec3 normal;
    Vec3 V0, V1, V2;
};

#endif