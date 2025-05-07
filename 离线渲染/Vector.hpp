
#pragma once
#ifndef VECTOR_H
#define VECTOR_H
#include <cmath>
#include <algorithm>

class Vec3
{
public:
    double x;
    double y;
    double z;
    double len;
    double ratio;
public:
    Vec3() {
        x = 0;
        y = 0;
        z = 0;
        len = 0;
        ratio = 1;
    }

    Vec3(double ax, double ay, double az)
    {
        x = ax;
        y = ay;
        z = az;
        ratio = 1;
       len = sqrtf(x * x + y * y + z * z);
    }

    Vec3(double ax, double ay, double az, double aratio)
    {
        x = ax;
        y = ay;
        z = az;
        ratio = ratio;
        len = sqrtf(x * x + y * y + z * z);

    }
    //start点指向end点的向量
    Vec3(Vec3 start, Vec3 end)
    {
        x = end.x - start.x;
        y = end.y - start.y;
        z = end.z - start.z;
        ratio = end.ratio - start.ratio;
        len = sqrtf(x * x + y * y + z * z);

    }

    //叉乘 this X input_vector
    Vec3 cross_multiply(const Vec3& input_vector)
    {
        double new_x = y * input_vector.z - z * input_vector.y;
        double new_y = z * input_vector.x - x * input_vector.z;
        double new_z = x * input_vector.y - y * input_vector.x;
        return  Vec3(new_x, new_y, new_z);
    }

    //点乘 this * input_vector
    double dot_multiply(const Vec3& input_vector)
    {
        return x * input_vector.x + y * input_vector.y + z * input_vector.z;
    }

    Vec3 normalize()
    {
        return Vec3(this->x / this->len, this->y / this->len, this->z / this->len);
    }
    Vec3 normalize_ratio()
    {
        if (this->ratio != 0)
            return Vec3(this->x / this->ratio, this->y / this->ratio, this->z, 1);
        return Vec3(this->x, this->y, this->z, this->ratio);
    }

    Vec3 operator*(double n)
    {
        double new_x = x * n;
        double new_y = y * n;
        double new_z = z * n;
        return Vec3(new_x, new_y, new_z);
    }
    Vec3 operator/(double n)
    {
        double new_x = x / n;
        double new_y = y / n;
        double new_z = z / n;
        return Vec3(new_x, new_y, new_z);
    }
    const Vec3& operator-(const Vec3& input_vector)
    {
        double new_x = x - input_vector.x;
        double new_y = y - input_vector.y;
        double new_z = z - input_vector.z;
        return Vec3(new_x, new_y, new_z);
    }

    const Vec3& operator-()
    {
        double new_x = -x;
        double new_y = -y;
        double new_z = -z;
        return Vec3(new_x, new_y, new_z);
    }
    const Vec3& operator-(double& a)
    {
        double new_x = x - a;
        double new_y = y - a;
        double new_z = z - a;
        return Vec3(new_x, new_y, new_z);
    }
    const Vec3& operator+()
    {
        double new_x = +x;
        double new_y = +y;
        double new_z = +z;
        return Vec3(new_x, new_y, new_z);
    }
    Vec3 operator+(const Vec3& input_vector)
    {
        double new_x = x + input_vector.x;
        double new_y = y + input_vector.y;
        double new_z = z + input_vector.z;
        return Vec3(new_x, new_y, new_z);
    }
    bool operator==(Vec3& input_vector)
    {
        if (fabs(input_vector.x - x) < 0.001 && fabs(input_vector.y - y) < 0.001 && fabs(input_vector.z - z) < 0.001)
        {
            return true;
        }
        else
            return false;
    }
};
class color
{
public:
    double r;
    double g;
    double b;
public:
    color() {}

    color(double r, double g, double b)
    {
        this->r = r;
        this->g = g;
        this->b = b;
    }

    color operator+(color rgb)
    {
        color temp;
        temp.r = this->r + rgb.r;
        temp.g = this->g + rgb.g;
        temp.b = this->b + rgb.b;
        return temp;
    }

    color operator-(color rgb)
    {
        color temp;
        temp.r = this->r - rgb.r;
        temp.g = this->g - rgb.g;
        temp.b = this->b - rgb.b;
        return temp;
    }

    color operator*(double rate)
    {
        color temp;
        temp.r = this->r * rate;
        temp.g = this->g * rate;
        temp.b = this->b * rate;
        return temp;
    }
    color operator*(const color& other) const {
        return color(r * other.r, g * other.g, b * other.b);
    }

    color operator/(double rate)
    {
        color temp;
        temp.r = this->r / rate;
        temp.g = this->g / rate;
        temp.b = this->b / rate;
        return temp;
    }

    void operator=(double t)
    {
        this->r = t;
        this->g = t;
        this->b = t;
    }
    color operator*(color& v)
    {
        color temp;
        temp.r = r * v.r;
        temp.g = g * v.g;
        temp.b = b * v.b;
        return temp;
    }
    color operator/(color& v)
    {
        color temp;
        temp.r = r / v.r;
        temp.g = g / v.g;
        temp.b = b / v.b;
        return temp;
    }
};
inline Vec3 Cal_Normal_3D(const Vec3& v1, const Vec3& v2, const Vec3& v3)
{
    //v1(n1,n2,n3);
    //平面方程: na * (x – n1) + nb * (y – n2) + nc * (z – n3) = 0 ;4
    double na = (v2.y - v1.y) * (v3.z - v1.z) - (v2.z - v1.z) * (v3.y - v1.y);
    double nb = (v2.z - v1.z) * (v3.x - v1.x) - (v2.x - v1.x) * (v3.z - v1.z);
    double nc = (v2.x - v1.x) * (v3.y - v1.y) - (v2.y - v1.y) * (v3.x - v1.x);
    Vec3 vn{ na,nb,nc };
    vn = vn.normalize();
    //平面法向量
    return vn;

}


#endif