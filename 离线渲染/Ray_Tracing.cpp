#include <cmath>
#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>
#include<fstream>
#include<sstream>
#include<GL/glut.h>
#include <random>
#include<algorithm>
#include<omp.h>
#include"Light.hpp"
#include"Vector.hpp"
#include"Material.hpp"

#define adress "D://home work//GAMES 101//PA7-1//PA7//Assignment7//models//cornellbox//shortbox.obj"
#define adress1 "D://home work//GAMES 101//PA7-1//PA7//Assignment7//models//cornellbox//floor.obj"
#define adress2 "D://home work//GAMES 101//PA7-1//PA7//Assignment7//models//cornellbox//left.obj"
#define adress3 "D://home work//GAMES 101//PA7-1//PA7//Assignment7//models//cornellbox//right.obj"
#define adress4 "D://home work//GAMES 101//PA7-1//PA7//Assignment7//models//cornellbox//tallbox.obj"
#define adress5 "D://home work//GAMES 101//PA7-1//PA7//Assignment7//models//cornellbox//light.obj"
#define adress6 "D://home work//GAMES 101//PA7-1//PA7//Assignment7//models//bunny//transformed_bunny.obj"

#define clipping_height 500		//�ü��߶�
#define clipping_width 500	//�ü����
#define clipping_polygon_point_count 4
#define IOR 1.5 // ������

int nearplane_width = clipping_width * 2;		//�Ӿ�����
int nearplane_height = clipping_height * 2;		//�Ӿ���߶�
int nearplane_distance = 1000; //�Ӿ����ƽ�����ӵ����
int farplane_distance = nearplane_distance + 2000; //�Ӿ���Զƽ�����ӵ����

float RussianRoulette = 0.9;//���̶ĸ���
int maxdepth = 8;//���ݹ�ֵ
double SSP = 64;//����������


class AABB {
public:
    Vec3 vmin; // ��Χ�е���С��
    Vec3 vmax; // ��Χ�е�����

    // ��ʼ������
    AABB()
    {
        vmin = Vec3(DBL_MAX, DBL_MAX, DBL_MAX);
        vmax = Vec3(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    }
    // ���°�Χ��
    void expand(const Vec3& point) {
        vmin.x = std::min(vmin.x, point.x);
        vmin.y = std::min(vmin.y, point.y);
        vmin.z = std::min(vmin.z, point.z);

        vmax.x = std::max(vmax.x, point.x);
        vmax.y = std::max(vmax.y, point.y);
        vmax.z = std::max(vmax.z, point.z);
    }
};
class Face
{
public:
    std::vector<Vec3> point;
    std::vector<Vec3> vt;
    std::vector<Vec3> vn;
};
class Body
{
public:
    std::vector<Face> faces;
    std::vector<Vec3> v;
    std::vector<Vec3> vt;//����
    std::vector<Vec3> vn;//���㷨����
    AABB aabb;
    Material* material;
    AreaLight* LIGHT;
    bool hasEmit = false;
};

class Ray
{
public:
    Vec3 orig;
    Vec3 dir; // �����������Ա������������ dirPoint
    Vec3 direction;
    double t;
    double tmax, tmin;
    Ray(Vec3 ori, Vec3 direct)
        : orig(ori), dir(direct), direction(ori, direct) // ���� Vec3 ����Ӧ�Ĺ��캯��
    {
        direction = direction.normalize();
        t = 0;
        tmax = -INT_MAX;
        tmin = INT_MAX;
    }
};

class Intersection
{
public:
    Intersection() {
        happened = false;
        coords = Vec3(0,0,0);
        normal = Vec3();
        distance = std::numeric_limits<double>::max();
        obj = nullptr;
        emit = color(0,0,0);
        m = nullptr;
    }
    bool happened;
    Vec3 coords;
    Vec3 tcoords;
    Vec3 normal;
    double distance;
    color emit;
    Body* obj;
    Face* face;
    Material* m;
    AreaLight* l = NULL;
};


int light_count = 2;
Vec3 light_position{ 200, 200, 0 };
Vec3 light_v0 = light_position + Vec3(100, 0, 0);
Vec3 light_v1 = light_position + Vec3(-100, 0, 0);
Vec3 light_v2 = light_position + Vec3(0, 100,0);

Vec3 light1_v0 = Vec3(600,200,100);
Vec3 light1_v1 =Vec3(600, 200, -100);
Vec3 light1_v2 = Vec3(600, 200, 100);

Vec3 light2_v0 = Vec3(343.0, 548.7 ,227.0);
Vec3 light2_v1 = Vec3(343.0, 548.7, 332.0);
Vec3 light2_v2 = Vec3(213.0, 548.7 ,332.0);

Vec3 light3_v0 = Vec3(343.0, 548.7, 227.0);
Vec3 light3_v1 = Vec3(343.0, 548.7, 332.0);
Vec3 light3_v2 = Vec3( 213.0, 548.7, 227.0);

color intensity{ 70,70,70};
color intensity1{ 300,300,300 };
color intensity2{  color(0.747f + 0.058f, 0.747f + 0.258f, 0.747f)*8.0f + color(0.740f + 0.287f,0.740f + 0.160f,0.740f) * 15.6 +  color(0.737f + 0.642f,0.737f + 0.159f,0.737f) * 18.4f};
color intensity3 = intensity2 * 4;
std::vector<AreaLight> lights;

AreaLight light1{ light_position, intensity ,light_v0 ,light_v1 ,light_v2 };
AreaLight light2{ intensity1 ,light1_v0 ,light1_v1 ,light1_v2 };
AreaLight light3{ intensity3 ,light2_v0 ,light2_v2 ,light2_v1 };
AreaLight light4{ intensity3 ,light3_v0 ,light3_v2 ,light3_v1 };
Body obj1[7];
int obj_count = 5;
Vec3 Camera{ 260, 250, -400 };

color colorbuffer[clipping_width][clipping_width];

void load_obj(const char* filename, Body* obj,Material *m)
{
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error" << std::endl;
        exit(1);
    }
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "v")
        {
            Vec3 temp;
            iss >> temp.x >> temp.y >> temp.z;
            temp.ratio = 1;
            obj->v.push_back(temp);

            ///����AABB
            obj->aabb.expand(temp);
        }
        else if (type == "vn")
        {
            Vec3 temp;
            iss >> temp.x >> temp.y >> temp.z;
            obj->vn.push_back(temp);
        }
        else if (type == "vt")
        {
            Vec3 temp;
            iss >> temp.x >> temp.y;
            obj->vt.push_back(temp);
        }
        else if (type == "f")
        {
            int temp[3];
            Face Temp;
            char slash;
            for (int i = 0; i < 3; i++)
            {
                //iss >> temp[0] >> slash >> temp[1] >> slash >> temp[2];
                iss >> temp[0] >>  temp[1] >>  temp[2];
                Temp.point.push_back(obj->v[temp[i] - 1]);
                if (obj->vn.size() != 0)
                {
                    Temp.vt.push_back(obj->vt[temp[1] - 1]);
                }
                if (obj->vn.size() != 0)
                {
                    Temp.vn.push_back(obj->vn[temp[2] - 1]);
                }
            }
            obj->faces.push_back(Temp);
        }
    }

    obj->material = m;
}

void load_obj(const char* filename, Body* obj, Material* m,AreaLight* l)
{
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error" << std::endl;
        exit(1);
    }
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "v")
        {
            Vec3 temp;
            iss >> temp.x >> temp.y >> temp.z;
            temp.ratio = 1;
            obj->v.push_back(temp);

            ///����AABB
            obj->aabb.expand(temp);
        }
        else if (type == "vn")
        {
            Vec3 temp;
            iss >> temp.x >> temp.y >> temp.z;
            obj->vn.push_back(temp);
        }
        else if (type == "vt")
        {
            Vec3 temp;
            iss >> temp.x >> temp.y;
            obj->vt.push_back(temp);
        }
        else if (type == "f")
        {
            int temp[3];
            Face Temp;
            char slash;
            for (int i = 0; i < 3; i++)
            {
                //iss >> temp[0] >> slash >> temp[1] >> slash >> temp[2];
                iss >> temp[0] >> temp[1] >> temp[2];
                Temp.point.push_back(obj->v[temp[i] - 1]);
                if (obj->vn.size() != 0)
                {
                    Temp.vt.push_back(obj->vt[temp[1] - 1]);
                }
                if (obj->vn.size() != 0)
                {
                    Temp.vn.push_back(obj->vn[temp[2] - 1]);
                }
            }
            obj->faces.push_back(Temp);
        }
    }

    obj->material = m;
    obj->LIGHT = l;
}
///

/////�Ƿ����Χ���ཻ
bool rayIntersectsAABB(Ray& ray, const AABB& box) {
    double tmin = -DBL_MAX;
    double tmax = DBL_MAX;

    // x �����
    if (ray.direction.x != 0) {
        double invDx = 1.0 / ray.direction.x;
        double tx0 = (box.vmin.x - ray.orig.x) * invDx;
        double tx1 = (box.vmax.x - ray.orig.x) * invDx;
        if (invDx < 0.0) std::swap(tx0, tx1);
        tmin = std::max(tmin, tx0);
        tmax = std::min(tmax, tx1);
    }
    else if (ray.orig.x < box.vmin.x || ray.orig.x > box.vmax.x) {
        return false; // �������ƽ����ԭ�㲻�ڷ�Χ�ڣ����ཻ
    }

    // y �����
    if (ray.direction.y != 0) {
        double invDy = 1.0 / ray.direction.y;
        double ty0 = (box.vmin.y - ray.orig.y) * invDy;
        double ty1 = (box.vmax.y - ray.orig.y) * invDy;
        if (invDy < 0.0) std::swap(ty0, ty1);
        tmin = std::max(tmin, ty0);
        tmax = std::min(tmax, ty1);
    }
    else if (ray.orig.y < box.vmin.y || ray.orig.y > box.vmax.y) {
        return false; // �������ƽ����ԭ�㲻�ڷ�Χ�ڣ����ཻ
    }

    // z �����
    if (ray.direction.z != 0) {
        double invDz = 1.0 / ray.direction.z;
        double tz0 = (box.vmin.z - ray.orig.z) * invDz;
        double tz1 = (box.vmax.z - ray.orig.z) * invDz;
        if (invDz < 0.0) std::swap(tz0, tz1);
        tmin = std::max(tmin, tz0);
        tmax = std::min(tmax, tz1);
    }
    else if (ray.orig.z < box.vmin.z || ray.orig.z > box.vmax.z) {
        return false; // �������ƽ����ԭ�㲻�ڷ�Χ�ڣ����ཻ
    }

    // �ж��Ƿ��ཻ
    return tmin <= tmax && tmax>0;
}

////�Ƿ�����������ƽ���ཻ������������ڸ������ε���������,dir��direction
bool rayTriangleIntersect(Vec3& v0, Vec3& v1, Vec3& v2,  Vec3& orig,
    Vec3& dir, float& tnear, float& u, float& v)
{
    bool isIn = false;
    Vec3 E1 = v1 - v0;
    Vec3 E2 = v2 - v0;
    Vec3 S = orig - v0;
    Vec3 S1 = dir.cross_multiply(E2);
    Vec3 S2 = S.cross_multiply(E1);
    float coeff = 1.0 / S1.dot_multiply(E1); // ��ͬϵ��
    float t = coeff * S2.dot_multiply(E2);
    float b1 = coeff * S1.dot_multiply(S);
    float b2 = coeff * S2.dot_multiply(dir);
    if (t >= 0 && b1 >= 0 && b2 >= 0 && (1 - b1 - b2) >= 0)
    {
        isIn = true;
        tnear = t;
        u = b1;
        v = b2;
    }
    return isIn;
}

Vec3 rayTriangleIntersect(Vec3& v0, Vec3& v1, Vec3& v2, Vec3& orig,
    Vec3& dir)
{
    bool isIn = false;
    Vec3 intersection;
    Vec3 E1 = v1 - v0;
    Vec3 E2 = v2 - v0;
    Vec3 S = orig - v0;
    Vec3 S1 = dir.cross_multiply(E2);
    Vec3 S2 = S.cross_multiply(E1);
    float coeff = 1.0 / S1.dot_multiply(E1); // ��ͬϵ��
    float t = coeff * S2.dot_multiply(E2);
    float b1 = coeff * S1.dot_multiply(S);
    float b2 = coeff * S2.dot_multiply(dir);
    if (t >= 0 && b1 >= 0 && b2 >= 0 && (1 - b1 - b2) >= 0)
    {
        isIn = true;
        intersection = orig + dir * t;
    }
    return intersection;
}

Intersection get_intersection(Ray& ray)
{
    Intersection inter;
    inter.happened = false;
    for (int i = 0; i < obj_count; i++)
    {
        if (rayIntersectsAABB(ray, obj1[i].aabb))
        {
            for (int j = 0; j < obj1[i].faces.size(); j++)
            {
                Vec3 v0 = obj1[i].faces[j].point[0];
                Vec3 v1 = obj1[i].faces[j].point[1];
                Vec3 v2 = obj1[i].faces[j].point[2];
                float tnear, u, v;
                tnear = 0;
                if (rayTriangleIntersect(v0, v2, v1, ray.orig, ray.direction, tnear, u, v))
                {
                    if (tnear > 0 && tnear < ray.tmin)
                    {
                        inter.happened = true;
                        ray.tmin = tnear;
                        inter.normal = Cal_Normal_3D(v0, v1, v2).normalize();
                        inter.obj = &obj1[i];
                        inter.face = &obj1[i].faces[j];
                        inter.m = inter.obj->material;
                        inter.coords = ray.orig;
                        Vec3 temp = (ray.direction * ray.tmin);
                        inter.coords = inter.coords + temp;
                        inter.emit = inter.m->m_emission;
                        //if (inter.m->hasEmission())
                        //{
                        //    inter.l = &light1;
                        //}
                    }
                }
            }
        }
    }

    return inter;
}

double get_distance(Vec3& v0,Vec3 &v1)
{
    
    return std::sqrt(powf(v0.x - v1.x, 2) + powf(v0.y - v1.y, 2) + powf(v0.z - v1.z, 2));

}
////�Է���Ĳ��������ڽ��㴦����wo�������һ�����ߣ�Ϊ�����wi,���ﰴ�����������

void Sample_Light(Intersection& pos, float& pdf,int count)
{
    // uniform sample on the hemisphere
    double emit = lights[count].GetArea();
    float x = std::sqrt(get_random_float()), y = get_random_float();//�ڹ�Դ�����ȡһ�㣬��Ϊ������
    pos.coords = lights[count].V0 * (1.0f - x) + lights[count].V1 * (x * (1.0f - y)) + lights[count].V2 * (x * y);
    pdf = 1 / emit;
    float theta = 2.0 * M_PI * get_random_float(), phi = M_PI * get_random_float();
    //dir -> {cos��,sin��*cos��,sin��*sin��}
    Vec3 dir(std::cos(phi), std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta));
    pos.normal = lights[count].normal;///dir.normalize(); //��������ߣ�
    pos.emit = lights[count].intensity;
    pos.l = &lights[count];
}

////�Թ�Դ���в���,��λ������
void sampleLight(Intersection& pos, float& pdf)
{
    ///���й�Դ
    float emit_area_sum = 0;

    for (int i = 0; i < light_count; i++)
    {
        emit_area_sum += lights[i].GetArea();
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (int i = 0; i < light_count; i++)
    {
        emit_area_sum +=lights[i].GetArea();
        if (p <= emit_area_sum)
        {
            Sample_Light(pos, pdf, i);
            return;
        }
    }
}
//������غ���
Vec3 reflect(Vec3& I,  Vec3& N)
{
    return I - N *2 * I.dot_multiply(N);
}
// �������䷽�򣬷����Ƿ�ɹ�����
bool refract(Vec3& I,Vec3& N, float eta_i, float eta_t, Vec3& refracted)
{
    float cosi = clamp(-1.0f, 1.0f, I.dot_multiply(N));
    float eta = eta_i / eta_t;
    float k = 1.0f - eta * eta * (1.0f - cosi * cosi);
    if (k < 0)
        return false; // ȫ�ڷ���
    refracted = I * eta + N * (eta * cosi - sqrt(k));
    return true;
}

// ʹ��Schlick���Ƽ������������ϵ��
float fresnel_schlick(Vec3& I,  Vec3& N, float eta_i, float eta_t)
{
    float cosi = clamp(-1.0f, 1.0f, I.dot_multiply(N));
    float R0 = powf((eta_i - eta_t) / (eta_i + eta_t), 2.0f);
    return R0 + (1 - R0) * powf(1 - fabsf(cosi), 5.0f);
}

color RayTracing(Ray& ray,int depth)
{
    //ֱ�ӹ��ղ���

    color L_dir(0, 0, 0);
    color L_indir(0, 0, 0);
    Intersection inter_obj;

    //inter_obj = get_intersection(ray);
    inter_obj = get_intersection(ray);
    Vec3 p = inter_obj.coords;//pΪ��������
    //if (!inter_obj.happened)
    //{
    //    // ���ñ�����ɫΪ�����
    //    return color(0.5f, 0.7f, 1.0f);
    //}
    if (!inter_obj.happened)//���û���������е����壬����0
       return L_dir + L_indir;
    if (inter_obj.m->hasEmission())
    {
        return inter_obj.m->m_emission;
    }

    switch (inter_obj.m->m_type)
    {
    case MaterialType::DIFFUSE: {
        //���Ч�ʣ��ӹ�Դ�������������㷴����ߺ͹�Դ�ཻ�����
        Intersection inter_light;
        float pdf_light;
        sampleLight(inter_light, pdf_light);
        Vec3 light_x = inter_light.coords;//������Դλ��
        Vec3 ws = (light_x - p);//��ɫ��ָ���Դ�ķ��򣬼�����ⷽ��
        ws = ws.normalize();
        Vec3 N = inter_obj.normal;//��ɫ�㷨��
        Vec3 wo = ray.direction;//������߷���
        //����ɫ��������ߣ�����ɫ�㷴��Ĺ��߱Ƚϣ��ж���ɫ�㷴��Ĺ����ܷ��䵽��Դ
        Ray objTOlight(p, light_x);//�������

        double d1 = get_distance(light_x, p);//��ɫ��͹�Դ�ľ���
        Intersection objtolight;
        objtolight = get_intersection(objTOlight);//���ֱ�ӹ��յķ�����Ƿ���ֱ���䵽��Դ
        Vec3 lightintersection;
        lightintersection = rayTriangleIntersect(inter_light.l->V0, inter_light.l->V1, inter_light.l->V2, objTOlight.orig, objTOlight.direction);
        double d2 = get_distance(p, objtolight.coords);

        ///����L_dir,ֱ�ӹ���
        if (d2 - d1 > -0.001)///���������ж��Ƿ����ڵ���
        {
            color Eval = inter_obj.m->eval(wo, ws, N);//������ߵ�ǿ��
            color emit = inter_light.emit;// ��ԴͶ�䵽��ɫ��Ĺ��ߵ�radiance������ֻ��һ����Դlight1
            //emit = light1.intensity;
            Vec3 NN = inter_light.normal;
            double cos_theta_obj = N.dot_multiply(ws);
            double cos_theta_light = NN.dot_multiply(-ws);
            double e = std::pow(d1, 2);

            L_dir = emit * Eval * cos_theta_obj * cos_theta_light / std::pow(d1, 2) / pdf_light;
            //printf(" ");
        }
        ////�ݹ�����ӹ��� 

        float P_RR = get_random_float(); // ���ö���˹���̶ģ���������ݹ飬������ѧ����һ��

        if (P_RR < RussianRoulette && depth < maxdepth)
        {
            Vec3 wi = inter_obj.m->sample(wo, N).normalize();//����ɫ����������������
            Vec3 dired = wi + p;
            Ray ray_objTOobj(p, dired);
            Intersection nextObjInter = get_intersection(ray_objTOobj);

            if (nextObjInter.happened == true && !nextObjInter.m->hasEmission())
            {
                color EVal = inter_obj.m->eval(wo, wi, N);//������ɫ��eval
                double Pdf = inter_obj.m->pdf(wo, wi, N);

                float cos_theta = wi.dot_multiply(N);
                if (Pdf > 0.001)
                {
                    Ray a(p, dired);
                    L_indir = RayTracing(a, depth + 1) * EVal * cos_theta / Pdf / RussianRoulette;

                }
            }
        }
        break;
    }
    case MaterialType::MIRROR: {
        float P_RR = get_random_float(); // ���ö���˹���̶ģ���������ݹ飬������ѧ����һ��
        if (P_RR > RussianRoulette)
        {
            return L_dir;
        }

        if (P_RR < RussianRoulette && depth < maxdepth)
        {
            Vec3 wi = inter_obj.m->sample(ray.direction, inter_obj.normal).normalize();//����ɫ����������������
            Vec3 dired = wi + p;
            Ray ray_objTOobj(p, dired);
            Intersection nextObjInter = get_intersection(ray_objTOobj);

            if (nextObjInter.happened == true)// && !nextObjInter.m->hasEmission())
            {
                double Pdf = inter_obj.m->pdf(ray.direction, wi, inter_obj.normal);
                color EVal = inter_obj.m->eval(ray.direction, wi, inter_obj.normal);//������ɫ��eval

                float cos_theta = wi.dot_multiply(inter_obj.normal);
                if (Pdf > 0.001)
                {
                    Ray a(p, dired);
                    L_indir = RayTracing(a, depth + 1) * EVal * cos_theta / Pdf / RussianRoulette;

                }
            }
        }
        break;
    }
    case MaterialType::Transparent:
    {
        color reflection_color(0, 0, 0);
        color refraction_color(0, 0, 0);
        Vec3 N = inter_obj.normal;
        // �жϹ����ǽ�����ʻ����뿪����
        bool entering = ray.direction.dot_multiply(N) < 0;
        Vec3 adjusted_N = N;
        float eta_i = 1.0f;       // ������ʵ������ʣ�Ĭ�Ͽ�����
        float eta_t = inter_obj.m->ior; // Ŀ����ʵ������ʣ����ʵ������ʣ�

        if (!entering)
        {
            // ���ߴӲ����ڲ��뿪������������
            std::swap(eta_i, eta_t);
            adjusted_N = -N;
        }

        // ���㷴�䷽��
        Vec3 reflection_dir = reflect(ray.direction, adjusted_N).normalize();

        // �������䷽��ʹ��˹��������
        Vec3 refraction_dir;
        bool has_refraction = refract(ray.direction, adjusted_N, eta_i, eta_t, refraction_dir);
        // ʹ��Schlick���Ƽ������������ϵ��
        float kr = fresnel_schlick(ray.direction, adjusted_N, eta_i, eta_t);

        // ����˹���̶ĸ���
        float P_RR = get_random_float();

        // �����µ���㣬�������ཻ
        Vec3 new_origin = p + adjusted_N * 1e-4f;

        if (P_RR < kr && depth < maxdepth)
        {
            // ׷�ٷ������
            Ray reflected_ray(new_origin, reflection_dir);
            reflection_color = RayTracing(reflected_ray, depth + 1) * kr / RussianRoulette;
        }
        if (has_refraction && depth < maxdepth)
        {
            // ׷���������
            Ray refracted_ray(new_origin, refraction_dir);

            refraction_color = RayTracing(refracted_ray, depth + 1) * (1.0f - kr) / RussianRoulette;
        }

        // ��Ϸ��������Ĺ���
        L_indir = L_indir + reflection_color + refraction_color;

        break;
    }
    }
    return L_dir + L_indir;
}
void drawTransparentSphere(float radius, int slices, int stacks, float alpha) {
    for (int i = 0; i < stacks; ++i) {
        float theta1 = i * M_PI / stacks;         // ��ǰ��ջ�ĽǶ�
        float theta2 = (i + 1) * M_PI / stacks;   // ��һ��ջ�ĽǶ�

        glBegin(GL_QUAD_STRIP); // ��������
        for (int j = 0; j <= slices; ++j) {
            float phi = j * 2 * M_PI / slices; // ��ǰ��Ƭ�ĽǶ�

            // ���� 1
            float x1 = radius * sin(theta1) * cos(phi);
            float y1 = radius * cos(theta1);
            float z1 = radius * sin(theta1) * sin(phi);
            glColor4f(0.2f, 0.6f, 1.0f, alpha); // ��ɫ͸��
            glVertex3f(x1, y1, z1);

            // ���� 2
            float x2 = radius * sin(theta2) * cos(phi);
            float y2 = radius * cos(theta2);
            float z2 = radius * sin(theta2) * sin(phi);
            glColor4f(0.2f, 0.6f, 1.0f, alpha); // ��ɫ͸��
            glVertex3f(x2, y2, z2);
        }
        glEnd();
    }
}
//���Ʋ���
float getRandomOffset(){
    static thread_local std::mt19937 generator(std::random_device{}());
    std::uniform_real_distribution<float> distribution(-0.5f, 0.5f);
    return distribution(generator);
}
void Pixelcolor(int x, int y, color& pixelcolor)
{
    Vec3 directed;
    directed.x = x;
    directed.y = y;
    directed.z = 0;

    // ����������
    color accumulatedColor(0, 0, 0); // �ۼ���ɫ

    for (int i = 0; i <SSP; i++) {
        float dx = getRandomOffset();
        float dy = getRandomOffset();
        directed.x = x + dx;
        directed.y = y + dy;
        Ray ray(Camera, directed);

        // ���ù���׷�ٺ��������г�������һ��
        color pixel = RayTracing(ray, 0) / SSP;
        accumulatedColor = accumulatedColor + pixel;
    }

    // ��������������ɫ���浽 framebuffer
    pixelcolor = accumulatedColor;
}
void Display( )
{
    glClearColor(1.f, 1.f, 1.f, 0.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBegin(GL_POINTS);
    for (int x = 0; x < clipping_height; ++x) {
        for (int y = 0; y < clipping_width; ++y) {

            glColor3f(colorbuffer[x][y].r, colorbuffer[x][y].g, colorbuffer[x][y].b);
            glVertex2f(clipping_height - x, y);
        }

    }
    glEnd();
    glutSwapBuffers();

}
// ����ͼ��
void Ca() {
    //omp_set_num_threads(15);
#pragma omp parallel for collapse(2)
    for (int x = 0; x < clipping_height; ++x) {
        for (int y = 0; y < clipping_width; ++y) {
            Pixelcolor(x, y, colorbuffer[x][y]);
        }
    }

}

void reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (w <= h)
        glOrtho(0, 0.5 * nearplane_width, 0, 0.5 * nearplane_height * (GLdouble)nearplane_height / (GLdouble)nearplane_width,
            nearplane_distance, farplane_distance); //������ӵ�
    else
        glOrtho(0, 0.5 * nearplane_width, 0, 0.5 * nearplane_height * (GLdouble)nearplane_width / (GLdouble)nearplane_height,
            nearplane_distance, farplane_distance);


    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(0, 0, nearplane_distance, 0, 0, 0, 0, 1, 0);
}

void keyboard(unsigned char key, int x, int y){
    switch (key)
    {
    case 'w':
    case 'W':
    {
        Camera.y += 10;
        glutPostRedisplay();
        break;
    }
    case 's':
    case 'S':
    {
        Camera.y -= 10;
        glutPostRedisplay();
        break;
    }
    case 'a':
    case 'A':
    {
        Camera.x += 10;
        glutPostRedisplay();
        break;
    }
    case 'd':
    case 'D':
    {
        Camera.x -= 10;
        glutPostRedisplay();
        break;
    }
    case 'n':
    case 'N':
    {
        Camera.z -= 10;
        glutPostRedisplay();
        break;
    }
    case 'f':
    case 'F':
    {
        Camera.z += 10;
        glutPostRedisplay();
        break;
    }
    }
}

int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(clipping_width, clipping_height);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Transparent Material Renderer");
    Material* red = new Material(MaterialType::DIFFUSE, color(0,0,0));
    red->Kd = color(0.63f, 0.065f, 0.05f);
    Material* green = new Material(MaterialType::DIFFUSE, color(0, 0, 0));
    green->Kd = color(0.14f, 0.45f, 0.091f);
    Material* white = new Material(MaterialType::DIFFUSE, color(0, 0, 0));
    white->Kd = color(0.725f, 0.71f, 0.68f);
    Material* light = new Material(MaterialType::DIFFUSE, (color(0.747f + 0.058f, 0.747f + 0.258f, 0.747f) * 8.0 + color(0.740f + 0.287f, 0.740f + 0.160f, 0.740f)* 15.6f + color(0.737f + 0.642f, 0.737f + 0.159f, 0.737f) * 18.4f));
    light->Kd = color(0.65f, 0.65f, 0.65f);
    Material* mirror = new Material(MaterialType::MIRROR, color(0, 0, 0));
    mirror->ior = 12.85;
    mirror->Ks = color(0.45, 0.45, 0.45);
    mirror->Kd = color(0.3, 0.3, 0.3);
    //sred->Kd = color(0.3, 0.3, 0.3);
    Material* glass = new Material(MaterialType::Transparent, color(0, 0, 0));
    glass->ior = 1.5;
    glass->Kd = color(0,0,0);
    glass->Ks = color(1, 1, 1);

    load_obj(adress, &obj1[0], white);
    load_obj(adress1, &obj1[1], white);
    load_obj(adress2, &obj1[2], red);
    load_obj(adress3, &obj1[3], green);
    load_obj(adress4, &obj1[4], white);
    load_obj(adress5, &obj1[5], light, &light1);
    //load_obj(adress6, &obj1[4], glass);
    //lights.push_back(light1);
    //lights.push_back(light2);
    
    lights.push_back(light3);
    lights.push_back(light4);
    Ca();
    glutReshapeFunc(reshape);
    glutDisplayFunc(Display); // ��ʾͼ��
    glutKeyboardFunc(keyboard);
    glutMainLoop();
    return 0;
}

