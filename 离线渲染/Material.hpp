#pragma once
#ifndef MATERIAL_H
#define MATERIAL_H
#include"Vector.hpp"
#include <random>
#define M_PI 3.1415926
enum MaterialType {
    DIFFUSE,
    MIRROR,
    Microfacet,
    Transparent
};
float clamp(const float& lo, const float& hi, const float& v)
{
    return std::max(lo, std::min(hi, v));
}
float get_random_float()
{
    static std::random_device dev;
    static std::mt19937 rng(dev());
    static std::uniform_real_distribution<float> dist(0.f, 1.f); // distribution in range [0，1]

    return dist(rng);
}
class Material 
{
public:
    MaterialType m_type;
    //Vector3f m_color;
    color m_emission;
    float ior;
    color Kd, Ks;
    float specularExponent;
    bool haslight = false;
    //Texture tex;

    Material(MaterialType t = MaterialType::DIFFUSE, color e = color(0,0,0)) {
        m_type = t;
        //m_color = c;
        m_emission = e;
    }
    // sample a ray by Material properties
    inline Vec3 sample(Vec3& wi, Vec3& N);
    inline Vec3 sample(Vec3& N);
    // given a ray, calculate the PdF of this ray
    inline double pdf(Vec3& wi, Vec3& wo,  Vec3& N);
    // given a ray, calculate the contribution of this ray
    inline color eval(Vec3& wi, Vec3& wo, Vec3& N);
    inline float DistributionGGX(Vec3& N, Vec3& H, const float& roughness);
    inline float GeometrySub(const float& NdotV, const float& roughness);
    inline float GeometrySmith(Vec3& N, Vec3& V, Vec3& L, const float& roughness);
    // Compute reflection direction
    Vec3 reflect(Vec3& I,  Vec3& N) 
    {
        return I -  N*2 * I.dot_multiply(N);
    }
    // 计算折射方向


    // Compute Fresnel equation
    //
    // \param I is the incident view direction
    //
    // \param N is the normal at the intersection point
    //
    // \param ior is the material refractive index
    //
    // \param[out] kr is the amount of light reflected
    ///菲涅尔方程系数
    void fresnel(Vec3& I, Vec3& N, const double& ior, double& kr)
    {
        float cosi = clamp(-1, 1, I.dot_multiply(N));
        float etai = 1, etat = ior;
        if (cosi > 0) { std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            kr = 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }
    float fresnel(Vec3& I, Vec3& N, const double& ior, double& kr,int i)
    {
        float cosi = clamp(-1, 1, I.dot_multiply(N));
        float etai = 1, etat = ior;
        if (cosi > 0) { std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            return 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            return (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }
    float fresnel(float cos_theta, float eta) 
    {
        float r0 = (1.0f - eta) / (1.0f + eta);
        r0 = r0 * r0; // 反射率的基础值
        return r0 + (1.0f - r0) * std::pow(1.0f - cos_theta, 5); // Schlick近似
    }

    Vec3 toWorld(const Vec3& a, Vec3& N) {
        Vec3 B, C;
        Vec3 resultx;
        if (std::fabs(N.x) > std::fabs(N.y)) {
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C = Vec3(N.z * invLen, 0.0f, -N.x * invLen);
        }
        else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vec3(0.0f, N.z * invLen, -N.y * invLen);
        }
        B = C.cross_multiply(N);
        resultx = B * a.x;
        return B * a.x + C * a.y + N * a.z;
    }
    bool hasEmission() const {
        if (m_emission.r*m_emission.r + m_emission.g * m_emission.g+m_emission.b * m_emission.b > 0.001) return true;
        else return false;
    }
    bool refract( Vec3& I, Vec3& N, double eta, Vec3& refracted) {
        double cosi = clamp(-1.0, 1.0, I.dot_multiply(N));
        double etai = 1.0, etat = eta;
        Vec3 n = N;
        if (cosi < 0) { cosi = -cosi; }
        else { std::swap(etai, etat); n = -N; }
        double etaRatio = etai / etat;
        double k = 1 - etaRatio * etaRatio * (1 - cosi * cosi);
        if (k < 0) {
            return false; // 全内反射
        }
        else {
            refracted = I * etaRatio + n * (etaRatio * cosi - sqrt(k));
            return true;
        }
    }
    Vec3 refract(Vec3& I, Vec3& N, const float& ior)
    {
        float cosi = clamp(-1, 1, I.dot_multiply(N));
        float etai = 1, etat = ior;
        Vec3 n = N;
        if (cosi < 0) { cosi = -cosi; }
        else { std::swap(etai, etat); n = -N; }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        if (k < 0) {
            return Vec3(0, 0, 0); // 当 k < 0 时，返回零向量
        }
        else {
            float G = (eta * cosi - sqrtf(k));
            return I * eta + n * G; // 当 k >= 0 时，计算并返回折射向量
        }

    }


};
float Material::DistributionGGX(Vec3& N, Vec3& H, const float& roughness)
{
    float a = roughness * roughness;
    float a2 = a * a;
    float cosnh = std::fmax(0.0f, N.dot_multiply(H));
    float cosnh2 = cosnh * cosnh;
    float x2 = (cosnh2 * (a2 - 1) + 1) * (cosnh2 * (a2 - 1) + 1);

    return a2 / (M_PI * x2);
}
float Material::GeometrySub(const float& NdotV, const float& roughness)
{
    float r = (roughness + 1.0);
    float k = (r * r) / 8.0;

    float nom = NdotV;
    float denom = NdotV * (1.0 - k) + k;
    return nom / denom;
}

float Material::GeometrySmith(Vec3& N, Vec3& V, Vec3& L, const float& roughness)
{
    float NdotV = std::fmax(N.dot_multiply(V), 0.0f);
    float NdotL = std::fmax(N.dot_multiply(L), 0.0f);
    return GeometrySub(NdotV, roughness) * GeometrySub(NdotL, roughness);
}
Vec3 Material::sample(Vec3& wi, Vec3& N)//入射光线方向和法向，返回出射方向
{
    switch (m_type) {
    case DIFFUSE:
    {
        // uniform sample on the hemisphere
        float x_1 = get_random_float(), x_2 = get_random_float();
        float z = std::fabs(1.0f - 2.0f * x_1);
        float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
        Vec3 localRay(r * std::cos(phi), r * std::sin(phi), z);
        return toWorld(localRay, N);

        break;
    }
    case MIRROR:
    {
        Vec3 reflected = reflect(wi, N);
        return reflect(wi, N);
        break;
    }
    case Microfacet:
    {
        return reflect(wi, N);
        break;
    }
    case Transparent:
    {
        // 使用菲涅尔方程决定反射或折射
        double kr;
        fresnel(wi, N, ior, kr);
        if (get_random_float() < kr) {
            // 反射
            return reflect(wi, N);
        }
        else {
            // 折射
            Vec3 refracted;
            double eta = (wi.dot_multiply(N) < 0) ? (1.0 / ior) : ior;
            if (refract(wi, N, eta, refracted)) {
                return refracted.normalize();
            }
            else {
                // 全内反射
                return reflect(wi, N);
            }
        }
    }
    }
}

double Material::pdf(Vec3& wi, Vec3& wo, Vec3& N) {
    switch (m_type)
    {
    case MaterialType::DIFFUSE:
    {
        // uniform sample probability 1 / (2 * PI)
        if (wo.dot_multiply(N) > 0.0f)
            return 0.5f / M_PI;
        else
            return 0.0f;

        break;
    }
    case MaterialType::MIRROR:
    {
        // PDF for specular reflection is Dirac delta, represented as 1 for computational purposes
        Vec3 reflected = reflect(wi, N);
        if (wo.dot_multiply(N) > 0.001){
            return 1.0f; // Dirac delta is approximated as 1
        }
        else {
            return 0.0f;
        }
        break;
    }
    case Microfacet: {
        Vec3 reflected = reflect(wi, N);
        if (reflected == wo) {
            return 1.0; // Dirac delta is approximated as 1
        }
        else {
            return 0.0;
        }
        break;
    }
    case MaterialType::Transparent:
    {
        double kr;
        fresnel(wi, N, ior, kr);
        Vec3 reflected = reflect(wi, N);
        Vec3 refracted;
        bool hasRefract = refract(wi, N, ior, refracted);
        if (hasRefract) {
            // 反射和折射的pdf
            return kr * 1.0 + (1 - kr) * 1.0; // Dirac delta 近似为1
        }
        else {
            // 只有反射存在
            return 1.0;
        }
    }
    }

}
color Material::eval(Vec3& wi, Vec3& wo, Vec3& N){

    switch (m_type)
    {
    case MaterialType::DIFFUSE:
    {
        float cosalpha = N.dot_multiply(wo);
        if (cosalpha > 0.0f) {
            color diffuse = Kd / M_PI;
            return diffuse;
        }
        else
            return color(0, 0, 0);
        break;
    }
    case MaterialType::MIRROR:
    {
        float cosalpha = N.dot_multiply(wo);
        //printf(" ");
        if (cosalpha > 0.001) {
            double F;
            fresnel(wi, N, ior, F);
            color mirror = { 1 / cosalpha, 1 / cosalpha, 1 / cosalpha };
            return  mirror* F;
        }
        else
            return color(0.0f,0,0);
        break;
    }
    case MaterialType::Transparent:
    {
        double kr;
        fresnel(wi, N, ior, kr);
        if (wi.dot_multiply(N) * wo.dot_multiply(N) > 0) {
            // 反射
            return color(kr, kr, kr);
        }
        else {
            // 折射
            return color(1.0 - kr, 1.0 - kr, 1.0 - kr);
        }
    }
    }
}


#endif