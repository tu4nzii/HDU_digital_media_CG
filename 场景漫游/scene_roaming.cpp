#define STB_IMAGE_IMPLEMENTATION
#define GLUT_DISABLE_ATEXIT_HACK
#define _CRT_SECURE_NO_WARNINGS

#include <GLFW/glfw3.h>
#include "stb_image.h"
#include <cmath>
#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>
#include<fstream>
#include<sstream>
#include "GLUT.H"

#define PI 3.1415926535
#define clipping_height 800		//裁剪高度
#define clipping_width 800	 //裁剪宽度
#define clipping_polygon_point_count 4  //裁剪多边形边数
#define lighting_height 3000    //光源Buffer
#define lighting_width 3000     //光源Buffer
#define model_count 10

char address[model_count][50];
int texWidth, texHeight, texChannels;

///这块要改路径!这块要改路径!这块要改路径!这块要改路径!这块要改路径!这块要改路径!
unsigned char* textureData = stbi_load("D:\\Project\\计算机图形学原理\\场景漫游\\场景漫游\\WoodQuarteredOrangeGrove001_COL_3K.jpg", &texWidth, &texHeight, &texChannels, 0);
void init_address()
{
    strcat(address[0], "D://Project//计算机图形学原理//场景漫游//场景漫游//model//backwall.obj");
    strcat(address[1], "D://Project//计算机图形学原理//场景漫游//场景漫游//model//windows.obj");
    strcat(address[2], "D://Project//计算机图形学原理//场景漫游//场景漫游//model//chair1.obj");
    strcat(address[3], "D://Project//计算机图形学原理//场景漫游//场景漫游//model//chair2.obj");
    strcat(address[4], "D://Project//计算机图形学原理//场景漫游//场景漫游//model//desk.obj");
    strcat(address[5], "D://Project//计算机图形学原理//场景漫游//场景漫游//model//floor.obj");
    strcat(address[6], "D://Project//计算机图形学原理//场景漫游//场景漫游//model//fridge.obj");
    strcat(address[7], "D://Project//计算机图形学原理//场景漫游//场景漫游//model//leftwall.obj");
    strcat(address[8], "D://Project//计算机图形学原理//场景漫游//场景漫游//model//mirror.obj");
    strcat(address[9], "D://Project//计算机图形学原理//场景漫游//场景漫游//model//bowl.obj");
}

///
int nearplane_width = 1000;		//视景体宽度
int nearplane_height = 1000;		//视景体高度
int nearplane_distance = 800; //视景体近平面与视点距离
int farplane_distance = nearplane_distance + 200; //视景体远平面与视点距离
int shadow_mapping_temp = 0;
int UV_SWITCH = 0;     //开关UV

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
    color operator*(color uv)
    {
        color temp;
        temp.r = this->r * uv.r;
        temp.g = this->g * uv.g;
        temp.b = this->b * uv.b;
        return temp;
    }
};

class my_homogeneous_point
{
public:
    double x;
    double y;
    double z;
    double ratio;
public:
    my_homogeneous_point() {}

    my_homogeneous_point(double x, double y, double z, double ratio = 1)
    {
        this->x = x;
        this->y = y;
        this->z = z;

        this->ratio = ratio;
    }

    my_homogeneous_point normalize()
    {
        if (this->ratio != 0)
            return my_homogeneous_point(this->x / this->ratio, this->y / this->ratio, this->z, 1);
        return my_homogeneous_point(this->x, this->y, this->z, this->ratio);
    }

    my_homogeneous_point operator*(double n)
    {
        double new_dx = this->x * n;
        double new_dy = this->y * n;
        double new_dz = this->z * n;
        return my_homogeneous_point(new_dx, new_dy, new_dz, this->ratio);
    }

    my_homogeneous_point operator/(double n)
    {
        double new_dx = this->x / n;
        double new_dy = this->y / n;
        double new_dz = this->z / n;
        return my_homogeneous_point(new_dx, new_dy, new_dz, this->ratio);
    }

    my_homogeneous_point operator+(my_homogeneous_point input_vector)
    {
        double new_dx = this->x + input_vector.x;
        double new_dy = this->y + input_vector.y;
        double new_dz = this->z + input_vector.z;
        return my_homogeneous_point(new_dx, new_dy, new_dz, this->ratio);
    }

    my_homogeneous_point operator-(my_homogeneous_point input_vector)
    {
        double new_dx = this->x - input_vector.x;
        double new_dy = this->y - input_vector.y;
        double new_dz = this->z - input_vector.z;
        return my_homogeneous_point(new_dx, new_dy, new_dz, this->ratio);
    }

    my_homogeneous_point operator-()
    {
        double new_dx = -this->x;
        double new_dy = -this->y;
        double new_dz = -this->z;
        return my_homogeneous_point(new_dx, new_dy, new_dz, this->ratio);
    }

};

class my_3Dvector
{
public:
    double dx;
    double dy;
    double dz;
    double len;
public:
    my_3Dvector() {}

    my_3Dvector(double x, double y, double z)
    {
        dx = x;
        dy = y;
        dz = z;

        len = sqrtf(powf(dx, 2) + powf(dy, 2) + powf(dz, 2));
    }

    //start点指向end点的向量
    my_3Dvector(my_homogeneous_point start, my_homogeneous_point end)
    {
        dx = end.x - start.x;
        dy = end.y - start.y;
        dz = end.z - start.z;

        len = sqrtf(powf(dx, 2) + powf(dy, 2) + powf(dz, 2));
    }

    //叉乘 this X input_vector
    my_3Dvector cross_multiply(my_3Dvector input_vector)
    {
        double new_dx = dy * input_vector.dz - dz * input_vector.dy;
        double new_dy = dz * input_vector.dx - dx * input_vector.dz;
        double new_dz = dx * input_vector.dy - dy * input_vector.dx;
        return  my_3Dvector(new_dx, new_dy, new_dz);
    }

    //点乘 this * input_vector
    double dot_multiply(my_3Dvector input_vector)
    {
        return dx * input_vector.dx + dy * input_vector.dy + dz * input_vector.dz;
    }

    my_3Dvector normalize()
    {
        return my_3Dvector(this->dx / this->len, this->dy / this->len, this->dz / this->len);
    }

    my_3Dvector operator*(int n)
    {
        double new_dx = dx * n;
        double new_dy = dy * n;
        double new_dz = dz * n;
        return my_3Dvector(new_dx, new_dy, new_dz);
    }

    my_3Dvector operator-(my_3Dvector input_vector)
    {
        double new_dx = dx - input_vector.dx;
        double new_dy = dy - input_vector.dy;
        double new_dz = dz - input_vector.dz;
        return my_3Dvector(new_dx, new_dy, new_dz);
    }

    my_3Dvector operator-()
    {
        double new_dx = -dx;
        double new_dy = -dy;
        double new_dz = -dz;
        return my_3Dvector(new_dx, new_dy, new_dz);
    }
};

class NET
{
public:
    double x_min;
    double z_min;
    double Delta_x;
    double Delta_z;

    my_homogeneous_point vn_min;
    my_homogeneous_point Delta_vn;

    my_homogeneous_point p_min;
    my_homogeneous_point Delta_p;

    my_homogeneous_point vt_min;
    my_homogeneous_point Delta_vt;

    int flag;
    double y_max;
    NET* next_x;
    NET* next_y;
public:
    NET()
    {
        this->flag = 0;
        this->next_x = NULL;
        this->next_y = NULL;
    }

    NET(my_homogeneous_point a, my_homogeneous_point b)
    {
        this->Delta_x = (b.x - a.x) / (b.y - a.y);
        this->Delta_z = (b.z - a.z) / (b.y - a.y);
        this->x_min = a.x;
        this->z_min = a.z;

        this->flag = 0;

        this->y_max = b.y;
        this->next_x = NULL;
        this->next_y = NULL;
    }

    NET(my_homogeneous_point a, my_homogeneous_point b, my_homogeneous_point vn_a, my_homogeneous_point vn_b)
    {
        this->Delta_x = (b.x - a.x) / (b.y - a.y);
        this->Delta_z = (b.z - a.z) / (b.y - a.y);
        this->x_min = a.x;
        this->z_min = a.z;

        this->flag = 1;

        this->Delta_vn = (vn_b - vn_a) / (b.y - a.y);
        this->vn_min = vn_a;

        this->y_max = b.y;
        this->next_x = NULL;
        this->next_y = NULL;
    }

    NET(my_homogeneous_point a, my_homogeneous_point b, my_homogeneous_point vn_a, my_homogeneous_point vn_b, my_homogeneous_point p_a, my_homogeneous_point p_b)
    {
        this->Delta_x = (b.x - a.x) / (b.y - a.y);
        this->Delta_z = (b.z - a.z) / (b.y - a.y);
        this->x_min = a.x;
        this->z_min = a.z;

        this->flag = 2;

        this->Delta_vn = (vn_b - vn_a) / (b.y - a.y);
        this->vn_min = vn_a;

        this->Delta_p = (p_b - p_a) / (b.y - a.y);
        this->p_min = p_a;

        this->y_max = b.y;
        this->next_x = NULL;
        this->next_y = NULL;
    }

    NET(my_homogeneous_point a, my_homogeneous_point b, my_homogeneous_point vn_a, my_homogeneous_point vn_b, my_homogeneous_point p_a, my_homogeneous_point p_b, my_homogeneous_point vt_a, my_homogeneous_point vt_b)
    {
        this->Delta_x = (b.x - a.x) / (b.y - a.y);
        this->Delta_z = (b.z - a.z) / (b.y - a.y);
        this->x_min = a.x;
        this->z_min = a.z;

        this->flag = 3;

        this->Delta_vn = (vn_b - vn_a) / (b.y - a.y);
        this->vn_min = vn_a;

        this->Delta_p = (p_b - p_a) / (b.y - a.y);
        this->p_min = p_a;

        this->Delta_vt = (vt_b - vt_a) / (b.y - a.y);
        this->vt_min = vt_a;

        this->y_max = b.y;
        this->next_x = NULL;
        this->next_y = NULL;
    }
};

typedef NET* LinkNET;

class AET
{
public:
    double x;
    double z;
    double Delta_x;
    double Delta_z;

    my_homogeneous_point vn;
    my_homogeneous_point Delta_vn;

    my_homogeneous_point p;
    my_homogeneous_point Delta_p;

    my_homogeneous_point vt;
    my_homogeneous_point Delta_vt;

    int flag;
    double y_max;
    AET* next;
public:
    AET()
    {
        this->next = NULL;
    }

    AET(LinkNET net)
    {
        this->x = net->x_min;
        this->z = net->z_min;

        this->Delta_x = net->Delta_x;
        this->Delta_z = net->Delta_z;

        if (net->flag >= 1)
        {
            this->vn = net->vn_min;
            this->Delta_vn = net->Delta_vn;
        }

        if (net->flag >= 2)
        {
            this->p = net->p_min;
            this->Delta_p = net->Delta_p;
        }

        if (net->flag >= 3)
        {
            this->vt = net->vt_min;
            this->Delta_vt = net->Delta_vt;
        }

        this->flag = net->flag;

        this->y_max = net->y_max;
        this->next = NULL;
    }

    void add()
    {
        this->x += this->Delta_x;
        this->z += this->Delta_z;

        if (this->flag >= 1)
        {
            this->vn = this->vn + this->Delta_vn;
        }

        if (this->flag >= 2)
        {
            this->p = this->p + this->Delta_p;
        }

        if (this->flag >= 3)
        {
            this->vt = this->vt + this->Delta_vt;
        }
    }
};

typedef AET* LinkAET;

typedef struct {
    std::vector<my_homogeneous_point> point;
    std::vector<my_homogeneous_point> vt;
    std::vector<my_homogeneous_point> vn;

}Face;

typedef struct {
    std::vector<Face> faces;
    std::vector<my_homogeneous_point> v;
    std::vector<my_homogeneous_point> vt;
    std::vector<my_homogeneous_point> vn;
}Body;

typedef struct
{
    my_homogeneous_point pre_p;
    double Depth;
}DP;

Body Obj[model_count];
Body obj1[1];
Body obj2[1];
Body project_obj[1];
Body lighting_project_obj[1];

my_homogeneous_point clipping_polygon[clipping_polygon_point_count];
std::vector<my_homogeneous_point> clipped_face;

int clipping_start_point_x = -clipping_width / 2;
int clipping_start_point_y = -clipping_height / 2;

int lighting_start_point_x = -lighting_width / 2;
int lighting_start_point_y = -lighting_height / 2;

color rgb = { 0,0,0 };
DP Depth_Buffer[clipping_height][clipping_width];
color Color_Buffer[clipping_height][clipping_width];
double Lighting_Depth_Buffer[lighting_height][lighting_width];

///
#define K_AMBIENT 0.5   //环境光系数
#define K_DIFFUSE 0.3  //漫反射系数
#define K_SPECULAR 0.4  //镜面反射系数
#define ALPHA_SMOOTH 2  //光滑度

color lighting_b = { 0.5,0.5,0.5 };     //背景光
color lighting_p = { 1,1,1 };     //点光
my_homogeneous_point eye(0, 0, nearplane_distance);    //视点
my_homogeneous_point point_lighting(200, 200, 800);  //点光位置

my_homogeneous_point camera[3];
my_homogeneous_point camera_position;

///
color Get_texture_color(double u, double v, unsigned char* textureData, int texWidth, int texChannels)
{
    int index = (v * texWidth + u) * texChannels;
    color a;
    a.r = textureData[index] / 255.0f;
    a.g = textureData[index + 1] / 255.0f;
    a.b = textureData[index + 2] / 255.0f;
    return a;
}

color UV(my_homogeneous_point p)
{
    color a(0, 0, 0);

    //p.x = p.x * 0.5 + 0.5;
    //p.y = p.y * 0.5 + 0.5;

    int texU = int(p.x * texWidth);
    int texV = int(p.y * texHeight);

    texU = (texU + texWidth) / 2 - 1;
    texV = (texV + texHeight) / 2 - 1;

    a = Get_texture_color(texU, texV, textureData, texWidth, texChannels);

    return a;
}
///
void setPixel2D(int x, int y)
{
    glBegin(GL_POINTS);
    glVertex2i(x, y);
    glEnd();
}

void DDALine(int x0, int y0, int x1, int y1)
{
    double dy = y1 - y0;
    double dx = x1 - x0;
    double xf = x0;
    double yf = y0;
    double k = dy / dx;
    if (k <= 1 && k >= -1)
    {
        if (xf > x1)
        {
            while (xf > x1)
            {
                //对当前点进行舍入，获取像素点坐标
                int xr = (int)(xf - 0.5);
                int yr = (int)(yf - 0.5);
                setPixel2D(xr, yr);

                //计算直线段上下一个点
                xf = xf - 1;
                yf = yf - k;
            }
        }
        else
        {
            while (xf < x1)
            {
                //对当前点进行舍入，获取像素点坐标
                int xr = (int)(xf + 0.5);
                int yr = (int)(yf + 0.5);
                setPixel2D(xr, yr);

                //计算直线段上下一个点
                xf = xf + 1;
                yf = yf + k;
            }
        }
    }
    else
    {
        k = 1 / k;
        if (yf > y1)
        {
            while (yf > y1)
            {
                //对当前点进行舍入，获取像素点坐标
                int xr = (floor)(xf - 0.5);
                int yr = (floor)(yf - 0.5);
                setPixel2D(xr, yr);

                //计算直线段上下一个点
                xf = xf - k;
                yf = yf - 1;
            }
        }
        else
        {
            while (yf < y1)
            {
                //对当前点进行舍入，获取像素点坐标
                int xr = (floor)(xf + 0.5);
                int yr = (floor)(yf + 0.5);
                setPixel2D(xr, yr);

                //计算直线段上下一个点
                xf = xf + k;
                yf = yf + 1;
            }
        }
    }
}

void DDA_Common(my_homogeneous_point* points, int point_count)
{
    for (int vindex = 0; vindex < point_count; vindex++)
    {
        my_homogeneous_point p1 = points[vindex];
        my_homogeneous_point p2 = points[(vindex + 1) % point_count];
        glBegin(GL_POINTS);
        DDALine(p1.x, p1.y, p2.x, p2.y);
        glEnd();
    }
}

void load_obj(const char* filename, Body* obj)
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
            my_homogeneous_point temp;
            iss >> temp.x >> temp.y >> temp.z;
            temp.ratio = 1;
            obj->v.push_back(temp);
        }
        else if (type == "vn")
        {
            my_homogeneous_point temp;
            iss >> temp.x >> temp.y >> temp.z;
            obj->vn.push_back(temp);
        }
        else if (type == "vt")
        {
            my_homogeneous_point temp;
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
                iss >> temp[0] >> slash >> temp[1] >> slash >> temp[2];
                Temp.point.push_back(obj->v[temp[0] - 1]);
                Temp.vt.push_back(obj->vt[temp[1] - 1]);
                Temp.vn.push_back(obj->vn[temp[2] - 1]);
            }
            obj->faces.push_back(Temp);
        }
    }
}

void init()
{
    init_address();

    for (int i = 0; i < model_count; i++)
    {
        load_obj(address[i], Obj + i);
        for (int j = 0; j < Obj[i].faces.size(); j++)
        {
            obj1->faces.push_back(Obj[i].faces[j]);
        }
    }

    for (int i = 0; i < clipping_height; i++)
    {
        for (int j = 0; j < clipping_width; j++)
        {
            Depth_Buffer[i][j].Depth = -INT_MAX;
            Color_Buffer[i][j] = 1;
        }
    }

    for (int i = 0; i < lighting_height; i++)
    {
        for (int j = 0; j < lighting_width; j++)
        {
            Lighting_Depth_Buffer[i][j] = -INT_MAX;
        }
    }

    clipping_polygon[0].x = clipping_start_point_x;
    clipping_polygon[0].y = clipping_start_point_y;
    clipping_polygon[0].z = 0;
    clipping_polygon[0].ratio = 1;

    clipping_polygon[1].x = -clipping_start_point_x;
    clipping_polygon[1].y = clipping_start_point_y;
    clipping_polygon[1].z = 0;
    clipping_polygon[1].ratio = 1;

    clipping_polygon[2].x = -clipping_start_point_x;
    clipping_polygon[2].y = -clipping_start_point_y;
    clipping_polygon[2].z = 0;
    clipping_polygon[2].ratio = 1;

    clipping_polygon[3].x = clipping_start_point_x;
    clipping_polygon[3].y = -clipping_start_point_y;
    clipping_polygon[3].z = 0;
    clipping_polygon[3].ratio = 1;
}

///
my_homogeneous_point matrix_multiply_vector(double matrix[][4], my_homogeneous_point input_v)
{
    my_homogeneous_point translated_v;
    translated_v.x = matrix[0][0] * input_v.x + matrix[0][1] * input_v.y + matrix[0][2] * input_v.z + matrix[0][3] * input_v.ratio;
    translated_v.y = matrix[1][0] * input_v.x + matrix[1][1] * input_v.y + matrix[1][2] * input_v.z + matrix[1][3] * input_v.ratio;
    translated_v.z = matrix[2][0] * input_v.x + matrix[2][1] * input_v.y + matrix[2][2] * input_v.z + matrix[2][3] * input_v.ratio;
    translated_v.ratio = matrix[3][0] * input_v.x + matrix[3][1] * input_v.y + matrix[3][2] * input_v.z + matrix[3][3] * input_v.ratio;
    return translated_v;
}

void my_traslate_homogeneous(Body* object, double tx, double ty, double tz)
{
    double translate_matrix[4][4];
    memset(translate_matrix, 0, sizeof(double) * 16);
    translate_matrix[0][0] = translate_matrix[1][1] = translate_matrix[2][2] = translate_matrix[3][3] = 1;//对角线赋值为1
    translate_matrix[0][3] = tx;
    translate_matrix[1][3] = ty;
    translate_matrix[2][3] = tz;

    for (int i = 0; i < object->faces.size(); i++)
    {
        for (int vIndex = 0; vIndex < 3; vIndex++)
        {
            my_homogeneous_point input_v = object->faces[i].point[vIndex];
            input_v = matrix_multiply_vector(translate_matrix, input_v); //平移矩阵作用到每个顶点
            object->faces[i].point[vIndex] = input_v;
        }
    }
}

void my_rotateY_homogeneous(Body* object, int deg) {
    double rad = deg * PI / 180;
    double rMat[4][4] = {
        {cos(rad), 0, sin(rad), 0},
        {0, 1, 0, 0},
        {-sin(rad), 0, cos(rad), 0},
        {0, 0, 0, 1}
    };

    for (int i = 0; i < object->faces.size(); i++)
    {
        for (int vIndex = 0; vIndex < 3; vIndex++)
        {
            my_homogeneous_point input_v = object->faces[i].point[vIndex];
            input_v = matrix_multiply_vector(rMat, input_v); //平移矩阵作用到每个顶点
            object->faces[i].point[vIndex] = input_v;
            input_v = object->faces[i].vn[vIndex];
            input_v = matrix_multiply_vector(rMat, input_v);
            object->faces[i].vn[vIndex] = input_v;
        }
    }
}

void projected(Body* obj, Body* result_obj, int nearplane)
{
    double P[4][4];
    memset(P, 0, sizeof(double) * 16);
    P[0][0] = 1;
    P[1][1] = 1;
    P[3][3] = 1;
    P[3][2] = -1.0 / nearplane;
    //P[3][2] = 0;

    result_obj->faces.clear();
    for (int i = 0; i < obj->faces.size(); i++)
    {
        Face Temp = obj->faces[i];
        for (int vIndex = 0; vIndex < 3; vIndex++)
        {
            my_homogeneous_point input_v = obj->faces[i].point[vIndex];
            double temp = input_v.z;
            input_v = matrix_multiply_vector(P, input_v);
            input_v = input_v.normalize();
            Temp.point[vIndex] = input_v;
            Temp.point[vIndex].z = temp;
        }
        result_obj->faces.push_back(Temp);
    }
};
///
bool determine_point_edge_position(my_homogeneous_point edge_start_point, my_homogeneous_point edge_end_point, my_homogeneous_point given_point)
{
    my_3Dvector face_normal(0, 0, 1);
    my_3Dvector startLine(edge_start_point, edge_end_point);
    my_3Dvector endLine(edge_start_point, given_point);
    my_3Dvector cross_vector = startLine.cross_multiply(endLine);
    return cross_vector.dot_multiply(face_normal) >= 0 ? false : true;
}

my_homogeneous_point compute_intersection(my_homogeneous_point clip_start,
    my_homogeneous_point clip_end,
    my_homogeneous_point cur_point,
    my_homogeneous_point next_point)
{
    my_homogeneous_point intersection;

    // 计算裁剪边和被裁剪边的向量
    float dx1 = clip_end.x - clip_start.x;
    float dy1 = clip_end.y - clip_start.y;

    float dx2 = next_point.x - cur_point.x;
    float dy2 = next_point.y - cur_point.y;

    // 计算向量叉积
    float denominator = dx1 * dy2 - dy1 * dx2;

    //防止平行情况（即没有交点）
    if (fabs(denominator) < 1e-5) // 用一个较小的值作为误差容忍
    {
        intersection.x = cur_point.x; // 返回一个默认点（也可以抛出错误或其他处理方式）
        intersection.y = cur_point.y;
        intersection.z = cur_point.z;
        intersection.ratio = 1;
        return intersection;
    }

    // 计算交点的参数 t 和 u
    float t = ((cur_point.x - clip_start.x) * dy2 - (cur_point.y - clip_start.y) * dx2) / denominator;

    // 计算交点位置
    intersection.x = clip_start.x + t * dx1;
    intersection.y = clip_start.y + t * dy1;
    // 使用 t 来计算 z 坐标
    intersection.z = cur_point.z + t * (next_point.z - cur_point.z); // 直接从 cur_point 到 next_point 插值
    //intersection.z = 0;
    intersection.ratio = 1;

    return intersection;
}

void Sutherland_Hodgman_Clipping(const my_homogeneous_point* clipping_polygon, std::vector<my_homogeneous_point> subject_polygon, int subject_count)
{
    int clip_count = 4;
    clipped_face.clear();
    std::vector<my_homogeneous_point> input_polygon = subject_polygon;
    std::vector<my_homogeneous_point> output_polygon;

    for (int i = 0; i < clip_count; i++)
    {
        my_homogeneous_point clip_start = clipping_polygon[i];
        my_homogeneous_point clip_end = clipping_polygon[(i + 1) % clip_count];

        output_polygon.clear();

        for (size_t j = 0; j < input_polygon.size(); j++)
        {
            my_homogeneous_point cur_point = input_polygon[j];
            my_homogeneous_point next_point = input_polygon[(j + 1) % input_polygon.size()];

            bool cur_position = !determine_point_edge_position(clip_start, clip_end, cur_point);
            bool next_position = !determine_point_edge_position(clip_start, clip_end, next_point);

            if (cur_position && next_position)
            {
                output_polygon.push_back(next_point);
            }
            else if (cur_position && !next_position)
            {
                my_homogeneous_point intersection = compute_intersection(clip_start, clip_end, cur_point, next_point);
                output_polygon.push_back(intersection);
            }
            else if (!cur_position && next_position)
            {
                my_homogeneous_point intersection = compute_intersection(clip_start, clip_end, cur_point, next_point);
                output_polygon.push_back(intersection);
                output_polygon.push_back(next_point);
            }
        }

        // 更新 input_polygon 为当前裁剪结果
        input_polygon = output_polygon;
    }

    // 将最终结果输出到 clipped_polygon
    clipped_face.insert(clipped_face.begin(), output_polygon.begin(), output_polygon.end());
}
///
bool judge_face(std::vector<my_homogeneous_point> points, my_3Dvector now_view)
{
    if (points.size() < 3)
    {
        return true;
    }
    my_3Dvector v1(points[0], points[1]);
    my_3Dvector v2(points[1], points[2]);
    my_3Dvector cross = v1.cross_multiply(v2);
    cross = cross.normalize();
    return cross.dot_multiply(now_view) < -1e-4 ? true : false;
}

void Linear_interpolation(LinkAET p, LinkAET k, int y, color face_color)
{
    AET temp1 = p->x < k->x ? *p : *k;
    AET temp2 = p->x < k->x ? *k : *p;

    //temp1.x = (int)temp1.x;
    //temp2.x = (int)temp2.x;
    temp1.x = (int)(temp1.x + 0.5 * (temp1.x >= 0 ? 1 : -1));
    temp2.x = (int)(temp2.x + 0.5 * (temp2.x >= 0 ? 1 : -1));

    double z_per_x = (temp2.z - temp1.z) / (temp2.x - temp1.x);
    double z_now = temp1.z;

    for (int i = temp1.x; i <= temp2.x; i++)
    {
        //确保访问的索引在 Depth_Buffer 和 Color_Buffer 的范围内
        int buffer_x = i - clipping_start_point_x;
        int buffer_y = y - clipping_start_point_y;
        if (buffer_x >= 0 && buffer_x < clipping_width && buffer_y >= 0 && buffer_y < clipping_height)
        {
            if (z_now > Depth_Buffer[buffer_x][buffer_y].Depth || Depth_Buffer[buffer_x][buffer_y].Depth == -INT_MAX)
            {
                Depth_Buffer[buffer_x][buffer_y].Depth = z_now;
                Color_Buffer[buffer_x][buffer_y] = face_color;
            }
        }
        z_now += z_per_x;
    }
}

void Buffer_fill(std::vector<my_homogeneous_point> polygon, color face_color)
{
    int polygon_count = polygon.size();
    for (int i = 0; i < polygon_count; i++)
    {
        polygon[i].y = (int)(polygon[i].y + 0.5 * (polygon[i].y < 0 ? -1 : 1));
    }

    // 确定多边形的y范围
    int max = polygon[0].y;
    int min = polygon[0].y;
    for (int i = 1; i < polygon_count; i++)
    {
        max = fmax(max, polygon[i].y);
        min = fmin(min, polygon[i].y);
    }

    //构造边表
    LinkNET Net = new NET;
    LinkNET Net_clean = Net;
    LinkNET p = Net;

    for (int i = min; i <= max; i++)
    {
        LinkNET q = new NET;
        p->next_y = q;
        p = p->next_y;
        for (int j = 0; j < polygon_count; j++)
        {
            if (polygon[j].y == i)
            {
                // 判断该点相邻的两点与该点的位置关系
                // 如果相邻两点均高于该店，则插入两个边的信息
                if (polygon[(j - 1 + polygon_count) % polygon_count].y > polygon[j].y && polygon[(j + 1) % polygon_count].y > polygon[j].y)
                {
                    LinkNET Temp1 = new NET(polygon[j], polygon[(j + 1) % polygon_count]);
                    LinkNET Temp2 = new NET(polygon[j], polygon[(j - 1 + polygon_count) % polygon_count]);

                    if (Temp1->Delta_x < Temp2->Delta_x)
                    {
                        Temp1->next_x = Temp2;
                        Temp2->next_x = NULL;
                        q->next_x = Temp1;
                    }
                    else
                    {
                        Temp2->next_x = Temp1;
                        Temp1->next_x = NULL;
                        q->next_x = Temp2;
                    }
                    q = q->next_x->next_x;
                }

                // 如果相邻两点只有一点高于该店，则插入一个边的信息
                else if (polygon[(j + 1) % polygon_count].y > polygon[j].y)
                {
                    LinkNET Temp = new NET(polygon[j], polygon[(j + 1) % polygon_count]);
                    q->next_x = Temp;
                    q = q->next_x;
                }

                // 如果相邻两点只有一点高于该店，则插入一个边的信息
                else if (polygon[(j - 1 + polygon_count) % polygon_count].y > polygon[j].y)
                {
                    LinkNET Temp = new NET(polygon[j], polygon[(j - 1 + polygon_count) % polygon_count]);
                    q->next_x = Temp;
                    q = q->next_x;
                }
            }
        }
    }

    // 初始化活性边表
    LinkAET Aet = new AET;
    LinkAET Aet_clean = Aet;

    // 填充多边形
    while (Net->next_y != NULL && Net->next_y->next_y != NULL) // 确保当前的Net节点存在，确保下一行的Net节点存在。
    {
        LinkNET q = Net->next_y->next_x;

        // 插入排序遍历AET并按照x坐标值递增的顺序插入Temp。如果x值相同，比较Delta_x来决定插入位置，以保持排序稳定性。
        while (q != NULL)
        {
            LinkAET p = Aet;
            LinkAET Temp = new AET(q);
            while (p->next != NULL)
            {
                if (p->next->x > Temp->x)
                {
                    Temp->next = p->next;
                    break;
                }
                else if (p->next->x == Temp->x)
                {
                    if (p->next->Delta_x < Temp->Delta_x)
                    {
                        Temp->next = p->next->next;
                    }
                    else
                    {
                        Temp->next = p->next;
                        break;
                    }
                }
                p = p->next;
            }
            p->next = Temp;
            q = q->next_x;
        }


        // 更新AET和绘制多边形
        LinkAET temp = Aet;

        while (temp->next != NULL && temp->next->next != NULL)
        {
            int flag = 1;
            LinkAET p = temp->next;
            LinkAET k = p->next;
            LinkAET TEMP1 = NULL;
            LinkAET TEMP2 = NULL;

            Linear_interpolation(p, k, min, face_color);

            p->add();
            k->add();

            min++;
            Net = Net->next_y;

            if (min == p->y_max && min == k->y_max)
            {
                TEMP1 = p;
                TEMP2 = k;
                temp->next = k->next;
                flag = 0;
            }
            else
            {
                if (min == k->y_max)
                {
                    TEMP2 = k;
                    p->next = k->next;
                    flag = 0;
                }
                if (min == p->y_max)
                {
                    TEMP1 = p;
                    temp->next = k;
                    flag = 0;
                }
            }
            if (TEMP1 != NULL)
            {
                delete TEMP1;
            }
            if (TEMP2 != NULL)
            {
                delete TEMP2;
            }
            if (flag)
                temp = temp->next->next;
        }
    }

    while (Aet_clean->next != NULL)
    {
        LinkAET TEMP = Aet_clean;
        Aet_clean = Aet_clean->next;
        delete TEMP;
    }
    delete Aet_clean;

    while (Net_clean->next_y != NULL)
    {
        LinkNET TEMP_X = Net_clean->next_x;
        LinkNET TEMP_Y = Net_clean;
        Net_clean = Net_clean->next_y;
        delete TEMP_Y;

        if (TEMP_X != NULL)
        {
            while (TEMP_X->next_x != NULL)
            {
                LinkNET TEMP_X_TEMP = TEMP_X;
                TEMP_X = TEMP_X->next_x;
                delete TEMP_X_TEMP;
            }
            delete TEMP_X;
        }
    }
    delete Net_clean;
    return;
}

void Show()
{
    for (int i = 0; i < clipping_height; i++)
    {
        for (int j = 0; j < clipping_width; j++)
        {
            double r = Color_Buffer[i][j].r;
            double g = Color_Buffer[i][j].g;
            double b = Color_Buffer[i][j].b;

            glBegin(GL_POINTS);

            glColor3f(r, g, b);
            glVertex2i(clipping_start_point_x + i, clipping_start_point_y + j);

            glEnd();
        }
    }

    for (int i = 0; i < clipping_height; i++)
    {
        for (int j = 0; j < clipping_width; j++)
        {
            Depth_Buffer[i][j].Depth = -INT_MAX;
            Color_Buffer[i][j] = 1;
        }
    }

    for (int i = 0; i < lighting_height; i++)
    {
        for (int j = 0; j < lighting_width; j++)
        {
            Lighting_Depth_Buffer[i][j] = -INT_MAX;
        }
    }
}
///
void init_camera(my_homogeneous_point start, my_homogeneous_point end, Body* obj)
{
    my_3Dvector temp1(end, { 0,0,0 });
    temp1 = temp1.normalize();
    my_3Dvector temp2(0, -1, 0);
    my_3Dvector temp3 = temp2.cross_multiply(temp1);
    temp3 = temp3.normalize();
    temp2 = temp1.cross_multiply(temp3);
    temp2 = temp2.normalize();

    camera[0].x = temp1.dx;
    camera[0].y = temp1.dy;
    camera[0].z = temp1.dz;
    camera[0].ratio = 1;//look at

    camera[1].x = temp2.dx;
    camera[1].y = temp2.dy;
    camera[1].z = temp2.dz;
    camera[1].ratio = 1;//updirect

    camera[2].x = temp3.dx;
    camera[2].y = temp3.dy;
    camera[2].z = temp3.dz;
    camera[2].ratio = 1;// look at X up direct

    //camera_position.x = end.x;
    //camera_position.y = end.y;
    //camera_position.z = end.z;

    my_3Dvector temp(start, end);

    camera_position.x = temp.dx;
    camera_position.y = temp.dy;
    camera_position.z = temp.dz;//视点位置

}

void init_camera()
{
    camera[0].x = 0;
    camera[0].y = 0;
    camera[0].z = -1;
    camera[0].ratio = 1;//look at

    camera[1].x = 0;
    camera[1].y = 1;
    camera[1].z = 0;
    camera[1].ratio = 1;//updirect

    camera[2].x = 1;
    camera[2].y = 0;
    camera[2].z = 0;
    camera[2].ratio = 1;// look at X up direct

    camera_position.x = 0;
    camera_position.y = 0;
    camera_position.z = nearplane_distance;
    camera_position.ratio = 1;//视点位置
}

my_homogeneous_point camera_change(my_homogeneous_point point)
{
    //装配生成视点变换矩阵 
    double camera_rotate_matrix[4][4];
    double camera_tanslate_matrix[4][4];
    memset(camera_rotate_matrix, 0, sizeof(double) * 16);
    memset(camera_tanslate_matrix, 0, sizeof(double) * 16);

    ////MR-1
    camera_rotate_matrix[0][0] = -camera[2].x;
    camera_rotate_matrix[0][1] = -camera[2].y;
    camera_rotate_matrix[0][2] = camera[2].z;

    camera_rotate_matrix[1][0] = -camera[1].x;
    camera_rotate_matrix[1][1] = -camera[1].y;
    camera_rotate_matrix[1][2] = camera[1].z;

    camera_rotate_matrix[2][0] = camera[0].x;
    camera_rotate_matrix[2][1] = camera[0].y;
    camera_rotate_matrix[2][2] = -camera[0].z;

    camera_rotate_matrix[3][3] = 1;
    ////MT-1
    camera_tanslate_matrix[0][0] = 1;
    camera_tanslate_matrix[1][1] = 1;
    camera_tanslate_matrix[2][2] = 1;
    camera_tanslate_matrix[0][3] = -camera_position.x;
    camera_tanslate_matrix[1][3] = -camera_position.y;
    camera_tanslate_matrix[2][3] = -camera_position.z;
    camera_tanslate_matrix[3][3] = 1;

    //进行平行矩阵MT-1

    my_homogeneous_point result;
    my_homogeneous_point input_v;
    input_v = point;
    input_v = input_v.normalize();
    input_v = matrix_multiply_vector(camera_tanslate_matrix, input_v);
    input_v = input_v.normalize();
    input_v = matrix_multiply_vector(camera_rotate_matrix, input_v);
    result = input_v.normalize();
    return result;
}

void obj_change(Body* obj1, Body* obj2)
{
    init_camera();
    init_camera(eye, point_lighting, obj1);
    for (int i = 0; i < obj1->faces.size(); i++)
    {
        Face Temp;
        for (int j = 0; j < obj1->faces[i].point.size(); j++)
        {
            my_homogeneous_point temp = camera_change(obj1->faces[i].point[j]);
            Temp.point.push_back(temp);
        }
        obj2->faces.push_back(Temp);
    }
}
///
color ambient_lighting()
{
    color result;
    result = lighting_b * K_AMBIENT;
    return result;
}

color diffuse_lighting(my_3Dvector vector_normal, my_3Dvector vector_light)
{
    color result;
    double diff = fmax(0, vector_light.dot_multiply(vector_normal));
    result = lighting_p * diff * K_DIFFUSE;
    return result;
}

color specular_lighting(my_3Dvector vector_reflect, my_3Dvector vector_eye)
{
    color result;
    //double spec = vector_reflect.dot_multiply(vector_eye);
    double spec = fmax(0, vector_reflect.dot_multiply(vector_eye));
    spec = powf(spec, ALPHA_SMOOTH);
    result = lighting_p * K_SPECULAR * spec;
    return result;
}

color Phong_lighting(std::vector<my_homogeneous_point>polygon)
{
    my_homogeneous_point point_mid = { 0,0,0,1 };   //中点代替面
    for (int i = 0; i < polygon.size(); i++)
    {
        point_mid = point_mid + polygon[i] / polygon.size();
    }

    my_3Dvector vector_light(point_mid, point_lighting);    //入射光向量
    vector_light = vector_light.normalize();

    my_3Dvector v1(polygon[0], polygon[1]);
    my_3Dvector v2(polygon[1], polygon[2]);
    my_3Dvector vector_normal = v1.cross_multiply(v2);      //法向量
    vector_normal = vector_normal.normalize();

    my_3Dvector vector_eye(point_mid, eye);     //视线向量
    vector_eye = vector_eye.normalize();

    my_3Dvector vector_reflect = -vector_light - (vector_light.cross_multiply(vector_normal)).cross_multiply(vector_normal) * 2;
    vector_reflect = vector_reflect.normalize();        //反射光向量
    vector_reflect = -vector_reflect;

    color l_a, l_d, l_s;

    //l_a = { 0,0,0 };      //测试用
    //l_d = { 0,0,0 };
    //l_s = { 0,0,0 };

    l_a = ambient_lighting();
    l_d = diffuse_lighting(vector_normal, vector_light);
    l_s = specular_lighting(vector_reflect, vector_eye);

    color result;
    result = l_a + l_d + l_s;

    return result;
}
///
color Blinn_Phong_lighting(my_homogeneous_point point, my_homogeneous_point vn, color uv_color)
{
    my_3Dvector vector_light(point, point_lighting);    //入射光向量
    vector_light = vector_light.normalize();

    my_3Dvector vector_normal(vn.x, vn.y, vn.z);

    my_3Dvector vector_eye(point, eye);     //视线向量
    vector_eye = vector_eye.normalize();

    my_3Dvector vector_reflect = -vector_light - (vector_light.cross_multiply(vector_normal)).cross_multiply(vector_normal) * 2;
    vector_reflect = vector_reflect.normalize();        //反射光向量
    vector_reflect = -vector_reflect;

    color l_a, l_d, l_s;

    //l_a = { 0,0,0 };      //测试用
    //l_d = { 0,0,0 };
    //l_s = { 0,0,0 };

    l_a = ambient_lighting();
    l_d = diffuse_lighting(vector_normal, vector_light);
    l_s = specular_lighting(vector_reflect, vector_eye);

    color result;

    if (UV_SWITCH % 2 == 1)
    {
        result = (l_a + l_d) * uv_color + l_s;
    }
    else
    {
        result = l_a + l_d + l_s;
    }

    return result;
}

void Linear_interpolation_vnvt(LinkAET p, LinkAET k, int y)
{
    AET temp1 = p->x < k->x ? *p : *k;
    AET temp2 = p->x < k->x ? *k : *p;

    //temp1.x = (int)temp1.x;
    //temp2.x = (int)temp2.x;
    temp1.x = (int)(temp1.x + 0.5 * (temp1.x >= 0 ? 1 : -1));
    temp2.x = (int)(temp2.x + 0.5 * (temp2.x >= 0 ? 1 : -1));

    double z_per_x = (temp2.z - temp1.z) / (temp2.x - temp1.x);
    double z_now = temp1.z;

    my_homogeneous_point vn_per_x = (temp2.vn - temp1.vn) / (temp2.x - temp1.x);
    my_homogeneous_point vn_now = temp1.vn;

    my_homogeneous_point p_per_x = (temp2.p - temp1.p) / (temp2.x - temp1.x);
    my_homogeneous_point p_now = temp1.p;

    my_homogeneous_point vt_per_x = (temp2.vt - temp1.vt) / (temp2.x - temp1.x);
    my_homogeneous_point vt_now = temp1.vt;

    for (int i = temp1.x; i <= temp2.x; i++)
    {
        //确保访问的索引在 Depth_Buffer 和 Color_Buffer 的范围内
        int buffer_x = i - clipping_start_point_x;
        int buffer_y = y - clipping_start_point_y;
        if (buffer_x >= 0 && buffer_x < clipping_width && buffer_y >= 0 && buffer_y < clipping_height)
        {
            if (z_now > Depth_Buffer[buffer_x][buffer_y].Depth || Depth_Buffer[buffer_x][buffer_y].Depth == -INT_MAX)
            {
                my_homogeneous_point now_point(i, y, z_now);
                color uv_color = UV(vt_now);
                color light_color = Blinn_Phong_lighting(now_point, vn_now, uv_color);
                Depth_Buffer[buffer_x][buffer_y].Depth = z_now;
                Depth_Buffer[buffer_x][buffer_y].pre_p = p_now;
                Color_Buffer[buffer_x][buffer_y] = light_color;
            }
        }
        z_now += z_per_x;
        vn_now = vn_now + vn_per_x;
        p_now = p_now + p_per_x;
    }
}

void Edge_Table_Blinn_Phong(Face face, Face pre_face)
{
    std::vector<my_homogeneous_point> polygon = face.point;
    std::vector<my_homogeneous_point> vn = face.vn;
    std::vector<my_homogeneous_point> vt = face.vt;
    std::vector<my_homogeneous_point> pre_p = pre_face.point;

    int polygon_count = polygon.size();
    for (int i = 0; i < polygon_count; i++)
    {
        polygon[i].y = (int)(polygon[i].y + 0.5 * (polygon[i].y < 0 ? -1 : 1));
    }

    // 确定多边形的y范围
    int max = polygon[0].y;
    int min = polygon[0].y;
    for (int i = 1; i < polygon_count; i++)
    {
        max = fmax(max, polygon[i].y);
        min = fmin(min, polygon[i].y);
    }

    //构造边表
    LinkNET Net = new NET;
    LinkNET Net_clean = Net;
    LinkNET p = Net;

    for (int i = min; i <= max; i++)
    {
        LinkNET q = new NET;
        p->next_y = q;
        p = p->next_y;
        for (int j = 0; j < polygon_count; j++)
        {
            if (polygon[j].y == i)
            {
                // 判断该点相邻的两点与该点的位置关系
                // 如果相邻两点均高于该店，则插入两个边的信息
                if (polygon[(j - 1 + polygon_count) % polygon_count].y > polygon[j].y && polygon[(j + 1) % polygon_count].y > polygon[j].y)
                {
                    LinkNET Temp1 = new NET(polygon[j], polygon[(j + 1) % polygon_count],
                        vn[j], vn[(j + 1) % polygon_count],
                        pre_p[j], pre_p[(j + 1) % polygon_count],
                        vt[j], vt[(j + 1) % polygon_count]);

                    LinkNET Temp2 = new NET(polygon[j], polygon[(j - 1 + polygon_count) % polygon_count],
                        vn[j], vn[(j - 1 + polygon_count) % polygon_count],
                        pre_p[j], pre_p[(j - 1 + polygon_count) % polygon_count],
                        vt[j], vt[(j - 1 + polygon_count) % polygon_count]);

                    if (Temp1->Delta_x < Temp2->Delta_x)
                    {
                        Temp1->next_x = Temp2;
                        Temp2->next_x = NULL;
                        q->next_x = Temp1;
                    }
                    else
                    {
                        Temp2->next_x = Temp1;
                        Temp1->next_x = NULL;
                        q->next_x = Temp2;
                    }
                    q = q->next_x->next_x;
                }

                // 如果相邻两点只有一点高于该店，则插入一个边的信息
                else if (polygon[(j + 1) % polygon_count].y > polygon[j].y)
                {
                    LinkNET Temp = new NET(polygon[j], polygon[(j + 1) % polygon_count],
                        vn[j], vn[(j + 1) % polygon_count],
                        pre_p[j], pre_p[(j + 1) % polygon_count],
                        vt[j], vt[(j + 1) % polygon_count]);
                    q->next_x = Temp;
                    q = q->next_x;
                }

                // 如果相邻两点只有一点高于该店，则插入一个边的信息
                else if (polygon[(j - 1 + polygon_count) % polygon_count].y > polygon[j].y)
                {
                    LinkNET Temp = new NET(polygon[j], polygon[(j - 1 + polygon_count) % polygon_count],
                        vn[j], vn[(j - 1 + polygon_count) % polygon_count],
                        pre_p[j], pre_p[(j - 1 + polygon_count) % polygon_count],
                        vt[j], vt[(j - 1 + polygon_count) % polygon_count]);
                    q->next_x = Temp;
                    q = q->next_x;
                }
            }
        }
    }

    // 初始化活性边表
    LinkAET Aet = new AET;
    LinkAET Aet_clean = Aet;

    // 填充多边形
    while (Net->next_y != NULL && Net->next_y->next_y != NULL) // 确保当前的Net节点存在，确保下一行的Net节点存在。
    {
        LinkNET q = Net->next_y->next_x;

        // 插入排序遍历AET并按照x坐标值递增的顺序插入Temp。如果x值相同，比较Delta_x来决定插入位置，以保持排序稳定性。
        while (q != NULL)
        {
            LinkAET p = Aet;
            LinkAET Temp = new AET(q);

            while (p->next != NULL)
            {
                if (p->next->x > Temp->x)
                {
                    Temp->next = p->next;
                    break;
                }
                else if (p->next->x == Temp->x)
                {
                    if (p->next->Delta_x < Temp->Delta_x)
                    {
                        Temp->next = p->next->next;
                    }
                    else
                    {
                        Temp->next = p->next;
                        break;
                    }
                }
                p = p->next;
            }
            p->next = Temp;
            q = q->next_x;
        }


        // 更新AET和绘制多边形
        LinkAET temp = Aet;

        while (temp->next != NULL && temp->next->next != NULL)
        {
            int flag = 1;
            LinkAET p = temp->next;
            LinkAET k = p->next;
            LinkAET TEMP1 = NULL;
            LinkAET TEMP2 = NULL;

            Linear_interpolation_vnvt(p, k, min);

            p->add();
            k->add();

            min++;
            Net = Net->next_y;

            if (min == p->y_max && min == k->y_max)
            {
                TEMP1 = p;
                TEMP2 = k;
                temp->next = k->next;
                flag = 0;
            }
            else
            {
                if (min == k->y_max)
                {
                    TEMP2 = k;
                    p->next = k->next;
                    flag = 0;
                }
                if (min == p->y_max)
                {
                    TEMP1 = p;
                    temp->next = k;
                    flag = 0;
                }
            }
            if (TEMP1 != NULL)
            {
                delete TEMP1;
            }
            if (TEMP2 != NULL)
            {
                delete TEMP2;
            }
            if (flag)
                temp = temp->next->next;
        }
    }

    while (Aet_clean->next != NULL)
    {
        LinkAET TEMP = Aet_clean;
        Aet_clean = Aet_clean->next;
        delete TEMP;
    }
    delete Aet_clean;

    while (Net_clean->next_y != NULL)
    {
        LinkNET TEMP_X = Net_clean->next_x;
        LinkNET TEMP_Y = Net_clean;
        Net_clean = Net_clean->next_y;
        delete TEMP_Y;

        if (TEMP_X != NULL)
        {
            while (TEMP_X->next_x != NULL)
            {
                LinkNET TEMP_X_TEMP = TEMP_X;
                TEMP_X = TEMP_X->next_x;
                delete TEMP_X_TEMP;
            }
            delete TEMP_X;
        }
    }
    delete Net_clean;
    return;
}

void Only_Lighting_Linear_interpolation(LinkAET p, LinkAET k, int y)
{
    AET temp1 = p->x < k->x ? *p : *k;
    AET temp2 = p->x < k->x ? *k : *p;

    //temp1.x = (int)temp1.x;
    //temp2.x = (int)temp2.x;
    temp1.x = (int)(temp1.x + 0.5 * (temp1.x >= 0 ? 1 : -1));
    temp2.x = (int)(temp2.x + 0.5 * (temp2.x >= 0 ? 1 : -1));

    double z_per_x = (temp2.z - temp1.z) / (temp2.x - temp1.x);
    double z_now = temp1.z;

    for (int i = temp1.x; i <= temp2.x; i++)
    {
        //确保访问的索引在 Depth_Buffer 和 Color_Buffer 的范围内
        int buffer_x = i - lighting_start_point_x;
        int buffer_y = y - lighting_start_point_y;
        if (buffer_x >= 0 && buffer_x < lighting_width && buffer_y >= 0 && buffer_y < lighting_height)
        {
            if (z_now > Lighting_Depth_Buffer[buffer_x][buffer_y] || Lighting_Depth_Buffer[buffer_x][buffer_y] == -INT_MAX)
            {
                Lighting_Depth_Buffer[buffer_x][buffer_y] = z_now;
            }
        }
        z_now += z_per_x;
    }
}

void Only_Lighting_Buffer_fill(std::vector<my_homogeneous_point> polygon)
{
    int polygon_count = polygon.size();
    for (int i = 0; i < polygon_count; i++)
    {
        polygon[i].y = (int)(polygon[i].y + 0.5 * (polygon[i].y < 0 ? -1 : 1));
    }

    // 确定多边形的y范围
    int max = polygon[0].y;
    int min = polygon[0].y;
    for (int i = 1; i < polygon_count; i++)
    {
        max = fmax(max, polygon[i].y);
        min = fmin(min, polygon[i].y);
    }

    //构造边表
    LinkNET Net = new NET;
    LinkNET Net_clean = Net;
    LinkNET p = Net;

    for (int i = min; i <= max; i++)
    {
        LinkNET q = new NET;
        p->next_y = q;
        p = p->next_y;
        for (int j = 0; j < polygon_count; j++)
        {
            if (polygon[j].y == i)
            {
                // 判断该点相邻的两点与该点的位置关系
                // 如果相邻两点均高于该店，则插入两个边的信息
                if (polygon[(j - 1 + polygon_count) % polygon_count].y > polygon[j].y && polygon[(j + 1) % polygon_count].y > polygon[j].y)
                {
                    LinkNET Temp1 = new NET(polygon[j], polygon[(j + 1) % polygon_count]);
                    LinkNET Temp2 = new NET(polygon[j], polygon[(j - 1 + polygon_count) % polygon_count]);

                    if (Temp1->Delta_x < Temp2->Delta_x)
                    {
                        Temp1->next_x = Temp2;
                        Temp2->next_x = NULL;
                        q->next_x = Temp1;
                    }
                    else
                    {
                        Temp2->next_x = Temp1;
                        Temp1->next_x = NULL;
                        q->next_x = Temp2;
                    }
                    q = q->next_x->next_x;
                }

                // 如果相邻两点只有一点高于该店，则插入一个边的信息
                else if (polygon[(j + 1) % polygon_count].y > polygon[j].y)
                {
                    LinkNET Temp = new NET(polygon[j], polygon[(j + 1) % polygon_count]);
                    q->next_x = Temp;
                    q = q->next_x;
                }

                // 如果相邻两点只有一点高于该店，则插入一个边的信息
                else if (polygon[(j - 1 + polygon_count) % polygon_count].y > polygon[j].y)
                {
                    LinkNET Temp = new NET(polygon[j], polygon[(j - 1 + polygon_count) % polygon_count]);
                    q->next_x = Temp;
                    q = q->next_x;
                }
            }
        }
    }

    // 初始化活性边表
    LinkAET Aet = new AET;
    LinkAET Aet_clean = Aet;

    // 填充多边形
    while (Net->next_y != NULL && Net->next_y->next_y != NULL) // 确保当前的Net节点存在，确保下一行的Net节点存在。
    {
        LinkNET q = Net->next_y->next_x;

        // 插入排序遍历AET并按照x坐标值递增的顺序插入Temp。如果x值相同，比较Delta_x来决定插入位置，以保持排序稳定性。
        while (q != NULL)
        {
            LinkAET p = Aet;
            LinkAET Temp = new AET(q);
            while (p->next != NULL)
            {
                if (p->next->x > Temp->x)
                {
                    Temp->next = p->next;
                    break;
                }
                else if (p->next->x == Temp->x)
                {
                    if (p->next->Delta_x < Temp->Delta_x)
                    {
                        Temp->next = p->next->next;
                    }
                    else
                    {
                        Temp->next = p->next;
                        break;
                    }
                }
                p = p->next;
            }
            p->next = Temp;
            q = q->next_x;
        }


        // 更新AET和绘制多边形
        LinkAET temp = Aet;

        while (temp->next != NULL && temp->next->next != NULL)
        {
            int flag = 1;
            LinkAET p = temp->next;
            LinkAET k = p->next;
            LinkAET TEMP1 = NULL;
            LinkAET TEMP2 = NULL;

            Only_Lighting_Linear_interpolation(p, k, min);

            p->add();
            k->add();

            min++;
            Net = Net->next_y;

            if (min == p->y_max && min == k->y_max)
            {
                TEMP1 = p;
                TEMP2 = k;
                temp->next = k->next;
                flag = 0;
            }
            else
            {
                if (min == k->y_max)
                {
                    TEMP2 = k;
                    p->next = k->next;
                    flag = 0;
                }
                if (min == p->y_max)
                {
                    TEMP1 = p;
                    temp->next = k;
                    flag = 0;
                }
            }
            if (TEMP1 != NULL)
            {
                delete TEMP1;
            }
            if (TEMP2 != NULL)
            {
                delete TEMP2;
            }
            if (flag)
                temp = temp->next->next;
        }
    }

    while (Aet_clean->next != NULL)
    {
        LinkAET TEMP = Aet_clean;
        Aet_clean = Aet_clean->next;
        delete TEMP;
    }
    delete Aet_clean;

    while (Net_clean->next_y != NULL)
    {
        LinkNET TEMP_X = Net_clean->next_x;
        LinkNET TEMP_Y = Net_clean;
        Net_clean = Net_clean->next_y;
        delete TEMP_Y;

        if (TEMP_X != NULL)
        {
            while (TEMP_X->next_x != NULL)
            {
                LinkNET TEMP_X_TEMP = TEMP_X;
                TEMP_X = TEMP_X->next_x;
                delete TEMP_X_TEMP;
            }
            delete TEMP_X;
        }
    }
    delete Net_clean;
    return;
}

void Shadow_Mapping_new()
{
    init_camera(eye, point_lighting, obj1);

    double P[4][4];
    memset(P, 0, sizeof(double) * 16);
    P[0][0] = 1;
    P[1][1] = 1;
    P[3][3] = 1;
    P[3][2] = -1.0 / nearplane_distance;

    for (int i = 0; i < clipping_height; i++)
    {
        for (int j = 0; j < clipping_width; j++)
        {
            if (Depth_Buffer[i][j].Depth == -INT_MAX)
                continue;
            my_homogeneous_point temp1 = Depth_Buffer[i][j].pre_p;
            my_homogeneous_point temp2 = camera_change(temp1);

            my_homogeneous_point input_v = temp2;
            double temp = input_v.z;
            input_v = matrix_multiply_vector(P, input_v);
            temp2 = input_v;
            temp2 = temp2.normalize();
            temp2.z = temp;

            int x = (int)(temp2.x + 0.5 * (temp2.x < 0 ? -1 : 1));
            int y = (int)(temp2.y + 0.5 * (temp2.y < 0 ? -1 : 1));
            int buffer_x = x - lighting_start_point_x;
            int buffer_y = y - lighting_start_point_y;

            if (buffer_x >= 0 && buffer_x < lighting_width && buffer_y >= 0 && buffer_y < lighting_height)
            {
                if (Lighting_Depth_Buffer[buffer_x][buffer_y] - temp2.z > 2)
                {
                    Color_Buffer[i][j] = rgb;
                }
            }
        }
    }
}

void SMT()
{
    obj2->faces.clear();
    obj_change(obj1, obj2);

    //for (int q = 0; q < obj2->faces.size(); q++)
    //{
    //    //if (judge_face(lighting_project_obj[0].faces[q].point,now_view))
    //    {
    //        //Edge_Table_Blinn_Phong(lighting_project_obj[0].faces[q]);
    //        //Buffer_fill(lighting_project_obj[0].faces[q].point, lighting_b);
    //        Only_Lighting_Buffer_fill(obj2->faces[q].point);
    //    }
    //}

    projected(obj2, lighting_project_obj, nearplane_distance);

    for (int q = 0; q < lighting_project_obj[0].faces.size(); q++)
    {
        //if (judge_face(lighting_project_obj[0].faces[q].point,now_view))
        {
            //Edge_Table_Blinn_Phong(lighting_project_obj[0].faces[q]);
            //Buffer_fill(lighting_project_obj[0].faces[q].point, lighting_b);
            Only_Lighting_Buffer_fill(lighting_project_obj[0].faces[q].point);
        }
    }

    Shadow_Mapping_new();
}
///
void keyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 'w':
    case 'W':
    {
        my_traslate_homogeneous(obj1, 0, 5, 0);
        glutPostRedisplay();
        break;
    }
    case 's':
    case 'S':
    {
        my_traslate_homogeneous(obj1, 0, -5, 0);
        glutPostRedisplay();
        break;
    }
    case 'a':
    case 'A':
    {
        my_traslate_homogeneous(obj1, -5, 0, 0);
        glutPostRedisplay();
        break;
    }
    case 'd':
    case 'D':
    {
        my_traslate_homogeneous(obj1, 5, 0, 0);
        glutPostRedisplay();
        break;
    }
    case 'n':
    case 'N':
    {
        my_traslate_homogeneous(obj1, 0, 0, 5);
        glutPostRedisplay();
        break;
    }
    case 'f':
    case 'F':
    {
        my_traslate_homogeneous(obj1, 0, 0, -5);
        glutPostRedisplay();
        break;
    }
    case 'p':
    case 'P':
    {
        if (shadow_mapping_temp++ % 2 == 0)
        {
            rgb = { 0,1,0 };
        }
        else
        {
            rgb = { 0,0,0 };
        }
        glutPostRedisplay();
        break;
    }
    case 'r':
    case 'R':
    {
        lighting_p.r += 0.1;
        glutPostRedisplay();
        break;
    }
    case 'U':
    case 'u':
    {
        UV_SWITCH++;
        glutPostRedisplay();
        break;
    }
    case 'Y':
    case 'y':
    {
        my_rotateY_homogeneous(obj1, 1);
        glutPostRedisplay();
        break;
    }
    case 27:
        exit(0);
        break;
    }
}

void reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (w <= h)
        glFrustum(-0.5 * nearplane_width, 0.5 * nearplane_width, -0.5 * nearplane_height * (GLdouble)nearplane_height / (GLdouble)nearplane_width, 0.5 * nearplane_height * (GLdouble)nearplane_height / (GLdouble)nearplane_width,
            nearplane_distance, farplane_distance); //相对于视点
    else
        glFrustum(-0.5 * nearplane_width, 0.5 * nearplane_width, -0.5 * nearplane_height * (GLdouble)nearplane_width / (GLdouble)nearplane_height, 0.5 * nearplane_height * (GLdouble)nearplane_width / (GLdouble)nearplane_height,
            nearplane_distance, farplane_distance);


    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(eye.x, eye.y, eye.z, 0, 0, 0, 0, 1, 0);
}

void display(void)
{
    glClearColor(1.f, 1.f, 1.f, 0.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    projected(obj1, project_obj, nearplane_distance);

    my_3Dvector now_view(eye, { 0,0,0 });
    now_view = now_view.normalize();
    for (int q = 0; q < project_obj[0].faces.size(); q++)
    {
        if (judge_face(project_obj[0].faces[q].point, now_view))
        {
            Sutherland_Hodgman_Clipping(clipping_polygon, project_obj[0].faces[q].point, 4);

            Edge_Table_Blinn_Phong(project_obj[0].faces[q], obj1[0].faces[q]);
        }
    }

    SMT();

    Show();

    glColor3f(0, 0, 0);
    DDA_Common(clipping_polygon, clipping_polygon_point_count);		//裁剪窗口

    glutSwapBuffers();

}

int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(nearplane_width, nearplane_height);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("光照");

    init();

    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    //glutMouseFunc(mouse);
    glutMainLoop();
    return 0;
}
