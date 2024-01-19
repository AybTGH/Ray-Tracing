/*

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.1415926535897932

#define albedo 1
#include <iostream>


static inline double sqr(double x) { return x * x; }

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
    }
    double& operator[](int i) { return coord[i]; }
    double operator[](int i) const { return coord[i]; }

    Vector& operator+=(const Vector& v) {
        coord[0] += v[0];
        coord[1] += v[1];
        coord[2] += v[2];
        return *this;
    }

    double norm2() const {
        return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
    }
    void normalize() {
        double norm = sqrt(norm2());
        coord[0] /= norm;
        coord[1] /= norm;
        coord[2] /= norm;

    }
    Vector getNormalized() {
        Vector result(*this);
        result.normalize();
        return result;
    }

    double coord[3];
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector& a, double b) {
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator*(double a, const Vector& b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

class Ray {
public:
    Ray(const Vector& o, const Vector& d) : origin(o), direction(d) {};
    Vector origin, direction;
};

class Sphere {
public:
    Sphere(const Vector& origin, double rayon) : O(origin), R(rayon) {};
    Vector O;
    double R;
};
/* bool intersection(const Ray& d, const Sphere& s) {
    // resout a*t*2 + b*t +c =c0

    double a = 1;
    double b = 2 * dot(d.direction, d.origin - s.O);
    double c = (d.origin - s.O).norm2() - s.R * s.R;

    double delta = b * b - 4 * a * c;
    if (delta < 0) return false;
    double t1 = (-b - sqrt(delta)) / 2 * a;
    double t2 = (-b + sqrt(delta)) / 2 * a;

    if (t2 > 0) return true;
    if (t2 > 0) return true;
    return false;
}
bool intersection(const Ray& d, const Sphere& s, Vector& P, Vector& N) {
    // resout a*t*2 + b*t +c =c0

    double a = 1;
    double b = 2 * dot(d.direction, d.origin - s.O);
    double c = (d.origin - s.O).norm2() - s.R * s.R;

    double delta = b * b - 4 * a * c;
    if (delta < 0) return false;
    double t1 = (-b - sqrt(delta)) / 2 * a;
    double t2 = (-b + sqrt(delta)) / 2 * a;

    if (t2 > 0) return true;
    double t;
    if (t1 > 0)
        t = t1;
    else
        t = t2;

    P = d.origin + t * d.direction;
    N = (P - s.O).getNormalized();
    return false;
}

int main() {
    int W = 512;
    int H = 512;
    double fov = 60 * M_PI / 100;
     
    Sphere s(Vector(0, 0, -55), 20);

    //Vector camera(0, 0, 55);
    Vector L(-10, 20, 40);
    double I = 1E8;

    std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector direction(j - W / 2 + 0.5, -i + H / 2 - 0.5, -W / (2 * tan(fov / 2)));
            direction.normalize();
            Ray r(Vector(0, 0, 0), direction);
            Vector P, N;
            bool has_inter = intersection(r, s, P,N);
            std::cout << "bool: " << std::boolalpha << has_inter << std::endl;
            if (has_inter) {
                Vector vecLum = L - P;
                double d2 = vecLum.norm2();
                vecLum.getNormalized();
                Vector color (albedo*(I * std::max(0., dot(vecLum, N) / (4 * M_PI * d2))));


                image[(i * W + j) * 3 + 0] = color[0];   // RED
                image[(i * W + j) * 3 + 1] = color[1];  // GREEN
                image[(i * W + j) * 3 + 2] = color[2];  // BLUE
            }
        }
    }
    stbi_write_png("image2.png", W, H, 3, &image[0], 0);

    return 0;
}
*/