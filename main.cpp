#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <cmath>

#define M_PI 3.1415926535897932


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
    Vector operator-() const {
        return Vector(-coord[0], -coord[1], -coord[2]);
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
    Sphere(const Vector& origin, double rayon, const Vector& couleur, bool mirror = false, bool transp = false) : O(origin), R(rayon), albedo(couleur), miroir(mirror), transparent(transp) {};
    Vector O;
    double R;
    Vector albedo;
    bool miroir;
    bool transparent;

    bool intersection(const Ray& d, Vector& P, Vector& N, double& t) const{
        // resout a*t*2 + b*t +c =c0

        double a = 1;
        double b = 2 * dot(d.direction, d.origin - O);
        double c = (d.origin - O).norm2() - R * R;

        double delta = b * b - 4 * a * c;
        if (delta < 0) return false;
        double t1 = (-b - sqrt(delta)) / 2 * a;
        double t2 = (-b + sqrt(delta)) / 2 * a;

        if (t2 < 0) return false;
        if (t1 > 0)
            t = t1;
        else
            t = t2;

        P = d.origin + t * d.direction;
        N = (P - O).getNormalized();
        return true;
    }
};


class Scene {
public:
    Scene() {};
    void addSphere(const Sphere& s) { spheres.push_back(s); }
    bool intersection(const Ray& d, Vector& P, Vector& N, int& sphere_id, double& min_t) const {
        bool has_inter = false;
        min_t = 1E99;
        for (int i = 0; i < spheres.size(); ++i) {
            Vector localP, localN;
            double t;
            bool local_has_inter = spheres[i].intersection(d, localP, localN, t);
            if (local_has_inter) {
                has_inter = true;
                if (t < min_t) {
                    min_t = t;
                    P = localP;
                    N = localN;
                    sphere_id = i;
                }
            }
        }
        return has_inter;
    }
    std::vector<Sphere> spheres;
    Vector position_lumiere;
    double intensite_lumiere;
};


Vector getColor(Ray &r, const Scene &s, int nbrebonds) {

    if (nbrebonds == 0) return Vector(0, 0, 0);

    Vector P, N;
    int sphere_id;
    double t;
    bool has_inter = s.intersection(r, P, N, sphere_id, t);

    Vector intensite_pix(0, 0, 0);
    if (has_inter) {

        if (s.spheres[sphere_id].transparent) {
            double n1 = 1;
            double n2 = 1.3;

            double k0 = pow((n1 - n2) / (n1 + n2), 2);
            double R = k0 + (1 - k0) * pow(1 - dot(N, -r.direction), 5);
            double T = 1 - R;

            Vector normale_pour_transparence(N);
            if (dot(r.direction, N) > 0) { // on sort de la sphere
                n1 = 1.3;
                n2 = 1;
                normale_pour_transparence = - N;
            }
            
            double radical = 1 - sqr(n1 / n2) * (1 - sqr(dot(normale_pour_transparence, r.direction)));
            if (radical > 0) {
                Vector direction_refraction = (n1 / n2) * (r.direction - dot(r.direction, normale_pour_transparence) * normale_pour_transparence) - normale_pour_transparence * sqrt(radical);
                Ray rayon_refracte(P - 0.01 * normale_pour_transparence, direction_refraction);
                intensite_pix = getColor(rayon_refracte, s, nbrebonds - 1);
            }
        }
        else if (s.spheres[sphere_id].miroir) {
            Vector direction_mirroir = r.direction - 2 * dot(N, r.direction) * N;
            Ray rayon_mirroir(P + 0.01 * N, direction_mirroir);
            intensite_pix = getColor(rayon_mirroir, s, nbrebonds - 1);

        }
        else {

            Ray ray_light(P + 0.01 * N, (s.position_lumiere - P).getNormalized());
            Vector P_light, N_light;
            int sphere_id_light;
            double t_light;
            bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, t_light);
            double d_light2 = (s.position_lumiere - P).norm2();
            if (has_inter_light && t_light * t_light < d_light2) {
                intensite_pix = Vector(0, 0, 0);
            }
            else {
                intensite_pix = s.spheres[sphere_id].albedo * (s.intensite_lumiere * std::max(0., dot((s.position_lumiere - P).getNormalized(), N)) / (s.position_lumiere - P).norm2());
            }
        }
    }
    return intensite_pix;
}

int main() {
    int W = 512;
    int H = 512;
    double fov = 60 * M_PI / 100;

    Sphere s1(Vector(-15, 0, -55), 10, Vector(1, 0, 0), false, true);
    Sphere s1bis(Vector(15, 0, -55), 10, Vector(1, 0, 0), true);
    Sphere s2(Vector(0, -2000 - 20, 0), 2000, Vector(1, 1, 1)); //sol
    Sphere s3(Vector(0, 2000 + 100, 0), 2000, Vector(1, 1, 1)); // plafond
    Sphere s4(Vector(-2000 - 50, 0, 0), 2000, Vector(0, 1, 0)); // mur gauche
    Sphere s5(Vector(2000 + 50, 0, 0), 2000, Vector(0, 0, 1)); // mur droit
    Sphere s6(Vector(0, 0, -2000 - 100), 2000, Vector(0, 1, 1)); // mur fond

    Scene s;
    s.addSphere(s1);
    s.addSphere(s1bis);
    s.addSphere(s2);
    s.addSphere(s3);
    s.addSphere(s4);
    s.addSphere(s5);
    s.addSphere(s6);
    s.position_lumiere = Vector(15, 60, -40);
    s.intensite_lumiere = 100000000;

    std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector direction(j - W / 2 + 0.5, -i + H / 2 - 0.5, -W / (2 * tan(fov / 2)));
            direction.normalize();
            Ray r(Vector(0, 0, 0), direction);

            Vector color = getColor(r, s, 5);

            // Flip vertically: use i instead of (H - i - 1) to reverse the pixel order
            image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(color[0], 1 / 2.2))); // RED
            image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(color[1], 1 / 2.2))); // GREEN
            image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(color[2], 1 / 2.2)));  // BLUE
        }
    }
    stbi_write_png("Results/image10.png", W, H, 3, &image[0], 0);

    return 0;
}