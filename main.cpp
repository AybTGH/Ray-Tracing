
#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <cmath>
#include <iostream>
#include <random>
#include <list>
#include <chrono> /
#define M_PI 3.1415926535897932

std::default_random_engine engine;
std::uniform_real_distribution<double> uniform(0, 1);


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

Vector operator*(const Vector& a, const Vector& b)
{
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector operator/(const Vector& a, const double& b)
{
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

Vector random_cos(const Vector& N)
{
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    Vector direction_aleatoire_local(cos(2 * M_PI * r1) * sqrt(1 - r2), sin(2 * M_PI * r1) * sqrt(1 - r2), sqrt(r2));
    Vector aleatoire(uniform(engine) - 0.5, uniform(engine) - 0.5, uniform(engine) - 0.5);
    Vector tangent1 = cross(N, aleatoire); tangent1.normalize();
    Vector tangent2 = cross(tangent1, N);
    return direction_aleatoire_local[2] * N + direction_aleatoire_local[0] * tangent1 + direction_aleatoire_local[1] * tangent2;
}


double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

class Ray {
public:
    Ray(const Vector& o, const Vector& d) : origin(o), direction(d) {};
    Vector origin, direction;
};

class Object {
public: Object() {};
      virtual bool intersect(const Ray& d, Vector& P, Vector& N, double& t, Vector& col) const = 0;
      Vector albedo;
      bool miroir;
      bool transp;
      bool sphere_txt;
};

class Triangle : public Object {
public:
    Triangle(const Vector& A, const Vector& B, const Vector& C, const Vector& couleur, bool mirror = false, bool transp = false)
        : A(A), B(B), C(C) {
        albedo = couleur;
        miroir = mirror;
        transp = transp;
    };

    bool intersect(const Ray& d, Vector& P_out, Vector& N_out, double& t_out, Vector& col) const {
        double alpha, beta, gamma;
        col = albedo;
        return intersect(d, P_out, N_out, t_out, alpha, beta, gamma);
    }

    bool intersect(const Ray& d, Vector& P_out, Vector& N_out, double& t_out, double& alpha, double& beta, double& gamma) const {
        Vector N_calc = cross(B - A, C - A).getNormalized();
        double t_calc = dot(C - d.origin, N_calc) / dot(d.direction, N_calc);
        if (t_calc <= 0) return false;

        Vector P_calc = d.origin + t_calc * d.direction;
        Vector u = B - A;
        Vector v = C - A;
        Vector w = P_calc - A;
        double m11 = u.norm2();
        double m22 = v.norm2();
        double m12 = dot(u, v);
        double detm = m11 * m22 - m12 * m12;

        double b11 = dot(w, u);
        double b21 = dot(w, v);
        double detb = b11 * m22 - b21 * m12;
        beta = detb / detm; // Barycentric coordinate related to B

        double g12 = b11;
        double g22 = b21;
        double detg = m11 * g22 - m12 * g12;
        gamma = detg / detm; // Barycentric coordinate related to C

        alpha = 1 - beta - gamma;
        if (alpha < 0 || alpha > 1) return false;
        if (beta < 0 || beta > 1) return false;
        if (gamma < 0 || gamma > 1) return false;

        // If the intersection is valid, set the output variables
        P_out = P_calc;
        N_out = N_calc;
        t_out = t_calc;
        return true;
    }
    Vector A, B, C;
};

class BoudingBox
{
public:
    bool intersect(const Ray& r) const {

        // les plans verticaux de la boîte englobante qui sont parallèles à l'axe yz
        double tx1 = (bbMin[0] - r.origin[0]) / r.direction[0];
        double tx2 = (bbMax[0] - r.origin[0]) / r.direction[0];
        double txMin = std::min(tx1, tx2), txMax = std::max(tx1, tx2);

        // les plans horizontaux parallèles à l'axe xz, en utilisant les coordonnées y.
        double ty1 = (bbMin[1] - r.origin[1]) / r.direction[1];
        double ty2 = (bbMax[1] - r.origin[1]) / r.direction[1];
        double tyMin = std::min(ty1, ty2), tyMax = std::max(ty1, ty2);

        //  les plans parallèles à l'axe xy, en utilisant les coordonnées z.
        double tz1 = (bbMin[2] - r.origin[2]) / r.direction[2];
        double tz2 = (bbMax[2] - r.origin[2]) / r.direction[2];
        double tzMin = std::min(tz1, tz2), tzMax = std::max(tz1, tz2);

        // tMin est le plus grand des points d'entrée sur les trois paires de plans,
        // et tMax est le plus petit des points de sortie.
        double tMin = std::max(txMin, std::max(tyMin, tzMin));
        double tMax = std::min(txMax, std::min(tyMax, tzMax));
        return (tMax > 0 && tMax > tMin);
    }

    Vector bbMin, bbMax;
};

class BVH {
public:
    int i0, i1; // Exemple de membres publics
    BoudingBox bbox; // Exemple d'un objet bbox de type BBox

    BVH* fg, * fd;
};

class TriangleIndices
{
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {};
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;    // indices within the uv coordinates array
    int ni, nj, nk;       // indices within the normals array
    int group;            // face group
};

class TriangleMesh : public Object
{
public:
    ~TriangleMesh() {}
    TriangleMesh(const Vector& albedo, bool mirror = false, bool transp = false)
    {
        this->albedo = albedo;
        miroir = mirror;
        transp = transp;
    };

    void readOBJ(const char* obj)
    {

        char matfile[255];
        char grp[255];

        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f))
        {
            char line[255];
            if (!fgets(line, 255, f))
                break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's')
            {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ')
            {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6)
                {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
                }
                else
                {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n')
            {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't')
            {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f')
            {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char* consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9)
                {
                    if (i0 < 0)
                        t.vtxi = vertices.size() + i0;
                    else
                        t.vtxi = i0 - 1;
                    if (i1 < 0)
                        t.vtxj = vertices.size() + i1;
                    else
                        t.vtxj = i1 - 1;
                    if (i2 < 0)
                        t.vtxk = vertices.size() + i2;
                    else
                        t.vtxk = i2 - 1;
                    if (j0 < 0)
                        t.uvi = uvs.size() + j0;
                    else
                        t.uvi = j0 - 1;
                    if (j1 < 0)
                        t.uvj = uvs.size() + j1;
                    else
                        t.uvj = j1 - 1;
                    if (j2 < 0)
                        t.uvk = uvs.size() + j2;
                    else
                        t.uvk = j2 - 1;
                    if (k0 < 0)
                        t.ni = normals.size() + k0;
                    else
                        t.ni = k0 - 1;
                    if (k1 < 0)
                        t.nj = normals.size() + k1;
                    else
                        t.nj = k1 - 1;
                    if (k2 < 0)
                        t.nk = normals.size() + k2;
                    else
                        t.nk = k2 - 1;
                    indices.push_back(t);
                }
                else
                {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6)
                    {
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        if (j0 < 0)
                            t.uvi = uvs.size() + j0;
                        else
                            t.uvi = j0 - 1;
                        if (j1 < 0)
                            t.uvj = uvs.size() + j1;
                        else
                            t.uvj = j1 - 1;
                        if (j2 < 0)
                            t.uvk = uvs.size() + j2;
                        else
                            t.uvk = j2 - 1;
                        indices.push_back(t);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3)
                        {
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            indices.push_back(t);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            if (k0 < 0)
                                t.ni = normals.size() + k0;
                            else
                                t.ni = k0 - 1;
                            if (k1 < 0)
                                t.nj = normals.size() + k1;
                            else
                                t.nj = k1 - 1;
                            if (k2 < 0)
                                t.nk = normals.size() + k2;
                            else
                                t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true)
                {
                    if (consumedline[0] == '\n')
                        break;
                    if (consumedline[0] == '\0')
                        break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3)
                    {
                        if (i0 < 0)
                            t2.vtxi = vertices.size() + i0;
                        else
                            t2.vtxi = i0 - 1;
                        if (i2 < 0)
                            t2.vtxj = vertices.size() + i2;
                        else
                            t2.vtxj = i2 - 1;
                        if (i3 < 0)
                            t2.vtxk = vertices.size() + i3;
                        else
                            t2.vtxk = i3 - 1;
                        if (j0 < 0)
                            t2.uvi = uvs.size() + j0;
                        else
                            t2.uvi = j0 - 1;
                        if (j2 < 0)
                            t2.uvj = uvs.size() + j2;
                        else
                            t2.uvj = j2 - 1;
                        if (j3 < 0)
                            t2.uvk = uvs.size() + j3;
                        else
                            t2.uvk = j3 - 1;
                        if (k0 < 0)
                            t2.ni = normals.size() + k0;
                        else
                            t2.ni = k0 - 1;
                        if (k2 < 0)
                            t2.nj = normals.size() + k2;
                        else
                            t2.nj = k2 - 1;
                        if (k3 < 0)
                            t2.nk = normals.size() + k3;
                        else
                            t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2)
                        {
                            if (i0 < 0)
                                t2.vtxi = vertices.size() + i0;
                            else
                                t2.vtxi = i0 - 1;
                            if (i2 < 0)
                                t2.vtxj = vertices.size() + i2;
                            else
                                t2.vtxj = i2 - 1;
                            if (i3 < 0)
                                t2.vtxk = vertices.size() + i3;
                            else
                                t2.vtxk = i3 - 1;
                            if (j0 < 0)
                                t2.uvi = uvs.size() + j0;
                            else
                                t2.uvi = j0 - 1;
                            if (j2 < 0)
                                t2.uvj = uvs.size() + j2;
                            else
                                t2.uvj = j2 - 1;
                            if (j3 < 0)
                                t2.uvk = uvs.size() + j3;
                            else
                                t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2)
                            {
                                if (i0 < 0)
                                    t2.vtxi = vertices.size() + i0;
                                else
                                    t2.vtxi = i0 - 1;
                                if (i2 < 0)
                                    t2.vtxj = vertices.size() + i2;
                                else
                                    t2.vtxj = i2 - 1;
                                if (i3 < 0)
                                    t2.vtxk = vertices.size() + i3;
                                else
                                    t2.vtxk = i3 - 1;
                                if (k0 < 0)
                                    t2.ni = normals.size() + k0;
                                else
                                    t2.ni = k0 - 1;
                                if (k2 < 0)
                                    t2.nj = normals.size() + k2;
                                else
                                    t2.nj = k2 - 1;
                                if (k3 < 0)
                                    t2.nk = normals.size() + k3;
                                else
                                    t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            }
                            else
                            {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1)
                                {
                                    if (i0 < 0)
                                        t2.vtxi = vertices.size() + i0;
                                    else
                                        t2.vtxi = i0 - 1;
                                    if (i2 < 0)
                                        t2.vtxj = vertices.size() + i2;
                                    else
                                        t2.vtxj = i2 - 1;
                                    if (i3 < 0)
                                        t2.vtxk = vertices.size() + i3;
                                    else
                                        t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                }
                                else
                                {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(f);
        build_bvh(&bvh, 0, indices.size());
    }
    void add_textures(const char* filename) {
        int w, h, c;
        unsigned char* texture = stbi_load(filename, &w, &h, &c, 3);
        if (texture != nullptr) {
            textures.push_back(texture);
            textures_widths.push_back(w);
            textures_heights.push_back(h);
            std::cout << "Texture ajoutée : " << filename << ", Largeur : " << w << ", Hauteur : " << h << std::endl;
        }
        else {
            std::cerr << "Erreur lors du chargement de la texture : " << filename << std::endl;
        }
    }

    BoudingBox build_bbox(int i0, int i1)
    {
        BoudingBox result;
        result.bbMax = vertices[indices[i0].vtxi];
        result.bbMin = vertices[indices[i0].vtxj];

        for (int i = i0; i < i1; i++) { // Indice de triangle
            result.bbMin[0] = std::min(result.bbMin[0], vertices[indices[i].vtxi][0]);
            result.bbMin[1] = std::min(result.bbMin[1], vertices[indices[i].vtxi][1]);
            result.bbMin[2] = std::min(result.bbMin[2], vertices[indices[i].vtxi][2]);
            result.bbMax[0] = std::max(result.bbMax[0], vertices[indices[i].vtxi][0]); // itérer pour i, j, k
            result.bbMax[1] = std::max(result.bbMax[1], vertices[indices[i].vtxi][1]);
            result.bbMax[2] = std::max(result.bbMax[2], vertices[indices[i].vtxi][2]);

            result.bbMin[0] = std::min(result.bbMin[0], vertices[indices[i].vtxj][0]);
            result.bbMin[1] = std::min(result.bbMin[1], vertices[indices[i].vtxj][1]);
            result.bbMin[2] = std::min(result.bbMin[2], vertices[indices[i].vtxj][2]);
            result.bbMax[0] = std::max(result.bbMax[0], vertices[indices[i].vtxj][0]); // itérer pour i, j, k
            result.bbMax[1] = std::max(result.bbMax[1], vertices[indices[i].vtxj][1]);
            result.bbMax[2] = std::max(result.bbMax[2], vertices[indices[i].vtxj][2]);

            result.bbMin[0] = std::min(result.bbMin[0], vertices[indices[i].vtxk][0]);
            result.bbMin[1] = std::min(result.bbMin[1], vertices[indices[i].vtxk][1]);
            result.bbMin[2] = std::min(result.bbMin[2], vertices[indices[i].vtxk][2]);
            result.bbMax[0] = std::max(result.bbMax[0], vertices[indices[i].vtxk][0]); // itérer pour i, j, k
            result.bbMax[1] = std::max(result.bbMax[1], vertices[indices[i].vtxk][1]);
            result.bbMax[2] = std::max(result.bbMax[2], vertices[indices[i].vtxk][2]);
        }
        return result;
    }
    void build_bvh(BVH* node, int i0, int i1) {
        node->i0 = i0;
        node->i1 = i1;
        node->bbox = build_bbox(node->i0, node->i1);
        Vector diag = node->bbox.bbMax - node->bbox.bbMin;

        int split_dim = (diag[0] > diag[1]) ? ((diag[0] > diag[2]) ? 0 : 2) : ((diag[1] > diag[2]) ? 1 : 2);
        double split_val = 0.5 * (node->bbox.bbMin[split_dim] + node->bbox.bbMax[split_dim]);

        int pivot = node->i0;
        for (int i = node->i0; i < node->i1; ++i) {
            double mid_point = (vertices[indices[i].vtxi][split_dim] + vertices[indices[i].vtxj][split_dim] + vertices[indices[i].vtxk][split_dim]) / 3.0;
            if (mid_point < split_val) {
                std::swap(indices[i], indices[pivot]);
                pivot++;
            }
        }

        node->fg = nullptr;
        node->fd = nullptr;
        if (pivot == node->i0 || pivot == node->i1 || node->i1 - node->i0 < 5) {
            return;
        }
        node->fg = new BVH();
        build_bvh(node->fg, node->i0, pivot);
        node->fd = new BVH();
        build_bvh(node->fd, pivot, node->i1);
    }

    void rescale(double scale) {
        for (int i = 0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] * scale; // Multiplie chaque composante du vecteur par 'scale'.
        }
        build_bvh(&bvh, 0, indices.size());
    }

    void translate(const Vector& translation) {
        for (int i = 0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] + translation;
        }
        build_bvh(&bvh, 0, indices.size());
    }

    void rotateY(double angle) {
        // Conversion de l'angle de degrés en radians
        double radians = angle * M_PI / 180.0;

        // Matrice de rotation autour de l'axe Y
        for (int i = 0; i < vertices.size(); i++) {
            double x = vertices[i][0];
            double y = vertices[i][1];
            double z = vertices[i][2];

            vertices[i][0] = x * cos(radians) + z * sin(radians);
            vertices[i][2] = -x * sin(radians) + z * cos(radians);
        }

        // Reconstruire la boîte englobante après rotation
        build_bvh(&bvh, 0, indices.size());
    }

    void rotateX(double angle) {
        double radians = angle * M_PI / 180.0;

        for (int i = 0; i < vertices.size(); i++) {
            double x = vertices[i][0];
            double y = vertices[i][1];
            double z = vertices[i][2];

            // Rotation autour de l'axe X
            vertices[i][1] = y * cos(radians) - z * sin(radians);
            vertices[i][2] = y * sin(radians) + z * cos(radians);
        }
        build_bvh(&bvh, 0, indices.size());

    }

    void rotateZ(double angle) {
        double radians = angle * M_PI / 180.0;

        for (int i = 0; i < vertices.size(); i++) {
            double x = vertices[i][0];
            double y = vertices[i][1];
            double z = vertices[i][2];

            // Rotation autour de l'axe Z
            vertices[i][0] = x * cos(radians) - y * sin(radians); // Nouvelle coordonnée x
            vertices[i][1] = x * sin(radians) + y * cos(radians); // Nouvelle coordonnée y
            // z reste inchangé
        }
        build_bvh(&bvh, 0, indices.size());
    }
    bool intersect(const Ray& r, Vector& P, Vector& N, double& t, Vector& col) const {
        t = 1E99;
        bool has_inter = false;

        if (!bvh.bbox.intersect(r))
            return false;

        std::list<const BVH*> l;
        l.push_front(&bvh);

        while (!l.empty()) {
            const BVH* current = l.front();
            l.pop_front();

            if (current->fg && current->fg->bbox.intersect(r)) {
                l.push_back(current->fg);
            }

            if (current->fd && current->fd->bbox.intersect(r)) {
                l.push_back(current->fd);
            }

            if (!current->fg) {
                for (int i = current->i0; i < current->i1; i++) {
                    // Récupère les points A, B et C du triangle courant.
                    int a = indices[i].vtxi;
                    int b = indices[i].vtxj;
                    int c = indices[i].vtxk;
                    Triangle tri(vertices[a], vertices[b], vertices[c], albedo, miroir, transp); // Assuming albedo, moroir, transparent are accessible here
                    Vector localP, localN;
                    double localt;
                    double alpha, beta, gamma;
                    if (tri.intersect(r, localP, localN, localt, alpha, beta, gamma)) {
                        has_inter = true;
                        if (localt < t) {
                            t = localt;
                            P = localP;
                            //N = localN;
                            N = alpha * normals[indices[i].ni] + beta * normals[indices[i].nj] + gamma * normals[indices[i].nk];
                            N.getNormalized();
                            int H = textures_heights[indices[i].group];
                            int W = textures_widths[indices[i].group];
                            Vector UV = alpha * uvs[indices[i].uvi] + beta * uvs[indices[i].uvj] + gamma * uvs[indices[i].uvk];
                            UV = UV * Vector(W, H, 0);
                            int uvx = static_cast<int>(UV[0] + 0.5) % W;
                            int uvy = static_cast<int>(UV[1] + 0.5) % H;
                            //  negative indices
                            // uvx = (uvx < 0) ? uvx + W : uvx;
                            // uvy = (uvy < 0) ? uvy + H : uvy;
                            uvx = (uvx < 0) ? uvx + W : uvx;
                            uvy = (uvy < 0) ? uvy + H : uvy;
                            uvy = H - uvy - 1;
                            int textureIndex = (uvy * W + uvx) * 3;
                            col = Vector(std::pow(textures[indices[i].group][(uvy * W + uvx) * 3] / 255., 2.2),
                                std::pow(textures[indices[i].group][(uvy * W + uvx) * 3 + 1] / 255., 2.2),
                                std::pow(textures[indices[i].group][(uvy * W + uvx) * 3 + 2] / 255., 2.2));
                            //try {
                            //    t = localt;
                            //    P = localP;
                            //    //N = localN;
                            //    N = alpha * normals[indices[i].ni] + beta * normals[indices[i].nj] + gamma * normals[indices[i].nk];
                            //    N.getNormalized();



                            //    Vector UV = alpha * uvs[indices[i].uvi] + beta * uvs[indices[i].uvj] + gamma * uvs[indices[i].uvk];
                            //    UV.normalize();
                            //    //std::cout << "Vecteur uv : (" << UV[0] << ", " << UV[1] << ", " << UV[2] << ")" << std::endl;

                            //    int x = fabs(UV[0]) * (textures_widths[indices[i].group] - 1);
                            //    int y = fabs(UV[1]) * (textures_heights[indices[i].group] - 1);

                            //    // Vérification des indices
                            //    if (x < 0 || y < 0 || x >= textures_widths[indices[i].group] || y >= textures_heights[indices[i].group]) {
                            //        throw std::out_of_range("Indices de texture en dehors des limites");
                            //    }

                            //    double c1 = textures[indices[i].group][(y * textures_widths[indices[i].group] + x) * 3 + 0] / 255.;
                            //    double c2 = textures[indices[i].group][(y * textures_widths[indices[i].group] + x) * 3 + 1] / 255.;
                            //    double c3 = textures[indices[i].group][(y * textures_widths[indices[i].group] + x) * 3 + 2] / 255.;
                            //    col = Vector(1, 0, 0);
                            //    std::cout << "Vecteur col : (" << col[0] << ", " << col[1] << ", " << col[2] << ")" << std::endl;
                            //}
                            //catch (const std::out_of_range& e) {
                            //    std::cerr << "Erreur lors de l'accès à la texture : " << e.what() << std::endl;
                            //    // Attribuer le vecteur (1, 0, 0) à col en cas d'erreur
                            //    col = Vector(1, 0, 0);
                            //}
                            //catch (const std::exception& e) {
                            //    std::cerr << "Erreur inattendue : " << e.what() << std::endl;
                            //}
                                //int tW = textures_widths[indices[i].group];
                                //int tH = textures_heights[indices[i].group];
                                //Vector UV = alpha * uvs[indices[i].uvi] + beta * uvs[indices[i].uvj] + gamma * uvs[indices[i].uvk];
                                //int x = UV[0] * (tW -1);
                                //int y = UV[1] * (tH - 1);

                                //std::cout << "Valeur de x : " << x << std::endl;
                                //std::cout << "Valeur de y : " << y << std::endl;
                                //// Vérification des indices de texture
                                //if (x < 0 || y < 0 || x >= textures_widths[indices[i].group] || y >= textures_heights[indices[i].group]) {
                                //    // Gérer l'erreur ici (par exemple, attribuer une valeur par défaut à la couleur de texture)
                                //    col = Vector(1, 0, 0); // Couleur par défaut (rouge)
                                //    std::cerr << "Erreur : Indices de texture en dehors des limites" << std::endl;

                                //    std::cout << "Valeur de x : " << x << std::endl;
                                //    std::cout << "Valeur de y : " << y << std::endl;
                                //}
                                //else {
                                //    // Accès à la texture avec les indices valides
                                //    double c1 = textures[indices[i].group][(y * textures_widths[indices[i].group] + x) * 3 + 0] / 255.;
                                //    double c2 = textures[indices[i].group][(y * textures_widths[indices[i].group] + x) * 3 + 1] / 255.;
                                //    double c3 = textures[indices[i].group][(y * textures_widths[indices[i].group] + x) * 3 + 2] / 255.;
                                //    col = Vector(c1, c2, c3);
                                //    //std::cout << "Vecteur col : (" << col[0] << ", " << col[1] << ", " << col[2] << ")" << std::endl;
                                //}

                        }



                    }
                }
            }
        }
        return has_inter;
    }

    std::vector<unsigned char*> textures;
    std::vector<int> textures_widths;
    std::vector<int> textures_heights;

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
private:
    BVH bvh;
};

class Sphere : public Object
{
public:
    Sphere(const Vector& centre, double Ray, const Vector& albedo, bool miroir = false, bool transp = false, bool sphere_txt = false) : O(centre), R(Ray)
    {
        this->albedo = albedo;
        this->miroir = miroir;
        this->transp = transp;
        this->sphere_txt = sphere_txt;
    }
    bool intersect(const Ray& d, Vector& P, Vector& N, double& t, Vector& col) const
    {
        // solves a*t^2 + b*t + c = 0
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
        if (sphere_txt)
        {
            Vector Pc = P - O;
            double phi = atan2(Pc[2], Pc[0]);
            double theta = acos(Pc[1] / R);
            double u = (phi + M_PI) / (2 * M_PI);
            double v = (M_PI - theta) / M_PI;
            int checkWidth = 20;
            int checkHeight = 20;
            int i = (int)(u * checkWidth);
            int j = (int)(v * checkHeight);
            if ((i + j) % 2 == 0) {
                col = Vector(1, 1, 1);
                // color = Vector(1, 0.2, 0.4); 
            }
            else {
                col = Vector(0, 0, 0);
                // Vector lightBlue(0x3D / 255.0, 0xDD / 255.0, 0xD6 / 255.0); 
                // color = lightBlue; 
            }
        }
        else {
            col = albedo;
        }

        return true;
    }
    Vector O;
    double R;
};
const int MAX_REBOUNDS = 5;

class Scene
{
public:
    Scene() {};

    void addObject(Object& o) { objects.push_back(&o); };

    std::vector<Object*> objects;
    Sphere* lumiere;
    double intensite_lumiere;
    bool intersect(const Ray& r, Vector& P, Vector& N, Vector& albedo, bool& mirror, bool& transp, double& t, int& objectId, Vector& col) const
    {
        t = 1E16;
        bool has_inter = false;
        for (int i = 0; i < objects.size(); i++)
        {
            Vector localP, localN, localCol;
            double localt;
            if (objects[i]->intersect(r, localP, localN, localt, localCol) && localt < t)
            {
                t = localt;
                has_inter = true;
                albedo = objects[i]->albedo;
                P = localP;
                N = localN;
                mirror = objects[i]->miroir;
                transp = objects[i]->transp;
                objectId = i;
                col = localCol;
            }
        }
        return has_inter;
    };
};
Vector getColor(const Ray& r, const Scene& s, int rebond, bool wasLastBounceDiffuse)
{
    double epsilon = 0.00001;
    Vector P, N, albedo, col;
    double t;
    bool mirror, transp;
    int objectId;
    bool inter = s.intersect(r, P, N, albedo, mirror, transp, t, objectId, col);
    Vector color(0, 0, 0);
    if (rebond > MAX_REBOUNDS)
        return Vector(0., 0., 0.);
    if (inter)
    {
        if (objectId == 0)
        {
            if (rebond == 0 || !wasLastBounceDiffuse)
            {
                return Vector(1., 1., 1.) * s.intensite_lumiere / (4 * M_PI * M_PI * s.lumiere->R * s.lumiere->R);
            }
            return Vector(0., 0., 0.);
        }
        if (mirror)
        {
            Vector direction_miroir = r.direction - 2 * dot(r.direction, N) * N;
            Ray Ray_miroir(P + epsilon * N, direction_miroir);
            return getColor(Ray_miroir, s, rebond + 1, false);
        }
        if (transp)
        {
            // calcul de la tarnsparence
            double n1 = 1, n2 = 1.3;
            Vector N_refract = N;
            double cosI = dot(N, r.direction);
            if (cosI > 0)
            { //on sort de la sphère
                std::swap(n1, n2);
                N_refract = -N;
            }

            double rad = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.direction, N_refract)));
            if (rad < 0) //Ray reflechie uniquement
            {
                Vector reflectedDir = r.direction - 2 * cosI * N;
                Ray reflectedRay(P + epsilon * N, reflectedDir);
                return getColor(reflectedRay, s, rebond + 1, false);
            }

            // implementation de la transmission de Fresnel
            double k0 = sqr(n1 - n2) / sqr(n1 + n2);
            double R = k0 + (1.0 - k0) * sqr(sqr(1 - std::abs(dot(N_refract, r.direction)))) * (1 - std::abs(dot(N_refract, r.direction)));
            double random = rand() % 100;
            if (random <= 100 * R)
            {
                Vector reflectedDir = r.direction - 2 * cosI * N;
                Ray reflectedRay(P + epsilon * N, reflectedDir);
                return getColor(reflectedRay, s, rebond + 1, false);
            }
            else
            {
                Vector t_T = n1 / n2 * (r.direction - dot(r.direction, N_refract) * N_refract);
                Vector t_N = -sqrt(rad) * N_refract;
                Vector refractedDir = t_T + t_N;
                return getColor(Ray(P - epsilon * N_refract, refractedDir), s, rebond + 1, false);
            }
        }
        else
        {
            // eclairage direct
            Vector axePO = (P - s.lumiere->O).getNormalized();
            Vector dir_aleatoire = random_cos(axePO);
            Vector point_aleatoire = dir_aleatoire * s.lumiere->R + s.lumiere->O;
            Vector wi = point_aleatoire - P;
            double d_light2 = sqrt(wi.norm2()); // recupere la  distance entre P et xprime
            wi = wi / d_light2; // normlise
            Vector P_prime, N_prime, albedo_prime, col_light;
            double t_light;
            int objectId;
            bool Mirror_prime, transp_prime;
            Ray Ray_prime(P + epsilon * N, wi);
            bool Inter_prime = s.intersect(Ray_prime, P_prime, N_prime, albedo_prime, Mirror_prime, transp_prime, t_light, objectId, col_light);
            // vérifier s'il n'est pas dans l'ombre (facteur de visibilité)
            if (Inter_prime && t_light < d_light2 * 0.99)
            {
                color = Vector(0., 0., 0.);
            }
            else
            {
                // calcul de l'éclairage direct en tenant compte de la source étendue avec la formule du cours
                Vector BRDF = col / M_PI;
                double R2 = sqr(s.lumiere->R);
                double proba = std::max(1E-8, dot(axePO, dir_aleatoire)) / (M_PI * s.lumiere->R);
                double J = std::max(0., dot(dir_aleatoire, -wi)) / (d_light2 * d_light2);
                double I = s.intensite_lumiere / (M_PI * R2);
                color = (I * BRDF * std::max(0., dot(wi, N)) * J) / proba;
            }

            // eclairage indirect
            Vector direction_aleatoire = random_cos(N);
            Ray rayon_aleatoire(P + epsilon * N, direction_aleatoire);
            color += col * getColor(rayon_aleatoire, s, rebond + 1, true);
        }
    }
    return color;
}



int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    float ini_time = clock();
    int W = 512;
    int H = 512;


    Vector C(0, 0, 55);
    Scene s;
    // sphere s1(Vector(-10, 35, -90), 18, Vector(1, 1, 1), false, true);

    Vector lightBlue(0x3D / 255.0, 0xDD / 255.0, 0xD6 / 255.0);     // #3DDDD6

    // x c'est la largeur
    // y c'est la hauteur
    // z c'est la profondeur

    Sphere s_lum(Vector(15, 40, 0), 15, Vector(1., 1., 1.));

    Sphere s0(Vector(-30, -2, -55), 15, Vector(1, 0, 0), false, true);  // Sphère à gauche
    Sphere s1(Vector(10, -2, -30), 20, Vector(1, 0, 0), true);  // Sphère à droite
    Sphere s1_(Vector(20, -2, -10), 10, Vector(1.0, 0.5, 0.5), false, false, true);  // Sphère à droite


    Sphere s2(Vector(0, -2000 - 20, 0), 2000, Vector(1, 1, 1), false, false); //sol
    Sphere s3(Vector(0, 2000 + 100, 0), 2000, Vector(1, 1, 1)); // plafond
    Sphere s4(Vector(-2000 - 50, 0, 0), 2000, Vector(0.5, 0.5, 0.5)); // mur gauche
    Sphere s5(Vector(2000 + 50, 0, 0), 2000, Vector(0.5, 0.5, 0.5)); // mur droit
    Sphere s6(Vector(0, 0, -2000 - 100), 2000, Vector(0, 1, 1)); // mur fond
    Sphere s7(Vector(0, 0, 2000 + 100), 2000, Vector(1, 1, 0)); // mur arrière caméra

    TriangleMesh m1(Vector(1., 1., 1.), false, false);
    m1.readOBJ("models/Chick.obj");

    //Triangle tri(Vector(-10, -2, -20), Vector(10, -10, -20), Vector(0, 10, -20), Vector(1, 0, 0));
    double scale1 = 0.2; 
    double scale2 = 0.2;
    m1.rescale(scale2);
    Vector translation2(-10, -10, -10);
    double angle = 90;
    m1.translate(translation2);
    m1.rotateX(-70);
    //m1.rotateY(-90);
    m1.add_textures("models/birdDiffuseMap.jpg");
    s.addObject(s_lum);
    //s.addTriangle(tri);
    //s.addObject(s1);
    s.addObject(s1_);
    //s.addObject(s0);
    s.addObject(s2);
    s.addObject(s3);
    s.addObject(s4);
    s.addObject(s5);
    s.addObject(s6);
    s.addObject(s7);
    s.addObject(m1);

    s.lumiere = &s_lum;
    s.intensite_lumiere = 2000000000;  // Brighter light source

    double fov = 60 * M_PI / 240;
    int nbrays = 100;
    int progress_counter = 0;
    double focus_distance = 55;
    std::vector<unsigned char> image(W * H * 3, 0);

#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            Vector color(0, 0, 0);
            for (int k = 0; k < nbrays; k++)
            {
                // methode de Box Muller 
                double r1 = uniform(engine), r2 = uniform(engine);
                double g1 = 0.25 * cos(2 * M_PI * r1) * sqrt(-2 * log(r2));
                double g2 = 0.25 * sin(2 * M_PI * r1) * sqrt(-2 * log(r2));
                // chaque Ray légèrement décalé par rapport au centre du pixek
                Vector u(j - W / 2 + g1 + 0.5, i - H / 2 + g2 - 0.5, -W / (2. * tan(fov / 2)));
                u = u.getNormalized();

                // depth of field
                double r3 = uniform(engine), r4 = uniform(engine);
                double g3 = 0.01 * cos(2 * M_PI * r3) * sqrt(-2 * log(r4));
                double g4 = 0.01 * sin(2 * M_PI * r3) * sqrt(-2 * log(r4));
                Vector destination = C + focus_distance * u;
                Vector new_C = C + Vector(g3, g4, 0);
                Vector new_u = (destination - new_C).getNormalized();
                Ray r(new_C, new_u);

                color += getColor(r, s, 5, false);
            }

            color = color / nbrays;

            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 1 / 2.2));
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 1 / 2.2));
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 1 / 2.2));
        }

#pragma omp atomic
        ++progress_counter;
#pragma omp critical
        {
            std::cout << "\rRendering progress: " << (progress_counter * 100.0 / H) << "%" << std::flush;
        }
        // Termine la mesure du temps et calcule la durée
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << std::endl << "Time elapsed so far: " << elapsed.count() << " seconds" << std::endl;
    }

    // wwriting the image on a png file
    stbi_write_png("Results/scene_finale.png", W, H, 3, &image[0], 0);

    // writing execution time on a txt file
    //std::ofstream execution_file;
    //execution_file.open("execution_time.txt");
    //execution_file << (clock() - ini_time) / CLOCKS_PER_SEC;
    //execution_file.close();

    return 0;
};
