#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <cmath>
#include "lbfgs.c"
#include <sstream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define VOLUME_AIR 0.7
#define VOLUME_FLUID 0.3

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }

    void operator+=(const Vector& b) {
        data[0] += b[0];
        data[1] += b[1];
        data[2] += b[2];
    }

    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector& a) {
    return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}


// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
    double area() const {
        if (vertices.size() < 3) return 0;
        double result = 0;
        for (int i = 0; i < vertices.size(); i++) {
            const Vector& A = vertices[i];
            const Vector& B = vertices[(i + 1) % vertices.size()];
            result += A[0] * B[1] - A[1] * B[0];
        }
        return result / 2;
    }

    double integrate_squared_distance(const Vector& P) const {
        if (vertices.size() < 3) return 0;
        double value = 0;
        for (int i=1; i<vertices.size() - 1; i++) {
            Vector triangle[3] = {vertices[0], vertices[i], vertices[i+1]};
            double local_value = 0;

            for (int k=0; k<3; k++) {
                for (int l=k; l<3; l++) {
                    local_value += dot(triangle[k] - P, triangle[l] - P);
                }
            }
            Vector e1 = triangle[1] - triangle[0];
            Vector e2 = triangle[2] - triangle[0];
            double area_triangle = 1/2 * abs(e1[1]*e2[0] - e1[0]*e2[1]); 
            value += local_value / 6 * area_triangle;
        }
        return value;
    }

    Vector centroid() const {
        Vector c(0, 0, 0);
        int N = vertices.size();
        for (int i = 0; i < N; i++) {
            c[0] += (vertices[i][0] + vertices[(i + 1) % N][0])* (vertices[i][0] * vertices[(i + 1) % N][1] - vertices[(i + 1) % N][0] * vertices[i][1]);
            c[1] += (vertices[i][1] + vertices[(i + 1) % N][1])* (vertices[i][0] * vertices[(i + 1) % N][1] - vertices[(i + 1) % N][0] * vertices[i][1]);
        }
        return c / (6 * area());
    }
    std::vector<Vector> vertices;
};