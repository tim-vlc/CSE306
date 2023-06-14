#include "vector.cpp"

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
    std::vector<Vector> vertices;
};