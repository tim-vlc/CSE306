#include "save_image.cpp"

class PowerDiagram {
public:
    PowerDiagram() {}

    Polygon disk;

    PowerDiagram(std::vector<Vector>& pts, const std::vector<double>& weights) {
        points = pts;
        this->weights = weights;
        const int N_disk = 50;
        disk.vertices.resize(N_disk);

        for (int i=0; i<N_disk; i++) {
            double t = i / (double)N_disk * M_PI * 2;
            disk.vertices[i][0] = cos(t);
            disk.vertices[i][1] = sin(t);
            disk.vertices[i][2] = 0;
        }
    }

    //Sutherland Hodgman Algorithm
    Polygon clip_polygon_by_bissector(const Polygon &poly, int index_0, int index_i, const Vector& P0, const Vector& Pi) const {
        Polygon result;
        Vector M = (P0 + Pi) / 2;
        Vector Mprime = M + (weights[index_0] - weights[index_i]) / (2. * (P0 - Pi).norm2()) * (Pi - P0);
        result.vertices.reserve(poly.vertices.size() + 1);

        for (int i=0; i<poly.vertices.size(); i++) {
            const Vector &A = (i==0)? poly.vertices[poly.vertices.size() - 1]:poly.vertices[i-1];
            const Vector &B = poly.vertices[i];
            double t = dot(Mprime-A, Pi - P0) / dot(B-A, Pi - P0);

            Vector P = A + t*(B-A);

            if ((B-P0).norm2() - weights[index_0] < (B-Pi).norm2() - weights[index_i]) { //B is inside
                if ((A-P0).norm2() - weights[index_0] > (A-Pi).norm2() - weights[index_i]) { //A outside
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else if ((A-P0).norm2() - weights[index_0] < (A-Pi).norm2() - weights[index_i]) { //A is inside
                result.vertices.push_back(P);
            }
        }
        return result;
    }

    Polygon clip_polygon_by_edge(const Polygon &poly, const Vector& u, const Vector& v) const {
        
        Polygon result;
        result.vertices.reserve(poly.vertices.size() + 1);
        Vector N(v[1]-u[1], -v[0]+u[0], 0);

        for (int i=0; i<poly.vertices.size(); i++) {
            const Vector &A = (i==0)? poly.vertices[poly.vertices.size() - 1]:poly.vertices[i-1];
            const Vector &B = poly.vertices[i];
            double t = dot(u-A, N) / dot(B-A, N);
            Vector P = A + t*(B-A);

            if (dot(B-u, N) < 0) { //B is inside
                if (dot(A-u, N) > 0) { //A outside
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else if (dot(A-u, N) < 0) { //A is inside
                result.vertices.push_back(P);
            }
        }
        return result;
    }

    Polygon intersect_with_disk(const Polygon& polygon, const Vector& center, double radius) {
        Polygon result(polygon);
        for (int i=0; i<disk.vertices.size(); i++) {
            const Vector& u = disk.vertices[i]*radius + center; 
            const Vector& v = disk.vertices[(i+1)%disk.vertices.size()] * radius + center;
            result = clip_polygon_by_edge(result, u, v);
        }
        return result;
    }

    Polygon compute_powerdiagram_cell(int idx) {
        Polygon result;
        result.vertices.resize(4);
        result.vertices[0] = Vector(0, 0, 0);
        result.vertices[1] = Vector(1, 0, 0);
        result.vertices[2] = Vector(1, 1, 0);
        result.vertices[3] = Vector(0, 1, 0);

        for (int i=0; i < points.size(); i++) {
            if (i == idx) {
                continue;
            }
            result = clip_polygon_by_bissector(result, idx, i, points[idx], points[i]);
        }
        result = intersect_with_disk(result, points[idx], sqrt(weights[idx] - weights[weights.size()-1]));

        return result;
    }

    void compute() {
        powerdiagram.resize(points.size());
        for (int i = 0; i < points.size(); i++) {
            powerdiagram[i] = compute_powerdiagram_cell(i);
        }
    }
    void save(std::string filename) {
        save_svg(powerdiagram, filename, "blue");
    }

    std::vector<Polygon> powerdiagram;
    std::vector<Vector> points;
    std::vector<Vector> velocities;
    std::vector<double> weights;
};