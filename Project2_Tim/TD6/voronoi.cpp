#include "vector.cpp"

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
    std::vector<Vector> vertices;
};  

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
    FILE* f = fopen(filename.c_str(), "w+"); 
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i=0; i<polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \""); 
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

class VoronoiDiagram {
public:
    VoronoiDiagram(std::vector<Vector>& pts) {
        points = pts;
    }

    //Sutherland Hodgman Algorithm
    Polygon clip_polygon_by_bissector(const Polygon &poly, const Vector& P0, const Vector& Pi) {
        Polygon result;
        for (int i=0; i<poly.vertices.size(); i++) {
            const Vector &A = (i==0)? poly.vertices[poly.vertices.size() - 1]:poly.vertices[i-1];
            const Vector &B = poly.vertices[i];
            Vector M = (P0 + Pi) / 2;
            double t = dot(M-A, Pi - P0) / dot(B-A, Pi - P0);

            Vector P = A + t*(B-A);

            if ((B-P0).norm2() < (B-Pi).norm2()) { //B is inside
                if ((A-P0).norm2() > (A-Pi).norm2()) { //A outside
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else if ((A-P0).norm2() < (A-Pi).norm2()) { //A is inside
                result.vertices.push_back(P);
            }
        }
        return result;
    }

    Polygon compute_voronoi_cell(int idx) {
        Polygon result;
        result.vertices.resize(4);
        result.vertices[0] = Vector(0, 0, 0);
        result.vertices[1] = Vector(0, 1, 0);
        result.vertices[2] = Vector(1, 1, 0);
        result.vertices[3] = Vector(1, 0, 0);

        for (int i=0; i < points.size(); i++) {
            if (i == idx) {
                continue;
            }
            result = clip_polygon_by_bissector(result, points[idx], points[i]);
        }

        return result;
    }

    void compute() {
        voronoi.resize(points.size());
        for (int i =0; i < points.size(); i++) {
            voronoi[i] = compute_voronoi_cell(i);
        }
    }
    void save(std::string filename) {
        save_svg(voronoi, filename, "white");
    }
    std::vector<Polygon> voronoi;
    std::vector<Vector> points;
};