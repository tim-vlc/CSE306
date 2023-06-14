#include "polygon.cpp"

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

class PowerDiagram {
public:
    PowerDiagram() {}
    PowerDiagram(std::vector<Vector>& pts, const std::vector<double>& weights) {
        points = pts;
        this->weights = weights;
    }

    //Sutherland Hodgman Algorithm
    Polygon clip_polygon_by_bissector(const Polygon &poly, int index_0, int index_i, const Vector& P0, const Vector& Pi) {
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

        return result;
    }

    void compute() {
        powerdiagram.resize(points.size());
        for (int i =0; i < points.size(); i++) {
            powerdiagram[i] = compute_powerdiagram_cell(i);
        }
    }
    void save(std::string filename) {
        save_svg(powerdiagram, filename, "white");
    }
    std::vector<Polygon> powerdiagram;
    std::vector<Vector> points;
    std::vector<double> weights;
};