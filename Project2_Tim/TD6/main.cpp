#include "voronoi.cpp"

int main() {
    std::vector<Vector> points(256);

    for (int i = 0; i < points.size(); i++) {
        points[i][0] = rand() / (double)RAND_MAX;
        points[i][1] = rand() / (double)RAND_MAX;
        points[i][2] = 0;
    }

    VoronoiDiagram voronoi(points);
    voronoi.compute();
    voronoi.save("voronoi.svg");
}