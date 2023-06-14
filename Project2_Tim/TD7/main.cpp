#include "ot.cpp"

int main() {
    std::vector<Vector> points(256);
    std::vector<double> lambdas(256);

    for (int i = 0; i < points.size(); i++) {
        points[i][0] = rand() / (double)RAND_MAX;
        points[i][1] = rand() / (double)RAND_MAX;
        points[i][2] = 0;
        lambdas[i] = 1. / points.size();
    }

    OT ot(points, lambdas);
    auto start = std::chrono::high_resolution_clock::now();
    ot.solve();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);

    std::cout << "Time: " << duration.count() << std::endl;

    ot.solution.save("ot.svg");
}