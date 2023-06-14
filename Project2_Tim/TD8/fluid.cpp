#include "ot.cpp"

class Fluid {
public:
    Fluid(int N) {
        particles.resize(N);
        for (int i=0; i<N; i++) {
            particles[i] = Vector(rand()/(double)RAND_MAX, rand()/(double)RAND_MAX, rand()/(double)RAND_MAX);
        }
        velocities.resize(N, Vector(0, 0, 0));
    }
    void stepFluid() {
        otsolver.pts = particles;
        otsolver.lambdas = std::vector<double>(particles.size(), 1/particles.size() * VOLUME_FLUID);
        otsolver.solve();
        double mass_particle = 200;
        double epsilon2 = 0.004*0.004;
        const double dt = 0.002;

        for (int i=0; i<particles.size(); i++) {
            Vector gravity = Vector(0, -9.81, 0) * mass_particle;
            Vector centroid = otsolver.solution.powerdiagram[i].centroid();
            Vector otForce = 1. / epsilon2 * (centroid - particles[i]);
            Vector forces = gravity + otForce;
            velocities[i] += dt / mass_particle * forces;
            particles[i] += dt * velocities[i];
        }
    }
    void runFluid() {
        for (int i=0; i < 1000; i++) {
            stepFluid();
            save_frame(otsolver.solution.powerdiagram, "animation", i);
        }
    }
    OT otsolver = OT();
    std::vector<Vector> points;
    std::vector<Vector> velocities;
    std::vector<Vector> particles;
    std::vector<double> weights;
};