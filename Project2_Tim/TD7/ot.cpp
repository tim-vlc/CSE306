#include "powerdiagram.cpp"

class OT {
public:
    OT(std::vector<Vector> &pts, const std::vector<double> &lambdas): pts(pts), lambdas(lambdas) {}

    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<OT*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        lbfgsfloatval_t fx = 0.0;

        for (int i = 0;i < n;i ++) {
            solution.weights[i] = x[i];
        }
        solution.compute();

        double s1 = 0;
        double s2 = 0;
        double s3 = 0;
        
        for (int i = 0;i < n;i ++) {
            double cell_area = solution.powerdiagram[i].area();
            g[i] = -(lambdas[i] - cell_area);
            s1 += solution.powerdiagram[i].integrate_squared_distance(solution.points[i]);
            s2 += lambdas[i] * x[i];
            s3 -= x[i] * cell_area;
        }
        fx = s1 + s2 + s3;


        return -fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<OT*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {   
        for (int i = 0;i < n;i ++) {
            solution.weights[i] = x[i];
        }
        solution.compute();

        double max_diff = 0;
        for (int i = 0;i < n;i ++) {
            double current_area = solution.powerdiagram[i].area();
            double desired_area = lambdas[i];
            max_diff = std::max(max_diff, std::abs(current_area - desired_area));
        }

        std::cout << "fx: " << fx << " max_diff= " << max_diff << "\t gnorm= " << gnorm << std::endl;

        return 0;
    }

    void solve() {
        solution.points = pts;
        solution.weights.resize(pts.size());
        std::fill(solution.weights.begin(), solution.weights.end(), 1.0);


        double fx = 0;
        // LBFGS solve
        int ret = lbfgs(pts.size(), &solution.weights[0], &fx, _evaluate, _progress, this, NULL);

        solution.compute();
    }

    PowerDiagram solution;
    std::vector<Vector> pts;
    std::vector<double> lambdas;
};