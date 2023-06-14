#include "triangle.cpp"

class Edge {
public:
    Edge(int a=0, int b=0) {
        if (a < b) {
            this->a = a;
            this->b = b;
        } else {
            this->a = b;
            this->b = a;
        }
    }


    int a, b;
};

bool operator<(const Edge& e1, const Edge& e2)  {
    int mine1 = std::min(e1.a, e1.b);
    int maxe1 = std::max(e1.a, e1.b);
    int mine2 = std::min(e2.a, e2.b);
    int maxe2 = std::max(e2.a, e2.b);
    return std::pair<int, int>(mine1, maxe1) < std::pair<int, int>(mine2, maxe2);
}