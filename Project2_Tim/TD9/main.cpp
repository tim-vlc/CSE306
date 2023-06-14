#include "edge.cpp"

int main() {
    TriangleMesh mesh;
    mesh.readOBJ("goethe.obj");

    mesh.center_scale();

    std::map< Edge, std::vector<int> > edge_to_triangle;
    std::map< int, std::vector<int> > vertex_to_triangle;

    for (int i=0; i<mesh.indices.size(); i++) {
        TriangleIndices triangle = mesh.indices[i];
        Edge e1(triangle.vtxi, triangle.vtxj);
        Edge e2(triangle.vtxj, triangle.vtxk);
        Edge e3(triangle.vtxk, triangle.vtxi);

        edge_to_triangle[e1].push_back(i);
        edge_to_triangle[e2].push_back(i);
        edge_to_triangle[e3].push_back(i);

        vertex_to_triangle[triangle.vtxi].push_back(i);
        vertex_to_triangle[triangle.vtxj].push_back(i);
        vertex_to_triangle[triangle.vtxk].push_back(i);
    }

    std::vector<bool> is_vertex_on_boundary(mesh.vertices.size(), false);
    
    std::vector<Edge> boundary_edges;
    for (auto it = edge_to_triangle.begin(); it != edge_to_triangle.end(); it++) {
        if (it->second.size() == 1) {
            boundary_edges.push_back(it->first);
            is_vertex_on_boundary[it->first.a] = true;
            is_vertex_on_boundary[it->first.b] = true;
        }
    }

    std::vector<Edge> ordered_boundary_edges(boundary_edges.size());
    int current_edge = 0;
    ordered_boundary_edges[0] = boundary_edges[current_edge];

    for (int i=1; i<boundary_edges.size(); i++) {
        for (int j=0; j<boundary_edges.size(); j++) {
            if (boundary_edges[j].a == ordered_boundary_edges[i-1].b) {
                ordered_boundary_edges[i] = boundary_edges[j];
                break;
            }
        }
    }
    for (int i=0; i<ordered_boundary_edges.size(); i++) {
        double theta = 2*M_PI*i/ (double)ordered_boundary_edges.size();
        Vector circle_vtx;
        circle_vtx[0] = cos(theta) * 0.5 + 0.5;
        circle_vtx[1] = sin(theta) * 0.5 + 0.5;
        circle_vtx[2] = 0;

        mesh.vertices[ordered_boundary_edges[i].a] = circle_vtx;
    }

    for (int iter = 0; iter < 1000; iter++) {
        std::vector<Vector> updated_vertices(mesh.vertices.size());
        updated_vertices = mesh.vertices;

        for (int i=0; i<mesh.vertices.size(); i++) {
            if (is_vertex_on_boundary[i]) {
                continue;
            }

            Vector avg_neighbors(0, 0, 0);
            int nb_neighbors = 0;

            for (int j=0; j<vertex_to_triangle[i].size(); j++) { // Loop over the neighboring triangles

                TriangleIndices triangle = mesh.indices[vertex_to_triangle[i][j]];
                Vector A = mesh.vertices[triangle.vtxi];
                Vector B = mesh.vertices[triangle.vtxj];
                Vector C = mesh.vertices[triangle.vtxk];

                if (triangle.vtxi != i) {
                    avg_neighbors += A;
                }
                if (triangle.vtxj != i) {
                    avg_neighbors += B;
                }
                if (triangle.vtxk != i) {
                    avg_neighbors += C;
                }
                nb_neighbors += 2;

            }
            updated_vertices[i] = avg_neighbors/nb_neighbors;
        }
        mesh.vertices = updated_vertices;
    }
    mesh.writeOBJ("flattened.obj");
}