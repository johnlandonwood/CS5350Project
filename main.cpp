#include <iostream>
#include "LinkedList.h"
#include <string>
#include <random>

#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <cmath>

// https://smu.hosted.panopto.com/Panopto/Pages/Viewer.aspx?id=170f1626-3bb5-4363-966f-afce002b29bb
// TODO: Create doubly linked list, first of ints, then of structs
// Change "yours" distribution to triangle distribution
// Option for final ordering: largest last vertex. should perform worse than smallest last.

using namespace std;

// Vertex struct to hold information for each vertex in the graph.
struct Vertex {
    int id;
    LinkedList edges;
    int degree;
    int deleted;
    int color;
    // DLL pointer
    // order deleted pointer
};

// Method to print the adjacency list.
void printAdjList(Vertex vertices[], int V) {
    for (int i = 0; i < V; i++) {
        cout << i << " --> ";
        vertices[i].edges.print();
        cout << " | Degree: " << vertices[i].degree << endl;
    }
}

// Method to check if an edge exists between two vertices.
// Returns true if the edge exists, or false if not.
// Time complexity can get bad here for large graphs - this function is O(2*n) I think
bool edgeExists(Vertex v1, Vertex v2) {
    bool v2_in_v1 = false;
    bool v1_in_v2 = false;

    Node *temp1 = v1.edges.head;
    Node *temp2 = v2.edges.head;

    // Search for v2 in v1's list
    while (temp1 != nullptr) {
        if (temp1->data == v2.id) {
            v2_in_v1 = true;
            break;
        }
        temp1 = temp1->next;
    }

    // Could get rid of this check
    // Since every time we are inserting, we have to insert it in both lists anyways
    // But keeping it in should improve performance on large graphs
//    if (temp1 == nullptr) {
//        return false;
//    }

    // Search for v1 in v2's list
    while (temp2 != nullptr) {
        if (temp2->data == v1.id) {
            v1_in_v2 = true;
            break;
        }
        temp2 = temp2->next;
    }

    return (v2_in_v1 && v1_in_v2);
}




//else if (DIST == "YOURS") { // Triangle distribution
//    // Can create a triangle distribution by adding two uniform values
//    int v1, v2;
//    std::random_device rd;
//    std::mt19937 rng(rd());
//    std::uniform_int_distribution<int> uniform_distribution(0,V-1);
//    int tri1, tri2, tri3, tri4;
//
//    tri1 = uniform_distribution(rng);
//    tri2 = uniform_distribution(rng);
//    v1 = tri1 + tri2;
//    tri3 = uniform_distribution(rng);
//    tri4 = uniform_distribution(rng);
//    v2 = tri3 + tri4;
//
//    std::map<int, int> hist{};
//    for (int n = 0; n != 10000; ++n) {
//        tri1 = uniform_distribution(rng);
//        tri2 = uniform_distribution(rng);
//        v1 = tri1 + tri2;
//        tri3 = uniform_distribution(rng);
//        tri4 = uniform_distribution(rng);
//        v2 = tri3 + tri4;
//        //                if (v1 <= V-1 && v2 <= V-1) {
//        ++hist[v1];
//        ++hist[v2];
//        //                }
//
//    }
//
//
//    for (auto [x, y] : hist)
//        std::cout << std::setw(2) << x << ' ' << std::string(y / 200, '*') << '\n';
//
//    //            while (edges_added_ctr < E) { // Until we have added E edges to the graph:
//    //                // Generate two pairs of random uniform numbers
//    //                // The vertex we select is the sum of the two random uniform numbers
//    //                tri1 = uniform_distribution(rng);
//    //                tri2 = uniform_distribution(rng);
//    //                v1 = tri1 + tri2;
//    //                tri3 = uniform_distribution(rng);
//    //                tri4 = uniform_distribution(rng);
//    //                v2 = tri3 + tri4;
//    //
//    //
//    //
//    //                // If the vertices generated are not out of bounds (since summing up values can exceed V-1):
//    //                if (v1 <= V-1 && v2 <= V-1) {
//    //                    // If the edge is not self-referential and does not already exist, add it to the adjacency list.
//    //                    if (v1 != v2 && !edgeExists(vertices[v1], vertices[v2])) {
//    //                        vertices[v1].edges.insert(v2);
//    //                        vertices[v1].degree++;
//    //                        vertices[v2].edges.insert(v1);
//    //                        vertices[v2].degree++;
//    //                        edges_added_ctr++;
//    //                    }
//    //                }
//    //            }
//}
//cout << "Edges added: " << edges_added_ctr << ", expected: " << E << endl;
//}


// TODO: Only doing the vertices array for now, need to add DLL support.
// Method to generate a graph with the specified command line arguments
// For complete cycles and graphs, the E and DIST parameters are unnecessary.
Vertex* generateGraph(Vertex vertices[], const int V, const int E, const string& G, const string& DIST) {
    int edges_added_ctr = 0;
    if (G == "COMPLETE") { // Generate a complete graph: |V| = V, |E| = (V * (V - 1))/2;
        int i, j;
        for (i = 0; i < V - 1; i++) {
            for (j = i + 1; j < V; j++) {
                if (i != j) {
                    vertices[i].edges.insert(j); // Add the edge to each vertex's adjacency list
                    vertices[i].degree++;
                    vertices[j].edges.insert(i);
                    vertices[j].degree++;
                    edges_added_ctr++;
                }
            }
        }
        cout << "Edges added: " << edges_added_ctr << ", expected: " << (V * (V - 1))/2 << endl;
    }
    else if (G == "CYCLE") { // Generate a cycle: |V| = V, |E| = |V|
        int i = 0;
        int j = 1;
        while (i < E) {
            if (j == E) {
                i = j - 1;
                j = 0;
                vertices[i].edges.insert(j); // Add the edge to each vertex's adjacency list
                vertices[i].degree++;
                vertices[j].edges.insert(i);
                vertices[j].degree++;
                edges_added_ctr++;
                break;
            }
            vertices[i].edges.insert(j); // Add the edge to each vertex's adjacency list
            vertices[i].degree++;
            vertices[j].edges.insert(i);
            vertices[j].degree++;
            edges_added_ctr++;
            i++;
            j++;
        }
        cout << "Edges added: " << edges_added_ctr << ", expected: " << V << endl;
    }
    else if (G == "RANDOM") { // Generate a random graph
        if (DIST == "UNIFORM") { // Uniform distribution
            std::random_device rd;
            std::mt19937 rng(rd());
            std::uniform_int_distribution<int> uniform_distribution(0,V-1);
            int v1, v2;
            while (edges_added_ctr < E) { // Until we have added E edges to the graph:
                v1 = uniform_distribution(rng); // Pick two random vertices, v1 and v2
                v2 = uniform_distribution(rng);
                // If the edge is not self-referential and does not already exist, add it to the adjacency list.
                if (v1 != v2 && !edgeExists(vertices[v1], vertices[v2])) {
                    vertices[v1].edges.insert(v2);
                    vertices[v1].degree++;
                    vertices[v2].edges.insert(v1);
                    vertices[v2].degree++;
                    edges_added_ctr++;
                }
            }
        }
        else if (DIST == "SKEWED") { // Skewed distribution
            // Generate skewed distribution array and select two indices from it for the edge
            // This has the same effect as a skewed distribution.
            int v1, v2;
            int size = (V * (V+1)) / 2;
            int skewed[size];
            int x, i;
            int start = 0;
            int end = V - 1;
            int num_times = V-1;
            // Generate skewed distribution array
            for (x = 0; x < V; x++) {
                for (i = start; i <= end; i++) {
                    skewed[i] = x;
                }
                num_times--;
                start = end + 1;
                end = start + num_times;
            }
            std::random_device rd;
            std::mt19937 rng(rd());
            std::uniform_int_distribution<int> uniform_distribution(0,size-1);
            int index1, index2;
            while (edges_added_ctr < E) { // Until we have added E edges to the graph:
                index1 = uniform_distribution(rng); // Generate two random indices and use the values there
                index2 = uniform_distribution(rng);  // for the vertices
                v1 = skewed[index1];
                v2 = skewed[index2];
                // If the edge is not self-referential and does not already exist, add it to the adjacency list.
                if (v1 != v2 && !edgeExists(vertices[v1], vertices[v2])) {
                    vertices[v1].edges.insert(v2);
                    vertices[v1].degree++;
                    vertices[v2].edges.insert(v1);
                    vertices[v2].degree++;
                    edges_added_ctr++;
                }
            }
        }
        else if (DIST == "YOURS") { // Normal distribution
            // https://www.baeldung.com/cs/uniform-to-normal-distribution
            // https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
            // Box-muller transform on uniform distribution

        }
    }

    return vertices;
}

int main(int argc, char** argv) {

    const int V = atoi(argv[1]); // MAX = 10,000
    const int E = atoi(argv[2]); // MAX = 2,000,000
    const string G = argv[3]; // COMPLETE | CYCLE | RANDOM (with DIST below)
    const string DIST = argv[4]; // UNIFORM | SKEWED | YOURS

    cout << V << " " << E << " " << G << " " << DIST << endl;

    // Array of structs to hold vertex information
    Vertex vertices[V];

    // Allocate and initialize V vertices
    for (int i = 0; i < V; i++) {
        Vertex v;
        v.id = i;
        v.edges.head = nullptr;
        v.degree = 0;
        v.color = -1;
        v.deleted = 0;
        // v.dll-> head = nullptr
        // v.order_deleted = nullptr
        vertices[i] = v;
    }

    // Generate the graph specified by the command line arguments.
    generateGraph(vertices, V, E, G, DIST);

    // Display the adjacency list describing the graph.
    printAdjList(vertices, V);


    cout << '\n' << "Done" << endl;
    return 0;
}
