#include <iostream>
#include "LinkedList.h"
#include <string>
#include <random>

// https://smu.hosted.panopto.com/Panopto/Pages/Viewer.aspx?id=170f1626-3bb5-4363-966f-afce002b29bb
// TODO: Create struct for vertices with 7 fields
// TODO: Create array of structs as described
// TODO: Create doubly linked list, first of ints, then of structs
// TODO: Redo normal distribution and implement skewed distribution
// Can redo normal distribution as the community package one?

using namespace std;

struct Vertex {
    int id;
    LinkedList *edges;
    int degree;
    int deleted;
    int color;
    // DLL pointer
    // order deleted pointer
};

// Method to print the adjacency list.
void printAdjList(LinkedList *adjList, int V) {
    for (int i = 0; i < V; i++) {
        cout << i << " --> ";
        adjList[i].print();
    }
}

// Return true if an edge between two vertices exists.
// Otherwise, return false.
// Time complexity can get bad here for large graphs - this function is O(2n) I think
// Since every time we are inserting, we have to insert it in both lists anyways
// Note: this will segfault if v2 >= V. If v1 >= V, it just returns false, idk why that is
bool edgeExists(LinkedList* adjList, int v1, int v2) {
    bool v2_in_v1 = false;
    bool v1_in_v2 = false;

    Node *temp1 = adjList[v1].head;
    Node *temp2 = adjList[v2].head;

    // Search for v2 in v1's list
    while (temp1 != nullptr) {
        if (temp1->data == v2) {
            v2_in_v1 = true;
            break;
        }
        temp1 = temp1->next;
    }
    if (temp1 == nullptr) {
        return false;
    }

    // Search for v1 in v2's list
    while (temp2 != nullptr) {
        if (temp2->data == v1) {
            v1_in_v2 = true;
            break;
        }
        temp2 = temp2->next;
    }

    return (v2_in_v1 && v1_in_v2);
}



// Method to generate a graph with the specified command line arguments
// For complete cycles and graphs, the E and DIST parameters are unnecessary.
LinkedList* generateGraph(LinkedList* adjList, const int V, const int E, const string& G, const string& DIST) {

    int edges_added_ctr = 0;
    if (G == "COMPLETE") { // Generate a complete graph: |V| = V, |E| = (V * (V - 1))/2;
        int i, j;
        for (i = 0; i < V - 1; i++) {
            for (j = i + 1; j < V; j++) {
                if (i != j) {
                    adjList[i].insert(j); // Add the edge to each vertex's adjacency list
                    adjList[j].insert(i);
                    edges_added_ctr++;
                }
            }
        }
        cout << "Edges added: " << edges_added_ctr << ", expected: " << (V * (V - 1))/2 << endl;
    }
    else if (G == "CYCLE") { // Generate a cycle: |V| = V, |E| = E
        int i = 0;
        int j = 1;
        while (i < E) {
            if (j == E) {
                i = j - 1;
                j = 0;
                adjList[i].insert(j);
                adjList[j].insert(i);
                edges_added_ctr++;
                break;
            }
            adjList[i].insert(j);
            adjList[j].insert(i);
            edges_added_ctr++;
            i++;
            j++;
        }
        cout << "Edges added: " << edges_added_ctr << ", expected: " << E << endl;
    }
    else if (G == "RANDOM") { // Generate a random graph
        if (DIST == "UNIFORM") { // Uniform distribution
            std::random_device rd;
            std::mt19937 rng(rd());
            std::uniform_int_distribution<int> uniform_distribution(0,V-1);
            int v1, v2;
            while (edges_added_ctr < E) { // Until we have added E edges to the graph:
                v1 = uniform_distribution(rng); // Generate two random vertices, v1 and v2
                v2 = uniform_distribution(rng);
                // If the edge is not self-referential and does not already exist, add it to the adjacency list.
                if (v1 != v2 && !edgeExists(adjList, v1, v2)) {
                    adjList[v1].insert(v2);
                    adjList[v2].insert(v1);
                    edges_added_ctr++;
                }
            }
        }
        else if (DIST == "SKEWED") { // Skewed distribution
            // use skewed community package
        }
        else if (DIST == "YOURS") { // C++ std::normal distribution
            int v1, v2 = 0;
            std::random_device rd;
            std::mt19937 rng(rd());
            int mid = V / 2;
            cout << "Mean: " << mid << ", stddev: " << mid / 2 << endl;
            // The normal distribution uses the middle value between 0 and V as the mean,
            // and uses the (mean / 2) as the standard deviation.
            std::normal_distribution<float> normal_distribution(mid, mid/2);
            while (edges_added_ctr < E) { // Until we have added E edges to the graph:
                v1 = normal_distribution(rng); // Generate two random vertices, v1 and v2
                v2 = normal_distribution(rng);
                // Normal distributions have a  chance to generate a number outside the intended bound [0, V-1].
                // Adjacency list cannot allow numbers outside those bounds to be added, so if we generate a vertex
                // outside those bounds, we set it to something else.
                if (v1 < 0 || v1 > (V-1)) { v1 = 0; } // Arbitrarily set v1 to 0.
                if (v2 < 0 || v2 > (V-1)) { v2 = V - 1; } // Arbitrarily set v2 to V - 1.
                // If the edge is not self-referential and does not already exist, add it to the adjacency list.
                if (v1 != v2 && !edgeExists(adjList, v1, v2)) {
                    adjList[v1].insert(v2);
                    adjList[v2].insert(v1);
                    edges_added_ctr++;
                }
            }
        }
        cout << "Edges added: " << edges_added_ctr << ", expected: " << E << endl;
    }

    return adjList;
}


int main(int argc, char** argv) {

    const int V = atoi(argv[1]); // MAX = 10,000
    const int E = atoi(argv[2]); // MAX = 2,000,000
    const string G = argv[3]; // COMPLETE | CYCLE | RANDOM (with DIST below)
    const string DIST = argv[4]; // UNIFORM | SKEWED | YOURS

    cout << V << " " << E << " " << G << " " << DIST << endl;

    // Adjacency list - array of singly linked lists.
//    LinkedList adjList[V] = {};
//
//    // Generate the graph specified by the command line arguments.
//    generateGraph(adjList, V, E, G, DIST);
//
//    // Display the generated graph.
//    printAdjList(adjList, V);

    // generate V vertex structures

//    Vertex vertices[V];
//
//    for (int i = 0; i < V; i++) {
//        Vertex v; // Initialize
//        v.id = i;
//        //v.edges->head = nullptr;
//        v.degree = 0;
//        v.color = -1;
//        v.deleted = 0;
//        // v.dll-> head = nullptr
//        // v.order_deleted = nullptr
//        vertices[i] = v;
//    }
//
//    // From panopto example:
//    // 0 --> 1
//    // 1 --> 0 --> 2 --> 3
//    // 2 --> 1 --> 3
//    // 3 --> 1 --> 2
//    vertices[0].edges->insert(1);
//    vertices[1].edges->insert(0);
//    vertices[1].edges->insert(2);
//    vertices[1].edges->insert(3);
//    vertices[2].edges->insert(1);
//    vertices[2].edges->insert(3);
//    vertices[3].edges->insert(1);
//    vertices[3].edges->insert(2);
//
//    for (int i = 0; i < V; i++) {
//        vertices[i].edges->print();
//    }

//    for (int i = 0; i < V; i++) {
//        cout << vertices[i].id << " ";
//    }

    return 0;
}
