#include <iostream>
#include "LinkedList.h"
#include <string>
#include <random>

using namespace std;

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
// TODO: Modify this so that if one list doesn't find it, then it has to not exist so just immediately return false?
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
// For complete cycles and graphs, the E and DIST parameters are not needed.
LinkedList* generateGraph(LinkedList* adjList, const int V, const int E, const string G, const string DIST) {

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
    else if (G == "RANDOM") { // Randomly generated graph
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
        else if (DIST == "YOURS") {
            // find simple distribution
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
    LinkedList adjList[V] = { };

    // Generate the graph specified by the command line arguments.
    generateGraph(adjList, V, E, G, DIST);

    // Display the generated graph.
    printAdjList(adjList, V);

    return 0;
}
