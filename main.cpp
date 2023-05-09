#include <iostream>
#include <string>
#include <random>
#include <climits>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <cmath>
#include <limits>
#include <utility>
#include <cmath>
#include <map>
#include "LinkedList.h"
#include "DoublyLinkedList.h"
using namespace std;
// https://smu.hosted.panopto.com/Panopto/Pages/Viewer.aspx?id=170f1626-3bb5-4363-966f-afce002b29bb


// Vertex struct to hold information for each vertex in the graph.
struct Vertex {
    int id; // Unique vertex ID
    LinkedList edges; // Edge list for connected vertices
    int degree; // Current degree of the vertex
    bool deleted; // True when vertex is deleted from graph, false otherwise
    int color; // Color value
    int degree_when_deleted; // Degree when deleted from graph (only applicable in smallest last vertex ordering)
    int original_degree; // Original degree of the vertex upon creation
    DoublyLinkedList* degree_DLL; // Pointer to DLL of vertices of same current degree
    LinkedList* ordering; // Pointer to ordering list
};

// Global array to hold, in order:
// total colors used, average original degree, maximum degree when deleted, and terminal clique size
double stats[4] = { 0.0, 0.0, 0.0, 0.0 };

// Global array to hold graph generation runtime and smallest last vertex ordering runtime if applicable
long runtimes[2] = {0, 0};

long runtime_data[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };


// Method to print the adjacency list.
void printAdjList(Vertex vertices[], int V) {
    for (int i = 0; i < V; i++) {
        cout << i << " --> ";
        vertices[i].edges.print();
        cout << " | degree: " << vertices[i].degree;
        cout << " | original_degree: " << vertices[i].original_degree;
        cout << endl;
    }
    cout << endl;
}

// Method to print the degree-indexed doubly linked lists.
void printDegreeDLLs(DoublyLinkedList degree_DLLs[], int V) {
    for (int i = 0; i < V; i++) {
        cout << i << ": " << &degree_DLLs[i] << " | ";
        degree_DLLs[i].print();
        cout << endl;
    }
    cout << endl;
}

// Method to write the final information of the graph and its coloring
void recordOrderingAndColoring(ofstream& output, Vertex vertices[], LinkedList ordering, int coloring_order[], int V, const string& ORDERING) {

    output << "Vertex, Color, Original degree";
    if (ORDERING == "SMALLEST_LAST")  {
        output << ", Degree when deleted, Total colors used, Avg. original degree, Max degree when deleted, Terminal clique size, Ordering runtime";
    }
    else {
        output << ", Total colors used, Avg. original degree, Graph generation";
    }
    output << '\n';
    for (int i = 0; i < V; i++) {
        output << vertices[i].id << ", " << vertices[i].color << ", " << vertices[i].original_degree;
        if (ORDERING == "SMALLEST_LAST") {
            output << ", " << vertices[i].degree_when_deleted;
        }
        if (i == 0) {
            output << setprecision(1) << ", " << stats[0] << ", ";
            output << setprecision(3) <<  stats[1];
            if (ORDERING == "SMALLEST_LAST") {
                output  << setprecision(1) << ", " << stats[2] << ", ";
                output << setprecision(1) << stats[3] << ", ";
                output << runtimes[1];
            }
        }
        output << '\n';
    }
}

// Method to check if an edge exists between two vertices.
// Returns true if the edge exists, or false if not.
// Time complexity can get bad here for large graphs - this function is O(2*n) I think
bool edgeExists(Vertex* v1, Vertex* v2) {
    bool v2_in_v1 = false;
    bool v1_in_v2 = false;

    Node *temp1 = v1->edges.head;
    Node *temp2 = v2->edges.head;

    // Search for v2 in v1's list
    while (temp1 != nullptr) {
        if (temp1->data == v2->id) {
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
        if (temp2->data == v1->id) {
            v1_in_v2 = true;
            break;
        }
        temp2 = temp2->next;
    }

    return (v2_in_v1 && v1_in_v2);
}

// Method to generate two random numbers following a normal distribution.
// Uses the Box-Muller transform.
// Implementation modified from https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform#Implementation.
std::pair<double, double> box_muller_transform(double mu, double sigma, int V) {
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    constexpr double two_pi = 2.0 * M_PI;
    double V_double = static_cast<double>(V);

    // Initialize the random uniform number generator in a range 0 to V-1
    static std::mt19937 rng(std::random_device{}());
    static std::uniform_int_distribution<int> runif(0, V-1);

    // Create two random numbers, ensuring u1 is greater than epsilon
    // The Box-Muller transform needs decimal numbers between 0 and 1 to work,
    // so we divide the random integers generated by V
    double u1, u2;
    // is this being a problem?
    while (u1 <= epsilon) {
        u1 = runif(rng)/V_double;
    }
    u2 = runif(rng)/V_double;

    // Compute z0 and z1
    auto mag = sigma * sqrt(-2.0 * log(u1));
    auto z0  = mag * cos(two_pi * u2) + mu;
    auto z1  = mag * sin(two_pi * u2) + mu;

    // Return the pair of random doubles, to be converted to ints for edge insertion
    return std::make_pair(z0, z1);
}

// Method to calculate mean of a population from [0, V-1].
double calculate_mean(int V) {
    double sum = 0.0;
    for (int i = 0; i < V; i++) {
        sum += static_cast<double>(i);
    }
    sum = sum / static_cast<double>(V);
    return sum;
}

// Method to calculate standard deviation of a population from [0, V-1].
double calculate_stddev(double mean, int V) {
    double variance = 0.0;
    for (int i = 0; i < V; i++) {
        variance += pow((i - mean), 2);
    }
    double stddev = sqrt(variance / static_cast<double>(V));
    return stddev;
}

// Method to generate a graph with the specified command line arguments
// For complete cycles and graphs, the E and DIST parameters are unnecessary.
void generateGraph(Vertex vertices[], DoublyLinkedList degree_DLLs[], const int V, const int E, const string& G, const string& DIST) {
    std::map<int, int> histogram {}; // For recording vertex chosen histograms
    auto time_start = std::chrono::high_resolution_clock::now(); // For runtime analysis data

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
        // cout << "Edges added: " << edges_added_ctr << ", expected: " << (V * (V - 1))/2 << endl;
    }
    else if (G == "CYCLE") { // Generate a cycle: |V| = V, |E| = |V|
        int i = 0;
        int j = 1;
        while (i < V) {
            if (j == V) {
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
        // cout << "Edges added: " << edges_added_ctr << ", expected: " << V << endl;
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
                if (v1 != v2 && !edgeExists(&vertices[v1], &vertices[v2])) {
                    vertices[v1].edges.insert(v2);
                    vertices[v1].degree++;
                    vertices[v2].edges.insert(v1);
                    vertices[v2].degree++;
                    edges_added_ctr++;
                    ++histogram[(v1)]; // Update vertex selection histogram
                    ++histogram[(v2)];
                }
            }
        }
        else if (DIST == "SKEWED") { // Skewed distribution
            // Generate skewed distribution array and select two indices from it for the edge
            // This has the same effect as a skewed distribution.
            int v1, v2;
            int size = (V * (V+1)) / 2;
            int* skewed = new int[size];
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
                if (v1 != v2 && !edgeExists(&vertices[v1], &vertices[v2])) {
                    vertices[v1].edges.insert(v2);
                    vertices[v1].degree++;
                    vertices[v2].edges.insert(v1);
                    vertices[v2].degree++;
                    edges_added_ctr++;
                    ++histogram[(v1)]; // Update vertex selection histogram
                    ++histogram[(v2)];
                }
            }
            delete[] skewed;
        }
        else if (DIST == "YOURS") { // Normal distribution
            // Uses the Box-Muller transform to turn a uniform int distribution into a normal distribution.
            // Calculate mean and standard deviation to be used in Box-Muller transform
            double mean = calculate_mean(V);
            double stddev = calculate_stddev(mean, V);
            int v1, v2;
            std::pair<double, double> edge_pair;
            while (edges_added_ctr < E) { // Until we have added E edges to the graph:
                // Call function to turn uniform int distribution into normal
                // We must re-call the function every time as we must re-seed the random number generator.
                edge_pair = box_muller_transform(mean, stddev, V);

                // Convert edge_pair doubles to ints
                v1 = floor(edge_pair.first);
                v2 = floor(edge_pair.second);

                // The normal distribution can generate numbers outside the intended bounds of [0, V-1],
                // so we must check that v1 and v2 are not outside the bounds as we cannot add them to the graph if so.
                // alternative: set arbitrarily to 0 and v-1
                if (((v1 < 0 || v1 > V-1) || (v2 < 0 || v2 > V-1))) {
                    continue;
                }
                else {
                    if (v1 != v2 && !edgeExists(&vertices[v1], &vertices[v2])) {
                          //cout << "Adding " << v1 << " -- > " << v2 << endl;
                        vertices[v1].edges.insert(v2);
                        vertices[v1].degree++;
                        vertices[v2].edges.insert(v1);
                        vertices[v2].degree++;
                        edges_added_ctr++;
                        ++histogram[(v1)]; // Update vertex selection histogram
                        ++histogram[(v2)];
                    }
                }
//                if (!(v1 < V && v1 >= 0) || !(v2 < V && v2 >= 0)) { // could change to or and continue?
//                    // If the edge is not self-referential and does not already exist, add it to the adjacency list.
////                    bool yo = edgeExists(vertices[v1], vertices[v2]);
////                    cout << "edgeExists(" << v1  << "," << v2 <<"): " << yo << endl;
////                    if (!yo) {
////                        cout << v1 << " --> " << v2 << " does not exist" << endl;
////                    }
//                    if (v1 != v2 && !edgeExists(&vertices[v1], &vertices[v2])) {
//                       //  cout << "Adding " << v1 << " -- > " << v2 << endl;
//                        vertices[v1].edges.insert(v2);
//                        vertices[v1].degree++;
//                        vertices[v2].edges.insert(v1);
//                        vertices[v2].degree++;
//                        edges_added_ctr++;
//                        ++histogram[(v1)]; // Update vertex selection histogram
//                        ++histogram[(v2)];
//                    }
////                    else {
////                        cout << v1 << " --> " << v2 << " are equal or edge already exists" << endl;
////                    }
//                }
//                else {
//                    cout << "Couldn't add " << v1 << " or " << v2 << " (bounds: [0, " << V-1 << "])" << endl;
//                }
            }
        }
//        cout << "Edges added: " << edges_added_ctr << ", expected: " << E << endl;
    }

    // Finish timing and record
    auto time_stop = std::chrono::high_resolution_clock::now();
    auto time_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(time_stop - time_start).count();
    // cout << "Graph generation time_elapsed: " << time_elapsed << endl;
    //cout << V << ", " <<  E << ", " << time_elapsed << endl;
    cout << time_elapsed << endl;
    runtimes[0] = time_elapsed;

    // Record an output .csv with the histogram of vertices chosen for conflict
    ofstream hist("histogram.csv");

    hist << "Vertex, Times chosen, Generation runtime" << '\n';
    for (int i = 0; i < V; i++) {
        if (histogram.count(i) != 0) {
            hist  << i << ", " << histogram[i];
            if (i == 0) {
                hist << ", " << runtimes[0];
            }
            hist << '\n';
        }
        else {
            hist  << i << ", " << 0 << '\n';
        }
    }
    hist.close();

    // Now that the graph is done generating, populate the degree-indeed doubly-linked list
    // with the IDs of each vertex that has the corresponding degree
    for (int i = 0; i < V; i++) {
        Vertex v = vertices[i];
        vertices[i].original_degree = v.degree; // Capture the original degree of the vertex, as degree may change while ordering
        degree_DLLs[v.degree].insert(v.id);
        vertices[v.id].degree_DLL = &degree_DLLs[v.degree];
    }

}

// Method to generate a smallest last vertex ordering.
void smallestLastVertexOrdering(Vertex vertices[], DoublyLinkedList degree_DLLs[], int V) {
    int max_degree_when_deleted = INT_MIN;
    int deleted_ctr = 0;
    int start = 0;
    int degree_just_deleted;
    Vertex v, c;
    auto time_start = std::chrono::high_resolution_clock::now(); // For runtime analysis data
    while (deleted_ctr != V) { // While the graph is not empty:
        cout << "----------------Start: " << start << "----------------"<< endl;
        DLLNode *curr = degree_DLLs[start].head; // Find vertex of smallest degree, v, at start index
        if (curr != nullptr) {
            v = vertices[curr->data];
            cout << "Vertex of smallest degree: " << v.id << endl;
            vertices[v.id].deleted = true; // Mark v as deleted
            vertices[v.id].degree_when_deleted = v.degree;
            if (vertices[v.id].degree_when_deleted > max_degree_when_deleted) { // Update max_degree_when_deleted if necessary
                max_degree_when_deleted = vertices[v.id].degree_when_deleted;
            }
            degree_just_deleted = v.degree;
            vertices[v.id].ordering->insert(v.id); // Add v to the ordering list
            vertices[v.id].degree_DLL->remove(v.id); // Remove v from its degree_DLL list
            cout << "Removing " << v.id << " from degree_DLLs[" << v.degree << "]" << endl;

            // For each vertex c connected to v:
            Node* temp = v.edges.head;
            while (temp != nullptr) {
                c = vertices[temp->data];
                cout << "----Adjusting " << c.id << "----" << endl;
                if (!c.deleted) { // If c has not already been deleted:
                    cout << "Removing " << c.id << " from degree_DLLs[" << c.degree << "]" << endl;
                    vertices[c.id].degree_DLL->remove(c.id); // Remove c from previous degree_DLL list

                    cout << "Degree of " << vertices[c.id].id << " decremented: " << vertices[c.id].degree;
                    vertices[c.id].degree--;  // Update c's degree
                    cout << " -> " << vertices[c.id].degree << endl;

                    cout << "degree_DLL pointer updated: " << vertices[c.id].degree_DLL;
                    vertices[c.id].degree_DLL = &degree_DLLs[vertices[c.id].degree]; // Update c's degree_DLL pointer
                    cout << " -> " << vertices[c.id].degree_DLL << endl;

                    cout << "Adding " << c.id << " to degree_DLL[" << vertices[c.id].degree << "]" << endl;
                    vertices[c.id].degree_DLL->insert(c.id); // Add current vertex to updated degree_DLL list
                    temp = temp->next;
                }
                else { // c has already been deleted, so ignore it
                    cout << c.id << " has already been deleted, skipping" << endl;
                    temp = temp->next;
                }
            }

            // Update index to start looking for vertices of smallest degree
            // Start looking at the degree_just_deleted -1, since no vertex can have its degree reduced by more than 1
            if (start != 0) {
                start = degree_just_deleted - 1;
            }

            deleted_ctr++;
            cout << endl;
            // printDegreeDLLs(degree_DLLs, V);
            cout << endl;
        }
        else { // No vertices at current degree list, move on
            cout << "No vertices of degree " << start << endl;
            start++;
        }
    }

    auto time_stop = std::chrono::high_resolution_clock::now();
    auto time_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(time_stop - time_start).count();
    cout << "SMALLEST_LAST time_elapsed: " << time_elapsed << endl << endl;
    // Record max degree when deleted and runtime
    stats[2] = static_cast<double>(max_degree_when_deleted);
    runtimes[1] = time_elapsed;
}

// Method to generate a smallest original degree last ordering.
void smallestOriginalDegreeLastOrdering(Vertex vertices[], DoublyLinkedList degree_DLLs[], int V) {
    // degree_DLLs, the array of doubly linked lists, already has each vertex stored by degree.
    // Therefore, we can make a smallest original degree last ordering by simply traversing all of the nodes
    // in degree_DLLs and adding them in order to the ordering list.
    Vertex v;
    for (int i = 0; i < V; i++) {
        DLLNode* curr = degree_DLLs[i].head;
        while (curr != nullptr) {
            v = vertices[curr->data];
            vertices[v.id].deleted = true;
            vertices[v.id].degree_when_deleted = v.degree;
            vertices[v.id].ordering->insert(curr->data);
            curr = curr->next;
        }
    }
}

// Method to generate a uniform random ordering.
void uniformRandomOrdering(Vertex vertices[], int V) {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniform_distribution(0,V-1);
    int rand;
    int ctr = 0;
    while (ctr != V) { // Until all edges have been added to the ordering:
        rand = uniform_distribution(rng); // Generate a random vertex id [0, V-1]
        if (!vertices[rand].deleted) { // If the vertex has not been added to the ordering yet,
            vertices[rand].deleted = true; // Mark it as deleted
            vertices[rand].degree_when_deleted = vertices[rand].degree;
            vertices[rand].ordering->insert(rand); // Add it to the ordering.
            ctr++;
        }
    }
}

// Method to generate a largest original degree last ordering.
void largestOriginalDegreeLastOrdering(Vertex vertices[], DoublyLinkedList degree_DLLs[], int V) {
    // degree_DLLs, the array of doubly linked lists, already has each vertex stored by degree.
    // Therefore, we can make a largest original degree last ordering by simply traversing all of the nodes
    // in degree_DLLs in reverse order and adding them to the ordering list.
    // This ordering is essentially the opposite of the smallest original degree last ordering.
    Vertex v;
    for (int i = V - 1; i >= 0; i--) {
        DLLNode* curr = degree_DLLs[i].head;
        while (curr != nullptr) {
            v = vertices[curr->data];
            vertices[v.id].deleted = true;
            vertices[v.id].degree_when_deleted = v.degree;
            vertices[v.id].ordering->insert(curr->data);
            curr = curr->next;
        }
    }
}

// Greedy coloring algorithm.
// For each vertex, picks the lowest integer available in the range [0, V-1] that is not adjacent to the vertex.
// Returns an int designating the highest color value used to keep track of the number of colors needed to color.
// Implementation modified from https://www.geeksforgeeks.org/graph-coloring-set-2-greedy-algorithm/
void colorVertices(Vertex vertices[], int coloring_order[], int V) {
    int max_color = INT_MIN;

    cout << "Coloring order: ";
    for (int i = 0; i < V; i++) {
        cout << coloring_order[i] << " ";
    }
    cout << endl;

    // Set the first vertex in the ordering to the first color
    // All other vertices in the graph already had their colors initialized to -1 when first allocating the structs
    vertices[coloring_order[0]].color = 0;
    cout << "----------------Colored " << coloring_order[0] << " with color 0----------------" << endl;

    // color_adjacent array denotes whether or not a color [0, V-1] is adjacent to the current vertex
    // If color_adjacent[i] == true, then that color has been assigned to an adjacent vertex,
    // so we cannot pick it for the current vertex.
    // Initialize all colors' adjacency to false
    bool color_adjacent[V];
    for (int i = 0; i < V; i++) {
        color_adjacent[i] = false;
    }

    // For each vertex in the graph (except the first one already colored):
    for (int i = 1; i < V; i++) {
        cout << "----------------Coloring " << vertices[coloring_order[i]].id << "----------------" << endl;
        // Check if any of the adjacent vertices have been colored
        Node* temp = vertices[coloring_order[i]].edges.head;
        while (temp != nullptr) {
            cout << "--------Adjacent: " << temp->data << "--------" << endl;
            if (vertices[temp->data].color != -1) { // If the adjacent vertex has already been colored,
                cout << vertices[temp->data].id << " has already been colored with color " << vertices[temp->data].color << endl;
                color_adjacent[vertices[temp->data].color] = true; // Mark its color as adjacent
            }
            else { // The adjacent vertex has not been colored, so we may ignore it
                cout << vertices[temp->data].id << " has not yet been colored, skipping " << endl;
            }
            temp = temp->next;
        }

        // Find the first non-adjacent color in the color_adjacent matrix
        // If an index in color_adjacent is false, that means that that color is not adjacent to the current vertex
        // Therefore, we take the first color we can find that is false and assign it to the vertex
        int c;
        for (c = 0; c < V; c++) {
            if (!color_adjacent[c]) {
                break;
            }
        }

        // Color v with the first non-adjacent color
        vertices[coloring_order[i]].color = c;
        cout << "Colored " << vertices[coloring_order[i]].id << " with color " << vertices[coloring_order[i]].color << endl;

        // Update max_color used if necessary
        if (c > max_color) {
            max_color = c;
        }

        // Reset all of the adjacent colors to false for next iteration
        temp = vertices[coloring_order[i]].edges.head;
        while (temp != nullptr) {
            if (vertices[temp->data].color != -1) {
                cout << vertices[temp->data].id << " was already colored; resetting color_adjacent["<< vertices[temp->data].color << "] to false" << endl;
                color_adjacent[vertices[temp->data].color] = false;
            }
            temp = temp->next;
        }
    }

    // Record total colors used
    stats[0] = static_cast<double>(max_color + 1);
}

int main(int argc, char** argv) {

    const int V = atoi(argv[1]); // MAX = 10,000
    const int E = atoi(argv[2]); // MAX = 2,000,000
    const string G = argv[3]; // COMPLETE | CYCLE | RANDOM (with DIST below)
    const string DIST = argv[4]; // UNIFORM | SKEWED | YOURS
    const string ORDERING = argv[5]; // SMALLEST_LAST| SMALLEST_ORIGINAL_DEGREE_LAST | RANDOM | LARGEST_ORIGINAL_DEGREE_LAST
    cout << V << " " << E << " " << G << " " << DIST << " " << ORDERING << endl;
//
//    // Array of structs to hold vertex information
//    Vertex vertices[V];
//    // Array of doubly linked lists to be indexed by degree. Each list contains integer IDs
//    // Each DLL node only holds an integer corresponding to the vertex's ID.
//    DoublyLinkedList degree_DLLs[V];
//    // Linked list to keep track of order vertices are deleted in
//    LinkedList ordering;
//
//    // Allocate and initialize V vertices
//    for (int i = 0; i < V; i++) {
//        Vertex v;
//        v.id = i;
//        v.edges.head = nullptr;
//        v.degree = 0;
//        v.color = -1;
//        v.deleted = false;
//        v.degree_when_deleted = -1;
//        v.original_degree = 0;
//        v.degree_DLL = nullptr;
//        v.ordering = &ordering;
//        vertices[i] = v;
//    }
//
//    // Allocate and initialize V doubly linked lists
//    for (int i = 0; i < V; i++) {
//        DoublyLinkedList DLL;
//        degree_DLLs[i] = DLL;
//    }
//
//    // Generate the graph specified by the command line arguments.
//    generateGraph(vertices, degree_DLLs, V, E, G, DIST);

//    // Display the adjacency list describing the graph.
//    printAdjList(vertices, V);
//
//    // Print the degree-indexed doubly linked list for the graph.
//    printDegreeDLLs(degree_DLLs, V);

//   //  Panopto example

//    vertices[0].edges.insert(1);
//    vertices[0].degree = 1;
//    vertices[0].original_degree = 1;
//    vertices[0].degree_DLL = &degree_DLLs[1];
//
//    vertices[1].edges.insert(0);
//    vertices[1].edges.insert(2);
//    vertices[1].edges.insert(3);
//    vertices[1].degree = 3;
//    vertices[1].original_degree = 3;
//    vertices[1].degree_DLL = &degree_DLLs[3];
//
//    vertices[2].edges.insert(1);
//    vertices[2].edges.insert(3);
//    vertices[2].degree = 2;
//    vertices[2].original_degree = 2;
//    vertices[2].degree_DLL = &degree_DLLs[2];
//
//    vertices[3].edges.insert(1);
//    vertices[3].edges.insert(2);
//    vertices[3].degree = 2;
//    vertices[3].original_degree = 2;
//    vertices[3].degree_DLL = &degree_DLLs[2];
//
//    degree_DLLs[1].insert(0);
//    degree_DLLs[2].insert(3);
//    degree_DLLs[2].insert(2);
//    degree_DLLs[3].insert(1);
//
//    printAdjList(vertices, V);
//    printDegreeDLLs(degree_DLLs, V);
//
    // Perform the ordering algorithm specified by the command line argument.
//    if (ORDERING == "SMALLEST_LAST") {
//        smallestLastVertexOrdering(vertices, degree_DLLs, V);
//    }
//    else if (ORDERING == "SMALLEST_ORIGINAL_DEGREE_LAST") {
//        smallestOriginalDegreeLastOrdering(vertices, degree_DLLs, V);
//    }
//    else if (ORDERING == "RANDOM") {
//        uniformRandomOrdering(vertices, V);
//    }
//    else if (ORDERING == "LARGEST_ORIGINAL_DEGREE_LAST") {
//        largestOriginalDegreeLastOrdering(vertices, degree_DLLs, V);
//    }
//
//    // Generate array for coloring order (reverse of ordering).
//    int coloring_order[V];
//    Node* temp = ordering.head;
//    for (int i = V-1; i >= 0; i--) {
//        if (temp != nullptr) {
//            coloring_order[i] = temp->data;
//        }
//        temp = temp->next;
//    }
//
//    // Perform greedy coloring algorithm on the graph and get the highest color value used.
//    colorVertices(vertices, coloring_order, V);
//    cout << endl;
//
//    // Calculate average original degree.
//    double sum = 0.0;
//    for (int i = 0; i < V; i++) {
//        sum += static_cast<double>(vertices[i].original_degree);
//    }
//    double avg_original_degree = sum / static_cast<double>(V);
//    stats[1] = avg_original_degree;
//
//    // Calculate size of terminal clique for smallest last vertex ordering.
//    if (ORDERING == "SMALLEST_LAST") {
//        int terminal_clique_size = 1;
//        for (int i = 1; i < V; i++) {
//            cout << "Comparing " << vertices[coloring_order[i-1]].degree_when_deleted << ", " << vertices[coloring_order[i]].degree_when_deleted << endl;
//            if (vertices[coloring_order[i]].degree_when_deleted <= vertices[coloring_order[i-1]].degree_when_deleted) {
//                break;
//            }
//            else {
//                terminal_clique_size++;
//            }
//        }
//        stats[3] = static_cast<double>(terminal_clique_size);
//    }
//
//    ofstream output("output.csv");
//    output << V << ", " << E << ", " << G << ", " << DIST << ", " << ORDERING << '\n';
//    // Record ordering, coloring order, and colors assigned to the graph.
//    recordOrderingAndColoring(output, vertices, ordering, coloring_order, V, ORDERING);
//
//    output.close();

    // TODO: Timing/data for ordering & coloring
    // SMALLEST_LAST
    // // Sparse graph
    // // Cycle
    // // Small Random uniform graph
    // // Small Random Skewed graph
    // // Small Random Normal graph
    // // Medium Random Uniform graph
    // // Medium Random Skewed graph
    // // Medium random Normal graph
    // // Large Random Uniform graph
    // // Large Random Skewed graph
    // // Large Random Normal graph
    // // Complete
    // SMALLEST_ORIGINAL_DEGREE_LAST
    // // Sparse graph
    // // Cycle
    // // Small Random uniform graph
    // // Small Random Skewed graph
    // // Small Random Normal graph
    // // Medium Random Uniform graph
    // // Medium Random Skewed graph
    // // Medium random Normal graph
    // // Large Random Uniform graph
    // // Large Random Skewed graph
    // // Large Random Normal graph
    // // Complete
    // UNIFORM RANDOM
    // // Sparse graph
    // // Cycle
    // // Small Random uniform graph
    // // Small Random Skewed graph
    // // Small Random Normal graph
    // // Medium Random Uniform graph
    // // Medium Random Skewed graph
    // // Medium random Normal graph
    // // Large Random Uniform graph
    // // Large Random Skewed graph
    // // Large Random Normal graph
    // // Complete
    // LARGEST_ORIGINAL_DEGREE_LAST
    // // Sparse graph
    // // Cycle
    // // Small Random uniform graph
    // // Small Random Skewed graph
    // // Small Random Normal graph
    // // Medium Random Uniform graph
    // // Medium Random Skewed graph
    // // Medium random Normal graph
    // // Large Random Uniform graph
    // // Large Random Skewed graph
    // // Large Random Normal graph
    // // Complete


    int V_values[10] = { 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000};
    int E_values[10] = { 200000, 400000, 600000, 800000, 1000000, 1200000, 1400000, 1600000, 1800000, 2000000 };
    int V_values2[10] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    int E_values2[10] = { 20, 40, 60, 80, 100, 120, 140, 160, 180, 200};

    int V_values3[10] = { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
    int E_values3[10] = { 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000};

    int complete_V_values[10] = { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };
    int complete_V_values2[10] = { 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000 };


    // TODO: if I want to run timing overnight, I can put this n loop inside another one from 0 to 4
    // TODO: and run everything that number of times, cout-ing a barrier for different trials

    int n;
    //for (n = 0; n < 10; n++) {
//        int V = V_values3[n];
//        int E = E_values3[n];
//        string G = "RANDOM";
//        string DIST = "YOURS";
        // cout << V << " " << E << " " << G << " " << DIST << endl;

        Vertex vertices[V];
        DoublyLinkedList degree_DLLs[V];
        LinkedList ordering;
        for (int i = 0; i < V; i++) {
            Vertex v;
            v.id = i;
            v.edges.head = nullptr;
            v.degree = 0;
            v.color = -1;
            v.deleted = false;
            v.degree_when_deleted = -1;
            v.original_degree = 0;
            v.degree_DLL = nullptr;
            v.ordering = &ordering;
            vertices[i] = v;
        }
        for (int i = 0; i < V; i++) {
            DoublyLinkedList DLL;
            degree_DLLs[i] = DLL;
        }

//        auto time_start = std::chrono::high_resolution_clock::now();
        generateGraph(vertices, degree_DLLs, V, E, G, DIST);
//        cout << "---------------------------------------------- ADJ LIST ----------------------------------------------" << endl;
//        printAdjList(vertices, V);
//        printDegreeDLLs(degree_DLLs, V);

//        auto time_stop = std::chrono::high_resolution_clock::now();
//        auto time_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(time_stop - time_start).count();
//        //runtimes[n] = time_elapsed;
        //runtime_data_file << V_values[n] << ", " << E_values[n] << ", " << runtimes[n] << '\n';
        //cout << V_values[n] << ", " << E_values[n] << ", " << time_elapsed << '\n';
    //}
//    runtime_data_file.close();




    cout << '\n' << "Done" << endl;
    // Sometimes <ntdll!RtlRaiseException> disassembler signal when freeing G?
    return 0;
}
