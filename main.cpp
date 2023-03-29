#include <iostream>
#include "LinkedList.h"
#include <string>
#include <random>
using namespace std;

LinkedList* generateGraph(LinkedList* adjList, const int V, const int E, const string G, const string DIST) {
    // to generate random graphs:
    // for i = 1 to E
    // generate v1 - random number between 0 and V - 1 (here is where uniform vs. skewed. vs yours comes in)
    // generate v2 - random number between 0 and V - 1 (here is where uniform vs. skewed. vs yours comes in)
    // the edge goes between those two vertices
    // add v2 to v1's adjacency list and v1 to v2's adjacency list

    // how to avoid a number being the same?
    // if we generate  0 twice in a row, we can't add 0 as an edge to itself
    // add check where if v2 == v1 after generation, then go again?
    // or just always generate two numbers - if they are equal, go again, if they are not equal add the edge


    int zero = 0;
    int one = 0;
    int two = 0;
    int three = 0;
//    for (int i = 0; i < 1000; i++) {
//        int rand = dist(rng);
//        if (rand == 0) zero++;
//        else if (rand == 1) one++;
//        else if (rand == 2) two++;
//        else if (rand == 3) three++;
//    }
//    cout << zero << " " << one << " " << two << " " << three << endl;

    //int rand = dist(rng);



    if (G == "COMPLETE") {
        // 0 --> 1
        // 0 --> 2
        // 0 --> 3
        // 1 --> 2
        // 1 --> 3
        // 2 --> 3
    }
    else if (G == "CYCLE") {
        // 0 --> 1
        // 1 --> 2
        // 2 --> 3
        // 3 --> 0
    }
    else if (G == "RANDOM") {
        if (DIST == "UNIFORM") {
            std::random_device rd;
            std::mt19937 rng(rd());
            std::uniform_int_distribution<int> uniform_distribution(0,V-1);
            int v1 = 0;
            int v2 = 0;
            int i = 0;
            while (i < E) { // How do we make sure that an edge does not already exist?
                i++;
                v1 = uniform_distribution(rng);
                v2 = uniform_distribution(rng);
                if (v2 == v1) {
                    cout << "Edge generation collision: " << v1 << " --> " << v2 << endl;
                    i--;
                }
                // else if () { // if edge does not already exist
                    // add edge to graph
                    cout << v1 << " --> " << v2 << endl;
                }
            }

//            for (int i = 0; i < E; i++) {
//                v1 = uniform_distribution(rng);
//                v2 = uniform_distribution(rng);
//                if (v2 == v1) {
//                    cout << "Edge generation collision: " << v1 << " --> " << v2 << endl;
//
//                }
//                cout << v1 << "-->" << v2 << endl;
//            }

        }
        else if (DIST == "SKEWED") {
            // use skewed community package
        }
        else if (DIST == "YOURS") {

        }
    }


    return adjList;
}

LinkedList* test (LinkedList* adjList) {

    cout << "In test" << endl;
    adjList[0].insert(5);

    return adjList;
}

int main(int argc, char** argv) {

    const int V = atoi(argv[1]); // MAX = 10,000
    const int E = atoi(argv[2]); // MAX = 2,000,000
    const string G = argv[3]; // COMPLETE | CYCLE | RANDOM (with DIST below)
    const string DIST = argv[4]; // UNIFORM | SKEWED | YOURS

    cout << V << " " << E << " " << G << " " << DIST << endl;

    // Adjacency list - array of lists
    LinkedList adjList[V] = { };
    generateGraph(adjList, V, E, G, DIST);

    // Create one list for each vertex
//    for (int i = 0; i < V; i++) {
//        LinkedList list;
//        for (int j = 0; j < 5; j++) {
//            list.insert(j);
//        }
//        adjList[i] = list;
//    }
//
//
//
//
//
//
//    for (int i = 0; i < V; i++) {
//        adjList[i].print();
//    }
//
//    test(adjList);


//    for (int i = 0; i < 5; i++) {
//        list0.insert(i);
//        list1.insert(i);
//        list2.insert(i);
//        list3.insert(i);
//        list4.insert(i);
//    }

    // generateCycle(int V)
    // generateComplete(int V)
    // generateRandom(int V, int E, string DIST)
    // method to generate complete graph
    // method to generate




    cout << "done";
    return 0;
}
