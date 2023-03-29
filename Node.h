#ifndef PROJECT_NODE_H
#define PROJECT_NODE_H

// Node code acquired from https://www.geeksforgeeks.org/program-to-implement-singly-linked-list-in-c-using-class/.

class Node {
public:
    int data;
    Node* next;

    Node() {
        data = 0;
        next = nullptr;
    }

    Node(int data) {
        this->data = data;
        this->next = nullptr;
    }
};

#endif //PROJECT_NODE_H
