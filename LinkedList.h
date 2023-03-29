#ifndef PROJECT_LINKEDLIST_H
#define PROJECT_LINKEDLIST_H
#include "Node.h"

// List code acquired from https://www.geeksforgeeks.org/program-to-implement-singly-linked-list-in-c-using-class/.

class LinkedList {
    Node* head;

public:
    // Default constructor
    LinkedList() {
        head = nullptr;
    }

    // Function to insert in order
    void insert(int);

    void print();
};

// In-order insertion code taken from https://www.youtube.com/watch?v=p0u_SFoZbl8.
void LinkedList::insert(int data) {
    Node* newNode = new Node(data);

    // Assign to head
    if (head == nullptr) {
        head = newNode;
        return;
    }

    // Insert node at end
    Node* temp = head;
    while (temp->next != nullptr) {
        temp = temp->next;
    }

    temp->next = newNode;
}

void LinkedList::print(){
    Node* temp = head;

    // Check for empty list.
    if (head == nullptr) {
        std::cout << "List empty" << std::endl;
        return;
    }

    while (temp != nullptr) {
        std::cout << temp->data << " ";
        temp = temp->next;
    }
    std::cout << std::endl;
}

#endif //PROJECT_LINKEDLIST_H
