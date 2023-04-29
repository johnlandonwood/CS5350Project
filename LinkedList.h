#ifndef PROJECT_LINKEDLIST_H
#define PROJECT_LINKEDLIST_H
#include "Node.h"

// List code acquired from https://www.geeksforgeeks.org/program-to-implement-singly-linked-list-in-c-using-class/.

class LinkedList {

public:
    Node* head;
    int listSize = 0;

    // Default constructor
    LinkedList() {
        head = nullptr;
    }

    void insert(int);
    int size();
    void print();
};

// Method to insert a new node at the end of the list
void LinkedList::insert(int data) {
    Node* newNode = new Node(data);

    // If list is empty, assign to head
    if (head == nullptr) {
        head = newNode;
        this->listSize++;
        return;
    }

    // Insert node at end
    Node* temp = head;
    while (temp->next != nullptr) {
        temp = temp->next;
    }

    temp->next = newNode;
    this->listSize++;
}

// Method to return size of the list
int LinkedList::size() {
    return this->listSize;
}

// Method to print the linked list
void LinkedList::print(){
    // Check for empty list.
    if (head == nullptr) {
        std::cout << "Empty.";
        return;
    }

    Node* temp = head;
    while (temp != nullptr) {
        std::cout << temp->data;
        if (temp->next != nullptr) {
            std::cout << " --> ";
        }
        temp = temp->next;
    }

}

#endif //PROJECT_LINKEDLIST_H
