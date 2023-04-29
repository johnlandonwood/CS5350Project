#ifndef PROJECT_DOUBLYLINKEDLIST_H
#define PROJECT_DOUBLYLINKEDLIST_H
#include "DLLNode.h"

// Doubly linked-list code modified from https://www.programiz.com/dsa/doubly-linked-list.

class DoublyLinkedList {

public:
    DLLNode* head;
    int listSize = 0;

    DoublyLinkedList() {
        head = nullptr;
    }

    void insert(int);
    int size();
    void print();
    void printReverse();
    void remove(int);

};

// Function to insert node at the head of the DLL
void DoublyLinkedList::insert(int data) {
    DLLNode* newNode = new DLLNode(data);

    newNode->next = head;
    newNode->prev = nullptr;

    if (head != nullptr) {
        head->prev = newNode;
    }

    head = newNode;
    this->listSize++;
}

// Function to remove a node from the DLL given the ID (data field) of the node to delete
void DoublyLinkedList::remove(int id_to_delete) {
    if (head == nullptr) {
        std::cout << "List empty, cannot delete." << std::endl;
        return;
    }

    // Traverse to the node to be deleted
    // TODO: Will this traversal make it worse than O(V+E)?
    DLLNode *temp = head;
    while (temp != nullptr) {
        if (temp->data == id_to_delete) {
            break;
        }
        temp = temp->next;
    }

    // If node to be deleted is head, update head pointer
    if (temp->data == head->data) {
        head = head->next;
    }
    // If node to be deleted is not the last node in the list, update next's prev pointer
    if (temp->next != nullptr) {
        temp->next->prev = temp->prev;
    }
    // If node to be deleted is not the head, update prev's next pointr
    if (temp->prev != nullptr) {
        temp->prev->next = temp->next;
    }

    this->listSize--;
}

// Method to return size of the DLL
int DoublyLinkedList::size() {
    return this->listSize;
}

// Method to print the DLL in order
void DoublyLinkedList::print() {
    if (head == nullptr) {
        std::cout << "Empty.";
        return;
    }

    DLLNode* temp = head;
    while (temp != nullptr) {
        std::cout << temp->data;
        if (temp->next != nullptr) {
            std::cout << " --> ";
        }
        temp = temp->next;
    }
    //std::cout << std::endl;
}

// Method to print the DLL in reverse order
void DoublyLinkedList::printReverse() {
    if (head == nullptr) {
        std::cout << "Empty.";
        return;
    }

    DLLNode* temp = head;
    while (temp->next != nullptr) {
        temp = temp->next;
    }

    while (temp != head) {
        std::cout << temp->data << " <-- ";
//        if (temp->prev != nullptr) {
//            std::cout << " <-- ";
//        }
        temp = temp->prev;
    }

    std::cout << std::endl;
}


#endif //PROJECT_DOUBLYLINKEDLIST_H
