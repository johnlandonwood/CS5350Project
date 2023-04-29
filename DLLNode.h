#ifndef PROJECT_DLLNODE_H
#define PROJECT_DLLNODE_H


class DLLNode {
public:
    int data;
    DLLNode* next;
    DLLNode* prev;

    DLLNode() {
        data = 0;
        next = nullptr;
        prev = nullptr;
    }

    DLLNode(int data) {
        this->data = data;
        this->next = nullptr;
        this->prev = nullptr;
    }

};


#endif //PROJECT_DLLNODE_H
