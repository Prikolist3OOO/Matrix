#include <stdio.h>
#include <stdlib.h>
struct node{
    int a;
    struct node* left;
    struct node* right;
};
struct node* init_node(int);
struct node* generate_BST(int*, int, int);
int main(){
    int array[] = {1, 3, 7, 8, 9, 10};
    struct node* derevo = generate_BST(array, 0, 5);
    
    return 0;
}
struct node* init_node(int a){
    struct node* t;
    t = (struct node*) malloc(sizeof(struct node));
    t->a = a;
    t->left = NULL;
    t->right = NULL;
    return t;
}
struct node* generate_BST(int arr[], int start, int end){
    if (start > end) return NULL;
    int mid = (start + end)/2;
    struct node* root = init_node(arr[mid]);
    root->left = generate_BST(arr, start, mid-1);
    root->right = generate_BST(arr, mid+1, end);
    return root;
}

void new_elem(int n, struct node* root){
    if (n == root->a) return;
    if(n < root->a){
        if (root->left == NULL){
            root->left = init_node(n);
        } else new_elem(n, root->left);
    } else if(root->right == NULL){
        root->right = init_node(n);
    } else new_elem(n, root->right);
}
