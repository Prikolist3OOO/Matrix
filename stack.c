#include <stdio.h>
#include <stdlib.h>

struct stack{
    int m[100];
    int top
};
void init(struct stack*);
void push(struct stack* , double);

int main(){

}

void init(struct stack *t) {
  t->top = -1;
}

void push(struct stack *t, double x) {
  if(t->top < 100) {
    t->m[t->top] = x;
    t->top++;
  }
}

double pop(struct stack *t) {
  double e;
  if((t->top) > 0)
  {
    t->top--;
    e = t->m[t->top];
    return e;
  }
}