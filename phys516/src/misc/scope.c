#include <stdio.h>
int i = 1; /* Global variable */

void function_a() {
  printf("In function_a: i = %d\n", i);
  i = i + 1;
}

int function_b() {
  int i = 3; /* Local variable overrides the global one */
  printf("In function_b: i = %d\n", i);
  i = i + 1;
  return i;
}

int main() {
  int li;
  function_a();
  li = function_b();
  printf("In main: global i = %d; local i from function b = %d\n", i, li);
  return 0;
}