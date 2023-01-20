#include <stdio.h>
#define NMAX 512
#define MAX 255
#define MAXLINE 1024

int main() {
  int image[NMAX][NMAX],nr,nc,r,c,n,work;
  FILE *f;
  char line[MAXLINE];

  f=fopen("Lenna512x512.pgm","r");

  fgets(line,MAXLINE,f);
  fgets(line,MAXLINE,f);
  fgets(line,MAXLINE,f);
  sscanf(line,"%d %d",&nc,&nr);
  fgets(line,MAXLINE,f);
  sscanf(line,"%d",&n);

  for (r=0; r<nr; r++) {
    for (c=0; c<nc; c++) {
      image[r][c] = (int)fgetc(f);
    }
  }

  fclose(f);

  f=fopen("Lenna256x256.pgm","w");

  fprintf(f,"P5\n");
  fprintf(f,"# Simple image test\n");
  fprintf(f,"%d %d\n",nc/2,nr/2);
  fprintf(f,"%d\n",MAX);

  for (r=0; r<nr/2; r++) {
    for (c=0; c<nc/2; c++) {
      work=( image[2*r  ][2*c  ]+image[2*r  ][2*c+1]
            +image[2*r+1][2*c  ]+image[2*r+1][2*c+1])/4;
      fputc((char)work,f);
    }
  }

  fclose(f);
}
