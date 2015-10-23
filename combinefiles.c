#include <stdio.h>
#include <stdlib.h>
#include <string.h>



int main(int argc, char **argv)
{

  int N;
  N = atoi(argv[1]);
  printf("Number of cores used=%d\n",N);
  char ch;
  FILE* res;
  res = fopen("solution_total.dat","w");
  char fname[200];
  FILE* fp;
  int i;
  for(i=0;i<N;i++)
     {
       sprintf(fname,"solution%d.dat",i);

       fp=fopen(fname,"r");
       while((ch = fgetc(fp))!=EOF)
             fputc(ch,res);
      fclose(fp);

     }
   fclose(res);

return 0;

