#include "mbimage.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>

void MBImage::save(const char *filename)
{
  char fName[8192];
  strcpy(fName,filename);
  int l=strlen(filename)-4;
  if(strcmp(&filename[l],".ppm"))
    strcat(fName,".ppm");
  printf("Saving to : %s\n",fName);

  FILE *out = fopen(fName, "wb");
  fprintf(out, "P6\n%d %d\n255\n", width, height);
  fwrite(&rgb[0], width * height * 3, 1, out);
  fflush(out);
  fclose(out);
}
