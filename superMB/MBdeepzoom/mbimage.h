#ifndef MBIMAGE_H
#define MBIMAGE_H

#include <vector>

// simple RGB24 image

class MBImage
{
public:
    typedef unsigned char uint8_t;
  int width;
  int height;
  std::vector<uint8_t> rgb;

  // construct empty image
  MBImage(int width, int height)
  : width(width)
  , height(height)
  , rgb(width * height * 3)
  { }

  // plot a point
  void plot(int x, int y, uint8_t r, uint8_t g, uint8_t b)
  {
    //TODO: check bounds?
    int k = (y * width + x) * 3;
      //if(x*x+y*y<100){rgb[k++] = 0;rgb[k++] = 0;rgb[k++] = 0; return;}

    rgb[k++] = r;
    rgb[k++] = g;
    rgb[k++] = b;
  }

  // plot a point with bound checking
  void plot1(int x, int y, uint8_t r, uint8_t g, uint8_t b)
  {
      if (x<0 || y<0 || x>=width || y>=height) return;
      plot(x,y,r,g,b);
  }

  //resize
  bool resize(int w, int h){
      if(w==width && h==height) return true;//nothing to do
      if(w<0 || w>16384 || h<0 || h>16384) return false;
      width = w; height = h;
      unsigned int newSize = width * height * 3;
      rgb.resize(newSize);
      return true;//OK!
  }

  //Is it Okay
  bool OK(){return (rgb.size() == (unsigned int)(width * height * 3));}

  // get pointer to the rendered image
  unsigned char * getPtr2Img(){return &rgb[0];}

  // save to PPM format
  void save(const char *filename);
};

#endif // MBIMAGE_H
