#include "HTMLColor.h"

HTMLColor::HTMLColor(double r, double g, double b)
{
  this->r = (unsigned char)(r * MAX_HTML_COLOR);
  this->g = (unsigned char)(g * MAX_HTML_COLOR);
  this->b = (unsigned char)(b * MAX_HTML_COLOR);
}

HTMLColor::HTMLColor(unsigned char r, unsigned char g, unsigned char b)
{
  this->r = r;
  this->g = g;
  this->b = b;
}

HTMLColor::~HTMLColor()
{
}

std::string HTMLColor::getHTMLCode() const
{
  std::string s = "";
  s += "#";
  s += this->toHex(this->r);
  s += this->toHex(this->g);
  s += this->toHex(this->b);
  
  return s;
}

std::string HTMLColor::toHex(unsigned char c) const
{
  unsigned first = c / 16;
  unsigned second = c % 16;
  
  std::string tmp = "";
  tmp += this->dec2hexChar(first);
  tmp += this->dec2hexChar(second);
  
  return tmp;
}

char HTMLColor::dec2hexChar(unsigned char c) const
{
  if (c <= 9)
    return (c + 48);
  else if (c <= 15)
    return (c + 55);
  else
  {
    std::cerr << "ERROR: Character is out of bounds!" << std::endl;
    exit(1);
    return 0;
  }
}
