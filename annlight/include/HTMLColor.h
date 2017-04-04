#ifndef _HTMLCOLOR_H_
#define _HTMLCOLOR_H_

class HTMLColor
{
  public:
    HTMLColor(double r, double g, double b);
    HTMLColor(unsigned char r, unsigned char g, unsigned char b);
    ~HTMLColor();
    
    std::string getHTMLCode() const;
    
  protected:
    static const unsigned char MAX_HTML_COLOR = 255;
    
    unsigned char r;
    unsigned char g;
    unsigned char b;
    
    std::string toHex(unsigned char c) const;
    char dec2hexChar(unsigned char c) const;
};

#endif
