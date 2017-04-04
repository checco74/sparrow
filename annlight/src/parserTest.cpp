#include "ParameterParser.h"
#include "ParserException.h"
#include <iostream>

int main (int argc, char** argv)
{
  ParameterParser* parser = new ParameterParser();

  try
  {
    parser->addParameter("name");
    parser->addParameter("alpha_plus", "0.01");
    parser->addParameter("test");

    std::cout << "before: " << std::endl;
    parser->writeCurrentParameters(std::cout);
    std::cout << std::endl;

    parser->readParameters("test.cfg");

    std::cout << "getter:" << std::endl;
    std::cout << "name = \"" << parser->getParameterValue("name") << "\"" << std::endl;
    std::cout << "alpha_plus = \"" << parser->getParameterValue("alpha_plus") << "\"" << std::endl;
    std::cout << "test = \"" << parser->getParameterValue("test") << "\"" << std::endl << std::endl;
    std::cout << "after: " << std::endl;
    parser->writeCurrentParameters(std::cout);
  }
  catch (ParserException e)
  {
    std::cerr << "Exception occured: " << e.what() << std::endl;
    exit(1);
  }

  delete parser;
}
