#include <string>
#include <map>
#include <iostream>

#ifndef _PARAMETERPARSER_H_
#define _PARAMETERPARSER_H_

/**
 * Class implementing a parameter parser capable of reading configuration files. In these configuration files the identifier and the value
 * of a parameter have to be seperated by a whitespace. Line comments are allowed in the configuration files.
 *
 * @brief A parser to read parameters from files.
 */
class ParameterParser
{
  public:
    ParameterParser();
    ~ParameterParser();

    void addParameter(std::string name, std::string defaultValue = "");
    void setCommentString(std::string commentString);

    void readParameters(std::string filename);
    void writeParameterFile(std::string filename);
    void writeCurrentParameters(std::ostream& stream);

    std::string getParameterValue(std::string name) const;

    /** @brief String representation for the internal maps representing a value which has not been set and for which no default value exists. */
    static const std::string NO_DEFAULT_VALUE_STRING;

    /** @brief The default string representing the beginning of a comment. */
    static const std::string DEFAULT_COMMENT_STRING;

  protected:
    std::string& trim(std::string& str) const;

    /** @brief A map associating a parameter name with a parameter value. */
    std::map<std::string, std::string> nameValueMap;

    /** @brief String representing the beginning of a comment. */
    std::string commentString;
};

#endif
