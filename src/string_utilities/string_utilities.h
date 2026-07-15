#ifndef KARMA_UTILS_STRING_H
#define KARMA_UTILS_STRING_H

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace KARMA {

template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}
class StringUtils {
public:
    static void StripComments   (std::string& str);
    static void StripWhiteSpaces(std::string& str);
    static void StripEmptyLines (std::string& str);
    static std::string ExtractInBetween(std::string& str,std::string& begin_tag,std::string &end_tag);
    static std::vector<std::string> Explode(std::string& str, char delim);
    static std::string int2str(const int &number);    
    static std::vector<int>    str2ivec (std::string& str);
    static std::vector<double> str2dvec (std::string& str);
    static std::vector<std::string> str2svec(std::string &str, std::string beginStr, std::string endStr);
    static std::vector<std::string> SplitButIgnoreEscapes(std::string data, std::string delimiter, char start, char end);
};

static inline void GetFormatSubstrings(std::vector<std::string>& subStrings, std::string templateStr)
{
    std::string delimiter = "{}";
    std::string templateStrCopy = templateStr;
    size_t pos = 0;
    std::string token;
    while ((pos = templateStrCopy.find(delimiter)) != std::string::npos)
    {
        token = templateStrCopy.substr(0, pos);
        subStrings.push_back(token);
        templateStrCopy.erase(0, pos + delimiter.length());
    }
    subStrings.push_back(templateStrCopy);
}

template <typename T> static inline void strunformat_recursive (std::string data, std::string& templateStr, std::vector<std::string>& subStrings, int& lev, T& t)
{
    std::string& startToken = subStrings[lev];
    std::string& endToken = subStrings[lev+1];
    size_t spos = data.find(startToken);
    if (spos==std::string::npos) return;
    size_t epos = data.find(endToken, spos+startToken.length());
    if (endToken=="") epos = data.length();
    if (epos==std::string::npos) return;
    size_t subStrStart = spos+startToken.length();
    size_t subStrLength = epos-spos-startToken.length();
    t = data.substr(subStrStart, subStrLength);
    lev++;
}

template <typename T, typename... Ts> static inline void strunformat_recursive (std::string data, std::string& templateStr, std::vector<std::string>& subStrings, int& lev, T& t, Ts&... ts)
{
    std::string& startToken = subStrings[lev];
    std::string& endToken = subStrings[lev+1];
    size_t spos = data.find(startToken);
    if (spos==std::string::npos) return;
    size_t epos = data.find(endToken, spos+startToken.length());
    if (endToken=="") epos = data.length();
    if (epos==std::string::npos) return;
    size_t subStrStart = spos+startToken.length();
    size_t subStrLength = epos-spos-startToken.length();
    t = data.substr(subStrStart, subStrLength);
    std::string dataSubStr = data.substr(epos, data.length()-epos);
    lev++;
    strunformat_recursive(dataSubStr, templateStr, subStrings, lev, ts...);
}

template <typename... Ts> static inline void strunformat (std::string data, std::string templateStr, Ts&... ts)
{
    std::vector<std::string> subStrings;
    GetFormatSubstrings(subStrings, templateStr);
    if ((sizeof...(Ts))!=(subStrings.size()-1))
    {
		std::cout << "Killing from file " << __FILE__ << ", line " << __LINE__ << ": mismatch in prototype delimiter." << std::endl;
        abort();
    }
    int lev = 0;
    strunformat_recursive(data, templateStr, subStrings, lev, ts...);
}

  
} // end namespace KARMA

#endif
