#include "string_utilities.h"
#include <locale>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

namespace KARMA {

void StringUtils::StripComments(string& str) {
    // This doesn't remove the lines
    // If the comment is in it's own line, this will leave an empty line
    size_t beg=0;
    while (beg!=string::npos) {
        beg = str.find("#");
        if (beg!=string::npos) str.erase(beg,str.find("\n",beg)-beg);
    }
    return;    
}

void StringUtils::StripWhiteSpaces(string& str) {
    // This doesn't remove line breaks
    string whitespaces(" \t\f\v\r");                                          
    size_t beg=0;
    while (beg!=string::npos) {                                               
        beg=str.find_first_of(whitespaces);                                  
        if (beg!=string::npos) str.erase(beg,1);                           
    }
    beg = 0;
    // remove line continuations
    std::string continuation = "\\\n";
    while (beg!=string::npos) {
        beg=str.find(continuation);
        if (beg!=string::npos) str.erase(beg,continuation.length());
    }
    return;
}

void StringUtils::StripEmptyLines(string& str) {
    size_t beg=0;                                                             
    // If the first line is empty, remove it
    if (str.substr(0,1)=="\n") str.erase(0,1);
    while (beg!=string::npos) {                                               
        beg=str.find("\n\n");                                  
        if (beg!=string::npos) str.erase(beg,1);                           
    }     
    return;
}    
    
string StringUtils::ExtractInBetween(string& str,string& begin_tag,string &end_tag) {
    size_t begin_pos = str.find(begin_tag);
    if (begin_pos==string::npos) return "";
    size_t end_pos   = str.find(end_tag,begin_pos);
    if (end_pos==string::npos) return "";
    // Shift begin_pos by the tag length
    begin_pos += begin_tag.length();
    return str.substr(begin_pos,end_pos-begin_pos);
}
string StringUtils::int2str(const int &number) {
    std::stringstream dummy;
    dummy << number;
    //return dummy as string
    return dummy.str();
}

vector<string> StringUtils::Explode(string &str, char delim) {
    vector<string> result;
    istringstream iss(str);

    for (string token; getline(iss, token, delim); ) {
        result.push_back(move(token));
    }
    return result;
}

std::vector<int> StringUtils::str2ivec(string &str) {
    string begin_tag = "(";
    string end_tag   = ")";
    string str_in = ExtractInBetween(str,begin_tag,end_tag);
    vector<string> str_exp = Explode(str_in,',');
    vector<int> result(str_exp.size());
    for (int i=0;i<str_exp.size();++i) {
        result[i] = stoi(str_exp[i]);
    }
    return result;
}

std::vector<double> StringUtils::str2dvec(string &str) {
    string begin_tag = "(";
    string end_tag   = ")";
    string str_in = ExtractInBetween(str,begin_tag,end_tag);
    vector<string> str_exp = Explode(str_in,',');
    vector<double> result(str_exp.size());
    for (int i=0;i<str_exp.size();++i) {
        result[i] = stod(str_exp[i]);
    }
    return result;
}

std::vector<std::string> StringUtils::str2svec(string &str, string beginStr, string endStr) {
    string begin_tag = beginStr;
    string end_tag   = endStr;
    string str_in = ExtractInBetween(str,begin_tag,end_tag);
    vector<string> str_exp = SplitButIgnoreEscapes(str_in,",",'[',']');
    return str_exp;
}
  
  std::vector<std::string> StringUtils::SplitButIgnoreEscapes(std::string data, std::string delimiter, char start, char end)
  {
	std::vector<std::string> subStrings;
	int level = 0;
	for (int i = 0; i < data.length(); i++)
      {
		if (data[i]==start) level++;
		if (data[i]==end) level--;
      }
	if (level!= 0)
      {
		std::cout << "Killing from file " << __FILE__ << ", line " << __LINE__ << ": mismatch in escape characters." << std::endl;
		std::cout << "Line:\n" << data << "\nEscapes:\n" << start << end << std::endl;
		abort();
      }
    std::string templateStrCopy = data;
    size_t pos = 0;
    std::string token;
    while ((pos = templateStrCopy.find(delimiter)) != std::string::npos)
      {
        token = templateStrCopy.substr(0, pos);
        subStrings.push_back(token);
        templateStrCopy.erase(0, pos + delimiter.length());
      }
    subStrings.push_back(templateStrCopy);
	int recombineStart = -1;
	int recombineEnd = -1;
	std::vector<std::string> subStringsRecombined;
	std::string line;
	for (int i = 0; i < subStrings.size(); i++)
      {
		line = subStrings[i];
		if (subStrings[i].find(start)!=std::string::npos) recombineStart = i;
		if (subStrings[i].find(end)!=std::string::npos) recombineEnd = i;
		if (!(recombineStart>=0 || recombineEnd >= 0))
          {
			subStringsRecombined.push_back(line);
          }
		if (recombineStart>=0 && recombineEnd >= 0)
          {
			line = "";
			for (int j = recombineStart; j <= recombineEnd; j++)
              {
				line += subStrings[j] + ((j==recombineEnd)?(""):(","));
              }
			recombineStart = -1;
			recombineEnd = -1;
			subStringsRecombined.push_back(line);
          }		
      }
	
    return subStringsRecombined;
  }
  
} // end of namespace KARMA
