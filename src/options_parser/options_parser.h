#ifndef KARMA_UTILS_OPTIONS_PARSER_H
#define KARMA_UTILS_OPTIONS_PARSER_H
//#define _GLIBCXX_USE_C99 1
#include <string>
#include <map>
#include <vector>
#include "file_io.h"
#include "string_utilities.h"
#include "PreprocessorContext.h"
#include <vector>
#include <iostream>
#include <iterator>
#ifdef KARMA_PTL
#include "PTL.h"
#endif
//
#define ENTRY_BOOL    1
#define ENTRY_INT     2
#define ENTRY_DOUBLE  3
#define ENTRY_STRING  4


namespace KARMA {

class OptionEntry {
protected:
    std::string value_;
public:
    std::string name;
    std::string GetValue(void) {return value_;}
    std::string GetPreProcessedValue(void) {return pContext->ResolveWithinContext(value_);}
#ifdef KARMA_PTL
    PTL::PropertyTree* ptlTreePtr=NULL;
#endif
    std::string completeString;
    OptionEntry (void) {value_=""; is_found_=false;pContext=NULL;}
    OptionEntry (std::string value,bool is_found=true) {value_=value; is_found_= (value_!=""); pContext=NULL;}
    void SetContext(PreprocessorContext* p) {pContext = p;}
    operator std::string();
    operator int ();
    operator std::vector<std::string> ();
    operator double ();
    operator bool ();
    operator std::vector<int> ();
    operator std::vector<double> ();
    bool is_found_;
    bool usePTL;
private:
    PreprocessorContext* pContext = NULL;

};

class OptionsSection {
protected:
    std::string name_;
    int  level_;
    bool is_found_;
    OptionsSection *parent_;
    std::string completeString;
#ifdef KARMA_PTL
    PTL::PropertyTree* ptlTreePtr=NULL;
#endif
    bool usePTL;
    std::string content_;
    PreprocessorContext* pContext = NULL;
    std::map<std::string,OptionsSection> children_;
    std::map<std::string,std::string> entry_string_;
    
    std::map<std::string,OptionEntry> entries_;
    
    void CreateChildSection(const std::string &name);
    
public:
    void OutputToFile(std::ofstream& myfile);
  OptionsSection& Section(const std::string &name);
  bool            IsFound(void);
  OptionEntry     Get(const std::string &name);
  void GetAllEntries(std::vector<std::string>* names);
  std::string& GetContent(void)
  {
      return this->content_;
  }
};

class OptionsParser: public OptionsSection {
private:
    std::string file_name_;
#ifdef KARMA_PTL
    PTL::PropertyTree ptlObj;
#endif
public:
    // Contructor
    OptionsParser(PreprocessorContext* p);
    int SetFile(const std::string &file_name);
    void OutputPreprocessedInput(std::string filename);
};
} // end namespace KARMA

#endif
