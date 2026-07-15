/*
 * nlfem - Nonlinear FEM Solver
 * Author  : Gokul G. Anugrah
 * Contact : gokulanugrah@gmail.com
 *
 * options_parser.cc — Structured input file parser.
 *   Parses the nested bracket input format used by nlfem.
 *   Supports preprocessing directives, arithmetic expressions,
 *   vector-valued entries, and optional PTL back-end.
 */
#include "options_parser.h"

using namespace std;

namespace KARMA {
    
OptionsParser::OptionsParser(PreprocessorContext* p) {
    name_    = "_ROOT_";
    level_   = 0;
    content_ = "";
    pContext = p;
    usePTL = false;
#ifdef KARMA_PTL
    ptlTreePtr = &ptlObj;
    PTL::Interactive i(p->GetArgC(), p->GetArgV(), &ptlObj);
#endif
    this->completeString = "";
    return;
}
void OptionsParser::OutputPreprocessedInput(std::string filename)
{
    std::ofstream myfile;
    myfile.open(filename);
    myfile << "#PREPROCESSED" << std::endl;
    for (auto& child:this->children_)
    {
        child.second.OutputToFile(myfile);
    }
    myfile.close();
}
bool EndsWith(const string& a, const string& b) {
    if (b.size() > a.size()) return false;
    return std::equal(a.begin() + a.size() - b.size(), a.end(), b.begin());
}
int OptionsParser::SetFile(const string &file_name) {
    usePTL = false;
#ifdef KARMA_PTL
    usePTL = EndsWith(file_name, ".ptl");
    if (usePTL)
    {
        ptlObj.ReadInputFileToTreeData(file_name);
        ptlObj.ResolveAllStrings();
    }
#endif
    if (!usePTL)
    {
        file_name_ = file_name;
        content_   = File_IO::FiletoString(file_name);
        if (pContext!=NULL) pContext->ParseDirectives(content_);
        StringUtils::StripWhiteSpaces(content_);
        StringUtils::StripComments(content_);
        StringUtils::StripEmptyLines(content_);
    }
    return 1;
}

void OptionsSection::OutputToFile(std::ofstream& myfile)
{
    std::string levelStrStart = "";
    std::string levelStrEnd = "";
    for (int i = 0; i < this->level_; i++)
    {
        levelStrStart += "[";
        levelStrEnd += "]";
    }
    if (this->entries_.size() == 0 && this->children_.size() == 0) return;
    myfile << levelStrStart << this->name_ << levelStrEnd << std::endl;
    for (auto& entry: this->entries_)
    {
        OptionEntry& val = entry.second;
        if (val.GetValue().length()>0)
        {
            myfile << val.name << " = " << val.GetPreProcessedValue() << std::endl;
        }
    }
    for (auto& child: this->children_)
    {
        OptionsSection& sec = child.second;
        sec.OutputToFile(myfile);
    }
    myfile << levelStrStart << "/" << this->name_ << levelStrEnd << std::endl;
}

void OptionsSection::CreateChildSection(const string &name) {
    OptionsSection new_section;
    new_section.name_   = name;
    new_section.level_  = level_+1;
    new_section.parent_ = this;
    new_section.is_found_ = false;
    if (level_==0) { new_section.completeString = name; }
    else { new_section.completeString = this->completeString + "." + name; }
    new_section.pContext = this->pContext;
    new_section.usePTL = this->usePTL;
#ifdef KARMA_PTL
    new_section.ptlTreePtr = this->ptlTreePtr;
#endif
        
    string sec_begin_tag = "";
    string sec_end_tag   = "";
    for (int i=0;i<new_section.level_;++i) {
        sec_begin_tag += "[";
        sec_end_tag   += "[";
    }
    sec_begin_tag += name;
    sec_end_tag   += "/"+name;
    for (int i=0;i<new_section.level_;++i) {
        sec_begin_tag += "]";
        sec_end_tag   += "]";
    }   
    
    if (usePTL)
    {
#ifdef KARMA_PTL
        new_section.content_ = "";
        new_section.is_found_ = ptlTreePtr->Query(new_section.completeString).Found();
        new_section.content_ = ptlTreePtr->Query(new_section.completeString).Content();
#endif
    }
    else
    {
        new_section.content_ = StringUtils::ExtractInBetween(content_,sec_begin_tag,sec_end_tag);
        if (new_section.content_ != "") {
            new_section.is_found_ = true;
        }
        StringUtils::StripEmptyLines(new_section.content_);
    }
    children_.insert(pair<string,OptionsSection>(name,new_section));
    return;
}

void OptionsSection::GetAllEntries(std::vector<std::string>* names)
{
    if (!usePTL)
    {
        std::vector<std::string> lines = StringUtils::Explode(this->GetContent(), '\n');
        for (int i = 0; i < lines.size(); i++)
        {
            string begin_tag = "=";
            string temp = "#"+lines[i];
            string temp2 = "#";
            string name      = StringUtils::ExtractInBetween(temp,temp2,begin_tag);
            names->push_back(name);
        }
    }
    else
    {
#ifdef KARMA_PTL
      auto qry = ptlTreePtr->Query(completeString);
      std::vector<std::string> stuff = qry.TreeSection()->GetTerminalSections();
      for (const auto& s:stuff) names->push_back(s);
#endif
    }
}

OptionsSection& OptionsSection::Section(const string &name) {
    // Check if the section is already created
    if (children_.find(name)==children_.end()) CreateChildSection(name);
    return children_[name];
}
  
  bool OptionsSection::IsFound(void) {
    return is_found_;
}

  OptionEntry OptionsSection::Get(const string &name) {
    // Check if the option is already parsed
    map<string,OptionEntry>::iterator it = entries_.find(name);
    if (it==entries_.end()) {
      // New entry
      string begin_tag = name+"=";
      string end_tag   = "\n";
      
      string value     = "";
      if (!usePTL) value = StringUtils::ExtractInBetween(content_,begin_tag,end_tag);
      else
      {
#ifdef KARMA_PTL
        auto qry = ptlTreePtr->Query(completeString + "." + name);
        value = qry.Content();
#endif
      }
      OptionEntry new_entry(value);
      new_entry.completeString = completeString + "." + name;
#ifdef KARMA_PTL
      new_entry.ptlTreePtr = this->ptlTreePtr;
#endif
      new_entry.usePTL = usePTL;
      new_entry.name = name;
      new_entry.SetContext(pContext);
      entries_.insert(pair<string,OptionEntry>(name,new_entry));
      return new_entry;
    } else {
      return it->second;
    }
  }
    OptionEntry::operator std::string()
    {
        if (!usePTL)
        {
            if(pContext) value_ = pContext->ResolveWithinContext(value_);
        }
        else
        {
#ifdef KARMA_PTL
              std::string output = ptlTreePtr->Query(completeString);
              return output;
#endif
        }
        return value_;
    }
    OptionEntry::operator int ()                 
    {
      try
      {
        if (!usePTL)
        {
            if(pContext) value_ = pContext->ResolveWithinContext(value_);
            return stoi(value_);
        }
        else
        {
#ifdef KARMA_PTL
              int output = ptlTreePtr->Query(completeString);
              return output;
#endif
        }
      }
      catch (...)
      {
          pContext->HaltKarma("Could not parse \"" + value_ + "\" as integer! Variable name: " + name);
      }
    }
    OptionEntry::operator std::vector<std::string> ()
    {
         if (!usePTL)
         {
             if(pContext) value_ = pContext->ResolveWithinContext(value_);
             return StringUtils::str2svec(value_, "(", ")");
         }
         else
         {
             return StringUtils::str2svec(value_, "[", "]");
         }
    }
    OptionEntry::operator double ()
    {
      try
      {
          if (!usePTL)
          {
              if(pContext) value_ = pContext->ResolveWithinContext(value_);
              return stod(value_);
          }
          else
          {
#ifdef KARMA_PTL
              double output = ptlTreePtr->Query(completeString);
              return output;
#endif
          }
      }
      catch (...)
      {
          pContext->HaltKarma("Could not parse \"" + value_ + "\" as double! Variable name: " + name);
      }    
    }
    OptionEntry::operator bool ()                
    {
        if (!usePTL)
        {
            if(pContext) value_ = pContext->ResolveWithinContext(value_);
            return (value_ == "true");
        }
        else
        {
#ifdef KARMA_PTL
            bool output = ptlTreePtr->Query(completeString);
            return output;
#endif
        }
    }
    OptionEntry::operator std::vector<int> ()
    {
      try
      {
          if (!usePTL)
          {
              if(pContext) value_ = pContext->ResolveWithinContext(value_);
              return StringUtils::str2ivec(value_);
          }
          else
          {
#ifdef KARMA_PTL
              std::vector<int> output = ptlTreePtr->Query(completeString);
              return output;
#endif
          }
      }
      catch (...)
      {
          pContext->HaltKarma("Could not parse \"" + value_ + "\" as integer list! Variable name: " + name);
      }
      
    }
    OptionEntry::operator std::vector<double> ()
    {
      try
      {
          if (!usePTL)
          {
              if(pContext) value_ = pContext->ResolveWithinContext(value_);
              return StringUtils::str2dvec(value_);
          }
          else
          {
#ifdef KARMA_PTL
            std::vector<double> output = ptlTreePtr->Query(completeString);
            return output;
#endif
          }
      }
      catch(...)
      {
          pContext->HaltKarma("Could not parse \"" + value_ + "\" as double list! Variable name: " + name);
      }
    }
} // end namespace KARMA
