#include "PreprocessorContext.h"
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "string_utilities.h"
#include "BuiltInFunctionList.h"
//wvn added
namespace KARMA
{
    PreprocessorContext::PreprocessorContext(int* pArgC_in, char*** pArgV_in, int procID_in)
    {
        pArgC = pArgC_in;
        pArgV = pArgV_in;
        procID = procID_in;
        //invocationStart = "${";
        //invocationEnd = "}";
        preProcessIndicator = '#';
        invocationSymbol = "$";
        functionInvocationSymbol = "@";
        invocationStart = "{";
        invocationEnd = "}";
        functionArgDelimiter = ',';
        dlmChar = ' ';
        CreateCLIDefs();
    }
    
    void PreprocessorContext::CreateCLIDefs(void)
    {
        for (int i = 0; i < *pArgC; i++)
        {
            char* strP = *((*pArgV)+i);
            std::string argString(strP);
            AddCLIDefIfParseable(argString);
        }
    }
    
    bool PreprocessorContext::DefinitionExists(std::string str)
    {
        std::map<std::string,std::string>::iterator it = definitions.find(str);
        return (it != definitions.end());
    }
    
    void PreprocessorContext::AddCLIDefIfParseable(std::string cliArg)
    {
        if (cliArg.length() < 2) return;
        if ((cliArg[0]!='-')||(cliArg[1]!='D')) return;
        std::string definition = cliArg.substr(2, cliArg.length()-2);
        size_t eq = definition.find("=");
        if (eq==std::string::npos) {CreateDefinition(definition, ""); return;}
        std::string key = definition.substr(0, eq);
        std::string val = definition.substr(eq+1, definition.length()-eq);
        ValidateKey(key);
        CreateDefinition(key, val);
    }
    
    bool PreprocessorContext::SymbolPrecedes(std::string str, std::string preceder, std::string test)
    {
        if (!StringContains(str, preceder)) return false;
        if (!StringContains(str, test)) return true;
        return (str.find(preceder) < str.find(test));
    }
    
    void PreprocessorContext::SplitFunction(std::string str, std::string* pre, std::string* post, std::string* func, std::vector<std::string>* args, int level, bool* success, std::string origLine)
    {
        size_t idxFuncSym = str.find(functionInvocationSymbol);
        *pre = str.substr(0, idxFuncSym);
        std::string noPre = str.substr(1+idxFuncSym, str.length()-idxFuncSym-1);
        if (!StringContains(noPre, invocationStart)) HaltKarma("Cannot find arguments for function call " + str);
        size_t iStart = noPre.find(invocationStart);
        size_t iEnd;
        int bracketLevel = 0;
        for (size_t i = iStart; i < noPre.length(); i++)
        {
            if (PositionIsStart(i, noPre, invocationStart)) bracketLevel++;
            if (PositionIsEnd(i, noPre, invocationEnd)) bracketLevel--;
            if (bracketLevel==0) {iEnd = i; break;}
        }
        std::string contents = noPre.substr(iStart+1, iEnd-iStart-1);
        *func = noPre.substr(0, iStart);
        *post = noPre.substr(iEnd+1, noPre.length()-iEnd-1);
        BuildArgs(args, contents, functionArgDelimiter, level, success, origLine);
    }
    
    void PreprocessorContext::BuildArgs(std::vector<std::string>* args, std::string line, char delimiterIn, int level, bool* success, std::string origLine)
    {
        std::vector<std::string> output;
        std::vector<size_t> positions;
        positions.push_back(0);
        int sectionLevel = 0;
        for (size_t i = 1; i < line.length()-1; i++)
        {
            if (PositionIsStart(i, line, invocationStart))  sectionLevel++;
            if (PositionIsEnd(i, line, invocationEnd)) sectionLevel--;
            if ((line[i]==delimiterIn) && (sectionLevel==0)) positions.push_back(i);
        }
        positions.push_back(line.length());
        for (int i = 0; i < positions.size()-1; i++)
        {
            size_t start = positions[i];
            if (i >0) start++;
            size_t end = positions[i+1];
            std::string elem = line.substr(start, end-start);
            if (elem.length()>0)
            {
                args->push_back(ResolveWithinContext(elem, level+1, success, origLine));
            }
        }
    }
    
    std::string PreprocessorContext::EvalFunction(std::string& func, std::vector<std::string>& args, std::string origLine)
    {
        if (BuiltIns::builtInFunctions.Exists(func))
        {
            try {return BuiltIns::builtInFunctions[func](args);}
            catch (BuiltIns::BuiltinException a)
            {
                HaltKarma("Error evaluating built-in function \"" + func + "\" from the following line:\n >> " + origLine + "\nMessage: " + a.what());
            }
        }
        else
        {
            HaltKarma("Cannot find built-in function \"" + func + "\"");
            return "NONE";
        }
    }
    
    bool PreprocessorContext::StringContains(std::string str, std::string c)
    {
        return (str.find(c) != std::string::npos);
    }
    
    void PreprocessorContext::ValidateKey(std::string key)
    {
        std::string acceptable = "1234567890_qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM";
        for (size_t i = 0; i < key.length(); i++)
        {
            if (acceptable.find(key[i])==std::string::npos) ErrorThrow("Found invalid character in key \"" + key + "\".");
        }
    }
    
    void PreprocessorContext::CreateDefinition(std::string key, std::string val)
    {
        if (DefinitionExists(key)) ErrorThrow("Found multiple definitions for \"" + key + "\".");
        definitions.insert({key, val});
    }
    
    std::string PreprocessorContext::GetDefinition(std::string key, bool* success)
    {
        return GetDefinition(key);
    }
    
    std::string PreprocessorContext::GetDefinition(std::string key)
    {
        if (!DefinitionExists(key)) ErrorThrow("No definition for \"" + key + "\".");
        return definitions[key];
    }
    
    void PreprocessorContext::HaltKarma(std::string message)
    {
        ErrorThrow(message);
    }
    
    void PreprocessorContext::ErrorThrow(std::string message)
    {
        std::cout << "------------------------------------------------------" << std::endl;
        std::cout << "Terminate called from " << __FILE__ << ", line " << __LINE__ << ":" << std::endl;
        std::cout << message << std::endl;
        std::cout << "------------------------------------------------------" << std::endl;
        abort();
    }
    
    std::string PreprocessorContext::ResolveWithinContext(std::string str)
    {
        bool dummy;
        return ResolveWithinContext(str, 0, &dummy, str);

    }
    
    std::string PreprocessorContext::ResolveWithinContext(std::string str_in, int level, bool* success, std::string origLine)
    {
        std::string str = str_in;
        std::string fullInvocationStart = invocationSymbol + invocationStart;
        if (str.length()==0) return str;
        if (StringContains(str, functionInvocationSymbol) && SymbolPrecedes(str, functionInvocationSymbol, invocationSymbol))
        {
            std::string pre, post, func;
            std::vector<std::string> args;
            SplitFunction(str, &pre, &post, &func, &args, level, success, origLine);
            std::string output = pre + ResolveWithinContext(EvalFunction(func, args, origLine), level+1, success, origLine) + ResolveWithinContext(post, level+1, success, origLine);
            if (procID == 0 && level == 0)
            {
                std::cout << "[P] Preprocessor resolution: " << str_in << " -> " << output << std::endl;
            }
            return output;
        }
        else
        {
            AssertBracketConsistency(str);
            size_t start, end;
            start = str.find(fullInvocationStart);
            if (start==std::string::npos)
            {
                return str;
            }
            int bracketLevel = 1;
            for (size_t i = start+fullInvocationStart.length(); i < str.length(); i++)
            {
                if (PositionIsStart(i, str)||PositionIsStart(i, str, invocationStart)) bracketLevel++;
                if (PositionIsEnd(i, str)||PositionIsEnd(i, str, invocationEnd)) bracketLevel--;
                if (bracketLevel==0) {end = i; break;}
            }
            std::string pre = str.substr(0, start);
            std::string med = str.substr(start+fullInvocationStart.length(), end-start-fullInvocationStart.length());
            std::string post = str.substr(end+invocationEnd.length(), str.length()-end-invocationEnd.length());
            std::string defn = GetDefinition(ResolveWithinContext(med, level+1, success, origLine), success);
            std::string after = ResolveWithinContext(post, level+1, success, origLine);
            std::string output = pre + defn + after;
            if (procID == 0 && level == 0)
            {
                std::cout << "[P] Preprocessor resolution: " << str_in << " -> " << output << std::endl;
            }
            return output;
        }
    }
    
    /*
    std::string PreprocessorContext::ResolveWithinContext(std::string str, int level)
    {        
        if (str.length()==0) return str;
        AssertBracketConsistency(str);
        size_t start, end;
        start = str.find(invocationStart);
        if (start==std::string::npos)
        {
            return str;
        }
        int bracketLevel = 1;
        for (size_t i = start+invocationStart.length(); i < str.length(); i++)
        {
            if (PositionIsStart(i, str)) bracketLevel++;
            if (PositionIsEnd(i, str)) bracketLevel--;
            if (bracketLevel==0) {end = i; break;}
        }
        std::string pre = str.substr(0, start);
        std::string med = str.substr(start+invocationStart.length(), end-start-invocationStart.length());
        std::string post = str.substr(end+invocationEnd.length(), str.length()-end-invocationEnd.length());
        std::string output = pre + GetDefinition(ResolveWithinContext(med, level+1)) + ResolveWithinContext(post, level+1);
        if (level==0 && procID==0)
        {
            std::cout << "[P] Preprocessor: " << str << " --> " << output << std::endl;
        }
        return output;
    }*/
    
    /*void PreprocessorContext::AssertBracketConsistency(std::string str)
    {
        int level = 0;
        for (size_t i = 0; i < str.length(); i++)
        {
            if (PositionIsStart(i, str)) level++;
            if (PositionIsEnd(i, str))
            {
                level--;
                if (level<0) ErrorThrow("Invocation \"" + str + "\" has inconsistent brackets.");
            }
        }
        if (level!=0) ErrorThrow("Invocation \"" + str + "\" has inconsistent brackets.");
    }*/
    
    void PreprocessorContext::AssertBracketConsistency(std::string str)
    {
        int level = 0;
        for (size_t i = 0; i < str.length(); i++)
        {
           if (PositionIsStart(i, str, invocationStart)) level++;
           if (PositionIsEnd(i, str, invocationEnd))
           {
               level--;
               if (level<0) HaltKarma("Invocation \"" + str + "\" has inconsistent brackets (level < 0).");
           }
        }
        if (level!=0) HaltKarma("Invocation \"" + str + "\" has inconsistent brackets (level != 0).");
    }
    
    bool PreprocessorContext::PositionIsStart(size_t i, std::string str, std::string toFind)
    {
        if (i < toFind.length()-1) return false;
        return (str.substr(i-toFind.length()+1, toFind.length())==toFind);
    }
    
    bool PreprocessorContext::PositionIsStart(size_t i, std::string str)
    {
        std::string fullInvocationStart = invocationSymbol + invocationStart;
        return PositionIsStart(i, str, fullInvocationStart);
    }

    bool PreprocessorContext::PositionIsEnd(size_t i, std::string str, std::string toFind)
    {
        if (i+toFind.length()-1 >= (str.length())) return false;
        return (str.substr(i, toFind.length()) == toFind);
    }
    bool PreprocessorContext::PositionIsEnd(size_t i, std::string str)
    {
        return PositionIsEnd(i, str, invocationEnd);
    }
    
    void PreprocessorContext::ParseDirectives(std::string & str)
    {
        std::vector<std::string> lines = StringUtils::Explode(str, '\n');
        std::string defineString = "#define ";
        for (int i = 0; i < lines.size(); i++)
        {
            size_t pos = lines[i].find(defineString);
            if (pos != std::string::npos)
            {
                std::string argument = lines[i].substr(pos+defineString.length(), lines[i].length()-defineString.length());
                size_t space = argument.find(" ");
                if (space == std::string::npos)
                {
                    StringUtils::StripWhiteSpaces(argument);
                    ValidateKey(argument);
                    CreateDefinition(argument, "");
                }
                else
                {
                    std::string newDef = argument.substr(0, space);
                    std::string newValue = argument.substr(space+1, argument.length()-space-1);
                    StringUtils::StripWhiteSpaces(newDef);
                    StringUtils::StripWhiteSpaces(newValue);
                    ValidateKey(newDef);
                    newValue = ResolveWithinContext(newValue);
                    CreateDefinition(newDef, newValue);
                }
            }
        }
    }
    
    PreprocessorContext::~PreprocessorContext(void)
    {
        definitions.clear();
    }
}
