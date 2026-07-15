#ifndef PREPROCESS_CONTEXT_H
#define PREPROCESS_CONTEXT_H
#include <string>
#include <map>
#include <vector>
namespace KARMA
{
    class PreprocessorContext
    {
        public:
            PreprocessorContext(int* pArgC_in, char*** pArgV_in, int procID_in);
            ~PreprocessorContext(void);
            void ParseDirectives(std::string & str);
            void CreateCLIDefs(void);
            bool DefinitionExists(std::string str);
            void CreateDefinition(std::string key, std::string val);
            void ErrorThrow(std::string message);
            void HaltKarma(std::string message);
            void AddCLIDefIfParseable(std::string cliArg);
            void ValidateKey(std::string key);
            std::string ResolveWithinContext(std::string str);
            std::string ResolveWithinContext(std::string str, int level, bool* success, std::string origLine);
            void AssertBracketConsistency(std::string str);
            bool PositionIsStart(size_t i, std::string str);
            bool PositionIsStart(size_t i, std::string str, std::string toFind);
            bool PositionIsEnd(size_t i, std::string str);
            bool PositionIsEnd(size_t i, std::string str, std::string toFind);
            std::string GetDefinition(std::string key, bool* success);
            std::string GetDefinition(std::string key);
            bool SymbolPrecedes(std::string str, std::string preceder, std::string test);
            void BuildArgs(std::vector<std::string>* args, std::string line, char delimiterIn, int level, bool* success, std::string origLine);
            bool StringContains(std::string str, std::string c);
            void SplitFunction(std::string str, std::string* pre, std::string* post, std::string* func, std::vector<std::string>* args, int level, bool* success, std::string origLine);
            std::string EvalFunction(std::string& func, std::vector<std::string>& args, std::string origLine);
            int GetArgC(void) {return *pArgC;}
            char** GetArgV(void) {return *pArgV;}
        private:
            int* pArgC;
            char*** pArgV;
            int procID;
            char functionArgDelimiter, dlmChar;
            std::map<std::string, std::string> definitions;
            std::string preProcessIndicator;
            std::string invocationSymbol;
            std::string functionInvocationSymbol;
            std::string invocationStart;
            std::string invocationEnd;
    };
}
#endif