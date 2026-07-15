#ifndef KARMA_UTILS_SCREEN_H
#define KARMA_UTILS_SCREEN_H

#include <vector>
#include <string>
#include "mpi_ops.h"
//
using namespace std;
namespace KARMA {
    
#define INFO_TAG    "\033[0;32m[I]\033[0m "
#define WARNING_TAG "\033[0;31m[W]\033[0m "
#define ERROR_TAG   "\033[1;31m[E]\033[0m "
#define INDENT_STR  "...."
#define LINE_CHAR   "="
#define LINE_LENGTH 80

class Screen {
    public:
        // Constructor takes in the mpi_rank
        static void Line          (void);
        static void Info          (const std::string& message, const int indent_level=0);
        static void Warning       (const std::string& message, const int indent_level=0);
        static void Error         (const std::string& message, const int indent_level=0);
        static void MasterInfo    (const std::string& message, const int indent_level=0);
        static void MasterWarning (const std::string& message, const int indent_level=0);
        static void MasterError   (const std::string& message, const int indent_level=0);
};
} // end namespace KARMA

#endif
