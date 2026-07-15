#ifndef KARMA_UTILS_FILE_IO_H
#define KARMA_UTILS_FILE_IO_H

#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include "screen.h"

namespace KARMA {
    
class File_IO {
    public:
        static inline bool FileExists (const std::string& file_name) {
            struct stat buffer;                                
            return (stat (file_name.c_str(), &buffer) == 0);
        };
        static inline void AssertFileExists (const std::string& file_name) {
            if (!FileExists(file_name)) Screen::MasterError("File \""+file_name+"\" does not exist");
        };
        static std::string FiletoString (const std::string& file_name);
};

} // end namespace KARMA

#endif