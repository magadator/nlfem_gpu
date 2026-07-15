#include "file_io.h"

using namespace std;

namespace KARMA {

std::string File_IO::FiletoString (const string& file_name) {
    AssertFileExists(file_name);
    ifstream f(file_name);
    stringstream buffer;
    buffer << f.rdbuf();
    return buffer.str();
}

} // end namespace KARMA
