#include "screen.h"
#include "string_utilities.h"
#include <iostream>
using namespace std;

namespace KARMA {
   
void Screen::Line(void) {
    if (MPI_Ops::Rank()!=0) return;
    for (int i=0;i<LINE_LENGTH;++i) cout << LINE_CHAR;
    cout << endl;
    return;
}

void Screen::Info(const std::string& message, const int indent_level) {
    cout << INFO_TAG;
    for (int i=0;i<indent_level;++i) cout << INDENT_STR;
    cout << " (proc=" + to_string(MPI_Ops::Rank()) + ") " << message << endl;
    return;
}

void Screen::Warning(const std::string& message, const int indent_level) {
    cout << WARNING_TAG;
    for (int i=0;i<indent_level;++i) cout << INDENT_STR;
    cout << " (proc=" + to_string(MPI_Ops::Rank()) + ") " << message << endl;
    return;
}   

void Screen::Error(const std::string& message, const int indent_level) {
    cout << ERROR_TAG;
    for (int i=0;i<indent_level;++i) cout << INDENT_STR;
    cout << " (proc=" + to_string(MPI_Ops::Rank()) + ") " << message << endl;
#ifdef KARMA_MPI
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#endif
    abort();
    return;
}   

void Screen::MasterInfo(const std::string& message, const int indent_level) {
    if (MPI_Ops::Rank()!=0) return;
    cout << INFO_TAG;
    if (indent_level!=0) {
        for (int i=0;i<indent_level;++i) cout << INDENT_STR;
        cout << " " << endl;
    }
    cout << message << endl;
    return;
}

void Screen::MasterWarning(const std::string& message, const int indent_level) {
    if (MPI_Ops::Rank()!=0) return;
    cout << WARNING_TAG;
    if (indent_level!=0) {
        for (int i=0;i<indent_level;++i) cout << INDENT_STR;
        cout << " " << endl;
    }
    cout << message << endl;
    return;
}

void Screen::MasterError(const std::string& message, const int indent_level) {
    if (MPI_Ops::Rank()!=0) return;
    cout << ERROR_TAG;
    if (indent_level!=0) {
        for (int i=0;i<indent_level;++i) cout << INDENT_STR;
        cout << " " << endl;
    }
    cout << message << endl;
#ifdef KARMA_MPI
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#endif
    abort();
    return;
}

} // end namespace KARMA
