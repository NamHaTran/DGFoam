#ifndef BCREADER_H
#define BCREADER_H

#include <string>

namespace readVectorBC {
    void u(std::string mode);
}

namespace readScalarBC {
    void p(std::string mode);
    void T(std::string mode);
}
#endif // BCREADER_H
