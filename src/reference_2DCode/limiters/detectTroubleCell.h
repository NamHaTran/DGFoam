#ifndef DETECTTROUBLECELL_H
#define DETECTTROUBLECELL_H
#include <vector>

namespace troubledCellDetector
{
    bool minmodDetector(std::vector<double> InputVector, double condition);

    bool PerssonPeraireDetector(int elem);
}
#endif // DETECTTROUBLECELL_H
