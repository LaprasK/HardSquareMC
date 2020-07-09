#ifndef Corr_Included
#define Corr_Included
#include <math.h>
#include <vector>
#include <string>
#include <fstream>

std::vector<double> orient_corr(const std::vector<std::vector<double> >& angles);

inline std::vector<std::vector<double> > load_orient(std::string filename);

#endif