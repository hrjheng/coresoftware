#ifndef DFINDERSPHENIX_DFINDERFUN4ALL_H
#define DFINDERSPHENIX_DFINDERFUN4ALL_H

#include "DFinderConfig.h"

#include <string>

class SubsysReco;

SubsysReco *makeDFinderSphenix(const std::string &name,
                               const DFinderConfig &config);

#endif
