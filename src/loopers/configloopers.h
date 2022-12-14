#ifndef CONFIGLOOPERS_H
#define CONFIGLOOPERS_H

class ConfigParser;

namespace Loopers {

class FineAlign;
class FineAlignDut;
class CoarseAlign;
class CoarseAlignDut;
class NoiseScan;
class Synchronize;

void configFineAlign(const ConfigParser& config, FineAlign& fineAlign);
void configFineAlign(const ConfigParser& config, FineAlignDut& fineAlign);
void configCoarseAlign(const ConfigParser& config, CoarseAlign& coarseAlign);
void configCoarseAlign(const ConfigParser& config, CoarseAlignDut& coarseAlign);
void configNoiseScan(const ConfigParser& config, NoiseScan& noiseScan);
void configSynchronize(const ConfigParser& config, Synchronize& sync);

}

#endif // CONFIGLOOPERS_H
