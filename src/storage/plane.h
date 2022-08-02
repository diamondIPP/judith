#ifndef PLANE_H
#define PLANE_H

#include <vector>
#include <map>
#include <iostream>
#include <sstream>

namespace Storage {

class Hit;
class Cluster;

/*******************************************************************************
 * Plane class contains the hits and clusters for one sensor plane as well as
 * a plane number.
 */

class Plane
{
private:
  unsigned int _planeNum;
  unsigned int _numHits;
  std::vector<Hit*> _hits;
  unsigned int _numClusters;
  std::vector<Cluster*> _clusters;
  std::map< std::string, std::vector<float>* > m_waveforms;
  std::map< std::string, float > m_wfSamplingWidths;

protected:
  void addHit(Hit* hit);
  void addCluster(Cluster* cluster);
  Plane(unsigned int planeNum);
  ~Plane() { ; }

public:
  Hit* getHit(unsigned int n) const;
  Cluster* getCluster(unsigned int n) const;

  std::vector<float>* getWaveform(std::string waveformName) const;
  float getWfSamplingWidth(std::string waveformName) const;
  void addWaveform(std::string waveformName, std::vector<float>* wf, float samplingWidth);

  void print();

  inline unsigned int getPlaneNum() const { return _planeNum; }
  inline unsigned int getNumHits() const { return _numHits; }
  inline unsigned int getNumClusters() const { return _numClusters; }
  inline std::map< std::string, std::vector<float>* > getWaveforms() const { return m_waveforms; }
  inline std::map< std::string, float > getWfSamplingWidths() const { return m_wfSamplingWidths; }

  friend class Event;
};

}

#endif // PLANE_H
