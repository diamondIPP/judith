#include "storageio.h"

#include <cassert>
#include <sstream>
#include <vector>
#include <iostream>

#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TBranch.h>
#include <TVectorD.h>

#include "event.h"
#include "track.h"
#include "plane.h"
#include "cluster.h"
#include "hit.h"

#ifndef VERBOSE
#define VERBOSE 1
#endif

using std::cout;
using std::endl;

namespace Storage {

void StorageIO::clearVariables()
{
  timeStamp = 0;
  frameNumber = 0;
  triggerOffset = 0;
  invalid = false;
  numHits = 0;
  numClusters = 0;
  numTracks = 0;
}

Event* StorageIO::readEvent(Long64_t n)
{
  /* Note: fill in reversed order: tracks first, hits last. This is so that
   * once a hit is produced, it can immediately recieve the address of its
   * parent cluster, likewise for clusters and track. */
  //cout << "read event" << endl;

  if (n >= _numEvents) throw "StorageIO: requested event outside range";

  if (_eventInfo &&_eventInfo->GetEntry(n) <= 0) throw "StorageIO: error reading event tree";
  if (_tracks && _tracks->GetEntry(n) <= 0) throw "StorageIO: error reading tracks tree";

  Event* event = new Event(_numPlanes);
  event->setTimeStamp(timeStamp);
  event->setFrameNumber(frameNumber);
  event->setTriggerOffset(triggerOffset);
  event->setTriggerInfo(triggerInfo);
  event->setInvalid(invalid);

  // Generate a list of track objects
  for (int ntrack = 0; ntrack < numTracks; ntrack++)
  {
    Track* track = event->newTrack();
    track->setOrigin(trackOriginX[ntrack], trackOriginY[ntrack]);
    track->setOriginErr(trackOriginErrX[ntrack], trackOriginErrY[ntrack]);
    track->setSlope(trackSlopeX[ntrack], trackSlopeY[ntrack]);
    track->setSlopeErr(trackSlopeErrX[ntrack], trackSlopeErrY[ntrack]);
    track->setCovariance(trackCovarianceX[ntrack], trackCovarianceY[ntrack]);
    track->setChi2(trackChi2[ntrack]);
  }

  for (unsigned int nplane = 0; nplane < _numPlanes; nplane++)
  {
    if (_hits.at(nplane) && _hits.at(nplane)->GetEntry(n) <= 0)
      throw "StorageIO: error reading hits tree";
    if (_clusters.at(nplane) && _clusters.at(nplane)->GetEntry(n) <= 0)
      throw "StorageIO: error reading clusters tree";

    // Insert waveforms into the plane
    // retrieve all of the branches from the tree
    if (m_waveformsTrees.size() == _numPlanes){
      if (m_waveformsTrees[nplane]){
        //cout << "plane: " << nplane << " waveforms" << endl;
        TObjArray* bArray = new TObjArray();
        bArray = m_waveformsTrees[nplane]->GetListOfBranches();
        int bSize = 0;
        // create a new std::vector<float> waveform and fill it
        std::vector<float>* temp_wf = new std::vector<float>;
        float temp_wf_sw = -1;
        std::string waveformName = "";
        // loop over all the branches
        for (Int_t i = 0; i < bArray->GetEntries(); i++ ) {
        // read the waveform branch for this plane and store its size
          bSize = ((TBranch*) bArray->At(i))->GetEntry(n);
          // read the name of the branch
          std::string bName( bArray->At(i)->GetName() );
          if (bSize < 10) {
            temp_wf_sw = WfSamplingWidths.at(bName);
waveformName = bName;
//cout << waveformName << " " << bName << endl;
            if (waveformName.compare("")) {
              event->getPlane(nplane)->addWaveform(waveformName, temp_wf, temp_wf_sw);
              //cout << "==== added " << waveformName << " with SW: " << temp_wf_sw << 
              //" size: " << temp_wf->size() << " wf pointer: " << temp_wf << endl;
            }
            else throw "error adding a waveform at storageio.cpp";
            waveformName = "";
            temp_wf = new std::vector<float>;
          }
          else {
            //cout << bName << " : " << bSize/sizeof(float) << endl;
            for (unsigned int i=0; i < bSize/sizeof(float); i++) {
              temp_wf->push_back( Waveforms.at(bName)[i] );
            }
            if (!waveformName.compare("")) waveformName = bName;
            else throw "error2 adding a waveform at storageio.cpp";
          }
        }   // loop over waveform branches
      }
    }

    // Generate the cluster objects
    for (int ncluster = 0; ncluster < numClusters; ncluster++)
    {
      Cluster* cluster = event->newCluster(nplane);
      cluster->setPix(clusterPixX[ncluster], clusterPixY[ncluster]);
      cluster->setPixErr(clusterPixErrX[ncluster], clusterPixErrY[ncluster]);
      cluster->setPos(clusterPosX[ncluster], clusterPosY[ncluster], clusterPosZ[ncluster]);
      cluster->setPosErr(clusterPosErrX[ncluster], clusterPosErrY[ncluster], clusterPosErrZ[ncluster]);

      // If this cluster is in a track, mark this (and the tracks tree is active)
      if (_tracks && clusterInTrack[ncluster] >= 0)
      {
        Track* track = event->getTrack(clusterInTrack[ncluster]);
        track->addCluster(cluster);
        cluster->setTrack(track);
      }
    }

    // Generate a list of all hit objects
    for (int nhit = 0; nhit < numHits; nhit++)
    {
      if (_noiseMasks && _noiseMasks->at(nplane)[hitPixX[nhit]][hitPixY[nhit]])
      {
        if (hitInCluster[nhit] >= 0)
          throw "StorageIO: tried to mask a hit which is already in a cluster";
        continue;
      }

      Hit* hit = event->newHit(nplane);
      hit->setPix(hitPixX[nhit], hitPixY[nhit]);
      hit->setPos(hitPosX[nhit], hitPosY[nhit], hitPosZ[nhit]);
      hit->setValue(hitValue[nhit]);
      hit->setTiming(hitTiming[nhit]);

      // If this hit is in a cluster, mark this (and the clusters tree is active)
      if (_clusters.at(nplane) && hitInCluster[nhit] >= 0)
      {
        Cluster* cluster = event->getCluster(hitInCluster[nhit]);
        cluster->addHit(hit);
      }
    }
  }

  return event;
}

void StorageIO::writeEvent(Event* event)
{
  if (_fileMode == INPUT) throw "StorageIO: can't write event in input mode";

  timeStamp = event->getTimeStamp();
  frameNumber = event->getFrameNumber();
  triggerOffset = event->getTriggerOffset();
  triggerInfo = event->getTriggerInfo();
  invalid = event->getInvalid();

  numTracks = event->getNumTracks();
  if (numTracks > MAX_TRACKS) throw "StorageIO: event exceeds MAX_TRACKS";

  // Set the object track values into the arrays for writing to the root file
  for (int ntrack = 0; ntrack < numTracks; ntrack++)
  {
    Track* track = event->getTrack(ntrack);
    trackOriginX[ntrack] = track->getOriginX();
    trackOriginY[ntrack] = track->getOriginY();
    trackOriginErrX[ntrack] = track->getOriginErrX();
    trackOriginErrY[ntrack] = track->getOriginErrY();
    trackSlopeX[ntrack] = track->getSlopeX();
    trackSlopeY[ntrack] = track->getSlopeY();
    trackSlopeErrX[ntrack] = track->getSlopeErrX();
    trackSlopeErrY[ntrack] = track->getSlopeErrY();
    trackCovarianceX[ntrack] = track->getCovarianceX();
    trackCovarianceY[ntrack] = track->getCovarianceY();
    trackChi2[ntrack] = track->getChi2();
  }

  for (unsigned int nplane = 0; nplane < _numPlanes; nplane++)
  {
    Plane* plane = event->getPlane(nplane);

    numClusters = plane->getNumClusters();
    if (numClusters > MAX_CLUSTERS) throw "StorageIO: event exceeds MAX_CLUSTERS";

    // Set the object cluster values into the arrays for writig into the root file
    for (int ncluster = 0; ncluster < numClusters; ncluster++)
    {
      Cluster* cluster = plane->getCluster(ncluster);
      clusterPixX[ncluster] = cluster->getPixX();
      clusterPixY[ncluster] = cluster->getPixY();
      clusterPixErrX[ncluster] = cluster->getPixErrX();
      clusterPixErrY[ncluster] = cluster->getPixErrY();
      clusterPosX[ncluster] = cluster->getPosX();
      clusterPosY[ncluster] = cluster->getPosY();
      clusterPosZ[ncluster] = cluster->getPosZ();
      clusterPosErrX[ncluster] = cluster->getPosErrX();
      clusterPosErrY[ncluster] = cluster->getPosErrY();
      clusterPosErrZ[ncluster] = cluster->getPosErrZ();
      clusterInTrack[ncluster] = cluster->getTrack() ? cluster->getTrack()->getIndex() : -1;
    }

    numHits = plane->getNumHits();
    if (numHits > MAX_HITS) throw "StorageIO: event exceeds MAX_HITS";

    // Set the object hit values into the arrays for writing into the root file
    for (int nhit = 0; nhit < numHits; nhit++)
    {
      Hit* hit = plane->getHit(nhit);
      hitPixX[nhit] = hit->getPixX();
      hitPixY[nhit] = hit->getPixY();
      hitPosX[nhit] = hit->getPosX();
      hitPosY[nhit] = hit->getPosY();
      hitPosZ[nhit] = hit->getPosZ();
      hitValue[nhit] = hit->getValue();
      hitTiming[nhit] = hit->getTiming();
      hitInCluster[nhit] = hit->getCluster() ? hit->getCluster()->getIndex() : -1;
    }
    
    // read all the waveforms from the plane
    std::map<std::string, std::vector<float>* > wfs = plane->getWaveforms();
    std::map<std::string, float > wfsSW = plane->getWfSamplingWidths();
    // loop over waveforms
    for (waveform_it iterator = wfs.begin(); iterator != wfs.end(); iterator++){
      //std::cout << iterator->first << " " << iterator->second->size() << std::endl;
      if(iterator->second->size() > MAX_WAVEFORM_POINTS) {
        throw "StorageO::WriteEvent: Number of waveform points exceeds MAX_WAVEFORM_POINTS";
      }
      else if(iterator->second->size() < 10)
        std::cout <<  "WARNING   StorageO::WriteEvent: Less then 10 points in a waveform" << std::endl;
      // store the sampling width
      WfSamplingWidths[iterator->first] = wfsSW.at(iterator->first);

      // create a key in Waveforms if it doesn't exist
      if(Waveforms.find(iterator->first)==Waveforms.end())
        Waveforms.insert( std::pair<std::string, float* > 
          ( iterator->first, new float[MAX_WAVEFORM_POINTS] ) );

      // copy plane->waveform content into the Waveform map
      float* Waveform_it = Waveforms.at(iterator->first);
      for(unsigned int i=0; i<iterator->second->size();i++)
        // first set the value of iter and then increase it to iter+1
        *Waveform_it++ = iterator->second->at(i); 

      // create a new waveform branch if it doesn't yet exist
      if( !m_waveformsTrees.empty() ) {
        if( !m_waveformsTrees.at(nplane)->GetBranch(iterator->first.c_str()) ) {
          // create title of the wf branch (i.e. branchName[xxx]/F)     
          std::stringstream branchType; branchType << iterator->first << 
          "[" << iterator->second->size() << "]/F";
          // create title of the wf sampling width branch (i.e. branchNameSW/F)     
          std::stringstream branchTypeSW; branchTypeSW << iterator->first << "SW/F";
          std::stringstream branchNameSW; branchNameSW << iterator->first << "SW";
          // create a new wf branch in the appropriate tree
          m_waveformsTrees.at(nplane)->Branch( iterator->first.c_str(), 
            Waveforms[iterator->first], branchType.str().c_str() );
          // create a new wf sampling width branch in the appropriate tree
          m_waveformsTrees.at(nplane)->Branch( branchNameSW.str().c_str(), 
            &WfSamplingWidths[iterator->first], branchTypeSW.str().c_str() );
        }
      }
    } // end loop over waveforms
 
    if (nplane >= _hits.size()) throw "StorageIO: event has too many planes for the storage";

    // Fill the plane by plane trees for this plane
    if (_hits.at(nplane)) _hits.at(nplane)->Fill();
    if (m_waveformsTrees.at(nplane)) m_waveformsTrees.at(nplane)->Fill();
    if (_clusters.at(nplane)) _clusters.at(nplane)->Fill();
  }

  // Write the track and event info here so that if any errors occured they won't be desynchronized
  if (_tracks) _tracks->Fill();
  if (_eventInfo) _eventInfo->Fill();

  _numEvents++;
}

void StorageIO::setNoiseMasks(std::vector<bool**>* noiseMasks)
{
  if (noiseMasks && _numPlanes != noiseMasks->size())
    throw "StorageIO: noise mask has more planes than will be read in";
  _noiseMasks = noiseMasks;
}

Long64_t StorageIO::getNumEvents() const
{
  assert(_fileMode != OUTPUT && "StorageIO: can't get number of entries in output mode");
  return _numEvents;
}

unsigned int StorageIO::getNumPlanes() const { return _numPlanes; }

Storage::Mode StorageIO::getMode() const { return _fileMode; }

StorageIO::StorageIO(const char* filePath, Mode fileMode, unsigned int numPlanes,
                     const unsigned int treeMask, const std::vector<bool>* planeMask, pointer_to_array BinWidth) :
  _filePath(filePath), _file(0), _fileMode(fileMode), _numPlanes(0), _numEvents(0),
  _noiseMasks(0)
{
  if      (fileMode == INPUT)  _file = new TFile(_filePath, "READ");
  else if (fileMode == OUTPUT) _file = new TFile(_filePath, "RECREATE");

  if (!_file) throw "StorageIO: file didn't initialize";

  _tracks = 0;
  _eventInfo = 0;

  timeStamp = 0;
  numHits = 0;
  numClusters = 0;
  numTracks = 0;

  // Plane mask holds a true for masked planes
  if (planeMask && fileMode == OUTPUT)
    throw "StorageIO: can't use a plane mask in output mode";

  // In output mode, create the directory structure and the relevant trees
  if (_fileMode == OUTPUT)
  {
    cout << "OUTPUT" << endl;
    if (planeMask && VERBOSE)
      cout << "WARNING :: StorageIO: disregarding plane mask in output mode";

    if (BinWidth) {
      // if width of bins in waveforms is specified,
      // write it to the file header.
      cout << BinWidth << endl;
      TVectorD binWidth0(1024);
      TVectorD binWidth1(1024);
      TVectorD binWidth2(1024);
      TVectorD binWidth3(1024);
      for(int i = 0; i < 1024; i++){
        binWidth0[i] = BinWidth[0][i];
        binWidth1[i] = BinWidth[1][i];
        binWidth2[i] = BinWidth[2][i];
        binWidth3[i] = BinWidth[3][i];
      }
      _file->cd();
      binWidth0.Write("binWidth0");
      binWidth1.Write("binWidth1");
      binWidth2.Write("binWidth2");
      binWidth3.Write("binWidth3");
    }

    _numPlanes = numPlanes;

    for (unsigned int nplane = 0; nplane < _numPlanes; nplane++)
    {
      std::stringstream ss;
      ss << "Plane" << nplane;

      TDirectory* dir = _file->mkdir(ss.str().c_str());
      dir->cd();

      TTree* hits = new TTree("Hits", "Hits");
      TTree* waveformsTreePl = new TTree("Waveforms", "Waveforms");
      TTree* clusters = new TTree("Clusters", "Clusters");

      _hits.push_back(hits);
      m_waveformsTrees.push_back(waveformsTreePl);
      _clusters.push_back(clusters);

      hits->Branch("NHits", &numHits, "NHits/I");
      hits->Branch("PixX", hitPixX, "HitPixX[NHits]/I");
      hits->Branch("PixY", hitPixY, "HitPixY[NHits]/I");
      hits->Branch("Value", hitValue, "HitValue[NHits]/I"); //Matevz 20141203 I to D
      hits->Branch("Timing", hitTiming, "HitTiming[NHits]/I");
      hits->Branch("HitInCluster", hitInCluster, "HitInCluster[NHits]/I");
//      hits->Branch("InCluster", hitInCluster, "HitInCluster[NHits]/I");
      hits->Branch("PosX", hitPosX, "HitPosX[NHits]/D");
      hits->Branch("PosY", hitPosY, "HitPosY[NHits]/D");
      hits->Branch("PosZ", hitPosZ, "HitPosZ[NHits]/D");

      clusters->Branch("NClusters", &numClusters, "NClusters/I");
      clusters->Branch("PixX", clusterPixX, "ClusterPixX[NClusters]/D");
      clusters->Branch("PixY", clusterPixY, "ClusterPixY[NClusters]/D");
      clusters->Branch("PixErrX", clusterPixErrX, "ClusterPixErrX[NClusters]/D");
      clusters->Branch("PixErrY", clusterPixErrY, "ClusterPixErrY[NClusters]/D");
      clusters->Branch("InTrack", clusterInTrack, "ClusterInTrack[NClusters]/I");
      clusters->Branch("PosX", clusterPosX, "ClusterPosX[NClusters]/D");
      clusters->Branch("PosY", clusterPosY, "ClusterPosY[NClusters]/D");
      clusters->Branch("PosZ", clusterPosZ, "ClusterPosZ[NClusters]/D");
      clusters->Branch("PosErrX", clusterPosErrX, "ClusterPosErrX[NClusters]/D");
      clusters->Branch("PosErrY", clusterPosErrY, "ClusterPosErrY[NClusters]/D");
      clusters->Branch("PosErrZ", clusterPosErrZ, "ClusterPosErrZ[NClusters]/D");
    }

    _file->cd();
    _tracks = new TTree("Tracks", "Track parameters");
    _eventInfo = new TTree("Event", "Event information");

    _eventInfo->Branch("TimeStamp", &timeStamp, "TimeStamp/l");
    _eventInfo->Branch("FrameNumber", &frameNumber, "FrameNumber/l");
    _eventInfo->Branch("TriggerOffset", &triggerOffset, "TriggerOffset/I");
    _eventInfo->Branch("TriggerInfo", &triggerInfo, "TriggerInfo/I");
    _eventInfo->Branch("Invalid", &invalid, "Invalid/O");

    _tracks->Branch("NTracks", &numTracks, "NTracks/I");
    _tracks->Branch("SlopeX", trackSlopeX, "TrackSlopeX[NTracks]/D");
    _tracks->Branch("SlopeY", trackSlopeY, "TrackSlopeY[NTracks]/D");
    _tracks->Branch("SlopeErrX", trackSlopeErrX, "TrackSlopeErrX[NTracks]/D");
    _tracks->Branch("SlopeErrY", trackSlopeErrY, "TrackSlopeErrY[NTracks]/D");
    _tracks->Branch("OriginX", trackOriginX, "TrackOriginX[NTracks]/D");
    _tracks->Branch("OriginY", trackOriginY, "TrackOriginY[NTracks]/D");
    _tracks->Branch("OriginErrX", trackOriginErrX, "TrackOriginErrX[NTracks]/D");
    _tracks->Branch("OriginErrY", trackOriginErrY, "TrackOriginErrY[NTracks]/D");
    _tracks->Branch("CovarianceX", trackCovarianceX, "TrackCovarianceX[NTracks]/D");
    _tracks->Branch("CovarianceY", trackCovarianceY, "TrackCovarianceY[NTracks]/D");
    _tracks->Branch("Chi2", trackChi2, "TrackChi2[NTracks]/D");
  }

  // In input mode,
  if (_fileMode == INPUT)
  {
    cout << "INPUT" << endl;
    if (_numPlanes && VERBOSE)
      cout << "WARNING :: StorageIO: disregarding specified number of planes" << endl;

    if(BinWidth && _file->GetKey("binWidth0") && _file->GetKey("binWidth1") && _file->GetKey("binWidth2") && _file->GetKey("binWidth3") ) {
      cout << "TEST" << endl;
      TVectorD* binWidth0 = (TVectorD*) _file->Get("binWidth0");
      TVectorD* binWidth1 = (TVectorD*) _file->Get("binWidth1");
      TVectorD* binWidth2 = (TVectorD*) _file->Get("binWidth2");
      TVectorD* binWidth3 = (TVectorD*) _file->Get("binWidth3");
      cout << "TEST2" << endl;
      for(int i = 0; i < 1024; i++){
        BinWidth[0][i] = (*binWidth0)[i];
        BinWidth[1][i] = (*binWidth1)[i];
        BinWidth[2][i] = (*binWidth2)[i];
        BinWidth[3][i] = (*binWidth3)[i];
      }     
    }

    _numPlanes = 0; // Determine num planes from file structure

    unsigned int planeCount = 0;
    while (true)
    {
      std::stringstream ss;
      ss << "Plane" << planeCount;

      // Try to get this plane's directory
      TDirectory* dir = 0;
      _file->GetObject(ss.str().c_str(), dir);
      if (!dir) break;

      planeCount++;

      if (planeMask && planeCount > planeMask->size())
        throw "StorageIO: plane mask is too small";

      if (planeMask && planeMask->at(planeCount - 1)) continue;

      TTree* hits;
      TTree* waveforms;
      TTree* clusters;

      _file->GetObject(ss.str().append("/Hits").c_str(), hits);
      _file->GetObject(ss.str().append("/Waveforms").c_str(), waveforms);
      _file->GetObject(ss.str().append("/Clusters").c_str(), clusters);

      _hits.push_back(hits);
      //m_waveformsTrees.push_back(waveforms);
      _clusters.push_back(clusters);
      _numPlanes++;

      if (hits)
      {
        hits->SetBranchAddress("NHits", &numHits, &bNumHits);
        hits->SetBranchAddress("PixX", hitPixX, &bHitPixX);
        hits->SetBranchAddress("PixY", hitPixY, &bHitPixY);
        hits->SetBranchAddress("Value", hitValue, &bHitValue);
        hits->SetBranchAddress("Timing", hitTiming, &bHitTiming);
        hits->SetBranchAddress("HitInCluster", hitInCluster, &bHitInCluster);
        hits->SetBranchAddress("PosX", hitPosX, &bHitPosX);
        hits->SetBranchAddress("PosY", hitPosY, &bHitPosY);
        hits->SetBranchAddress("PosZ", hitPosZ, &bHitPosZ);
      }
      
      if (waveforms) {   //INPUT mode is not implemented for waveforms in this 
        cout << "waveforms" << endl;
        //if (0) {
        // Add this tree to the current plane
        m_waveformsTrees.push_back(waveforms);
        // retrieve all of the branches from the tree
        TObjArray* bArray = waveforms->GetListOfBranches();
        // loop over all the branches
        for (Int_t i = 0; i < bArray->GetEntries(); i++ ) {
          // retreive the name of the waveform branch
          const char* bName = bArray->At(i)->GetName();
          const float bSize = ( (TBranch*)(bArray->At(i)) )->GetTotalSize()/
          ( (TBranch*)(bArray->At(i)) )->GetEntries();
          cout << "waveform name: " << bName << " , size: " << bSize << endl;
          if (bSize > 10) {
            // create a new float* key in Waveforms map
            Waveforms.insert( std::pair<std::string, float* > 
              (std::string(bName), new float[MAX_WAVEFORM_POINTS] ) );
            // set address of the branch
            waveforms->SetBranchAddress(bName, Waveforms.at( std::string(bName) ) );
          }
          else {
            // create a new float* key in Waveforms map
            WfSamplingWidths.insert( std::pair<std::string, float > 
              (std::string(bName), float(0) ) );
            // set address of the branch
            waveforms->SetBranchAddress(bName, &WfSamplingWidths.at( std::string(bName) ) );
          }
        }   // loop over waveform branches
      }

      if (clusters)
      {
        clusters->SetBranchAddress("NClusters", &numClusters, &bNumClusters);
        clusters->SetBranchAddress("PixX", clusterPixX, &bClusterPixX);
        clusters->SetBranchAddress("PixY", clusterPixY, &bClusterPixY);
        clusters->SetBranchAddress("PixErrX", clusterPixErrX, &bClusterPixErrX);
        clusters->SetBranchAddress("PixErrY", clusterPixErrY, &bClusterPixErrY);
        clusters->SetBranchAddress("InTrack", clusterInTrack, &bClusterInTrack);
        clusters->SetBranchAddress("PosX", clusterPosX, &bClusterPosX);
        clusters->SetBranchAddress("PosY", clusterPosY, &bClusterPosY);
        clusters->SetBranchAddress("PosZ", clusterPosZ, &bClusterPosZ);
        clusters->SetBranchAddress("PosErrX", clusterPosErrX, &bClusterPosErrX);
        clusters->SetBranchAddress("PosErrY", clusterPosErrY, &bClusterPosErrY);
        clusters->SetBranchAddress("PosErrZ", clusterPosErrZ, &bClusterPosErrZ);
      }
    }

    _file->GetObject("Tracks", _tracks);
    _file->GetObject("Event", _eventInfo);

    if (_eventInfo)
    {
      _eventInfo->SetBranchAddress("TimeStamp", &timeStamp, &bTimeStamp);
      _eventInfo->SetBranchAddress("FrameNumber", &frameNumber, &bFrameNumber);
      _eventInfo->SetBranchAddress("TriggerOffset", &triggerOffset, &bTriggerOffset);
      _eventInfo->SetBranchAddress("TriggerInfo", &triggerInfo, &bTriggerInfo);
      _eventInfo->SetBranchAddress("Invalid", &invalid, &bInvalid);
    }

    if (_tracks)
    {
      _tracks->SetBranchAddress("NTracks", &numTracks, &bNumTracks);
      _tracks->SetBranchAddress("SlopeX", trackSlopeX, &bTrackSlopeX);
      _tracks->SetBranchAddress("SlopeY", trackSlopeY, &bTrackSlopeY);
      _tracks->SetBranchAddress("SlopeErrX", trackSlopeErrX, &bTrackSlopeErrX);
      _tracks->SetBranchAddress("SlopeErrY", trackSlopeErrY, &bTrackSlopeErrY);
      _tracks->SetBranchAddress("OriginX", trackOriginX, &bTrackOriginX);
      _tracks->SetBranchAddress("OriginY", trackOriginY, &bTrackOriginY);
      _tracks->SetBranchAddress("OriginErrX", trackOriginErrX, &bTrackOriginErrX);
      _tracks->SetBranchAddress("OriginErrY", trackOriginErrY, &bTrackOriginErrY);
      _tracks->SetBranchAddress("CovarianceX", trackCovarianceX, &bTrackCovarianceX);
      _tracks->SetBranchAddress("CovarianceY", trackCovarianceY, &bTrackCovarianceY);
      _tracks->SetBranchAddress("Chi2", trackChi2, &bTrackChi2);
    }
  }

  if (_numPlanes < 1) throw "StorageIO: didn't initialize any planes";

  // Delete trees as per the tree flags
  if (treeMask)
  {
    for (unsigned int nplane = 0; nplane < _numPlanes; nplane++)
    {
      if (treeMask & Flags::HITS) { delete _hits.at(nplane); _hits.at(nplane) = 0; }
      if (treeMask & Flags::CLUSTERS) { delete _clusters.at(nplane); _clusters.at(nplane) = 0; }
    }
    if (treeMask & Flags::TRACKS) { delete _tracks; _tracks = 0; }
    if (treeMask & Flags::EVENTINFO) { delete _eventInfo; _eventInfo = 0; }
  }

  assert(_hits.size() == _clusters.size() && "StorageIO: varying number of planes");

  _numEvents = 0;

  if (_fileMode == INPUT)
  {
    Long64_t nEventInfo = (_eventInfo) ? _eventInfo->GetEntriesFast() : 0;
    Long64_t nTracks = (_tracks) ? _tracks->GetEntriesFast() : 0;
    Long64_t nHits = 0;
    Long64_t nClusters = 0;
    for (unsigned int nplane = 0; nplane < _numPlanes; nplane++)
    {
      if (_hits.at(nplane)) nHits += _hits.at(nplane)->GetEntriesFast();
      if (_clusters.at(nplane)) nClusters += _clusters.at(nplane)->GetEntriesFast();
    }

    if (nHits % _numPlanes || nClusters % _numPlanes)
      throw "StorageIO: number of events in different planes mismatch";

    nHits /= _numPlanes;
    nClusters /= _numPlanes;

    if (!_numEvents && nEventInfo) _numEvents = nEventInfo;
    if (!_numEvents && nTracks) _numEvents = nTracks;
    if (!_numEvents && nHits) _numEvents = nHits;
    if (!_numEvents && nClusters) _numEvents = nClusters;

    if ((nEventInfo && _numEvents != nEventInfo) ||
        (nTracks && _numEvents != nTracks) ||
        (nHits && _numEvents != nHits) ||
        (nClusters && _numEvents != nClusters))
      throw "StorageIO: all trees don't have the same number of events";
  }
  cout << "end IO" << endl;
}

StorageIO::~StorageIO()
{
  if (_file && _fileMode == OUTPUT)
  {
    _file->Write();
    delete _file;
  }
}

}
