// R__LOAD_LIBRARY($ROOTSYS/test/libEvent.so)

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <TROOT.h>
#include <TVector3.h>
//#include <TRandom3.h>

#include "TPC3DCell.h"
// #include "TPC3DCellID.h"



// struct to be filled with hits
struct hit {
  double x, y, z, t;
  double xReco, yReco, zReco, tReco;
  int trackID, pdg;
  Bool_t registered;
  UChar_t eventID;
  UShort_t eventIndex;
  TVector3 recoPos;

  TPC3DCellID cellID;

  double getRecoTime() const { return tReco; }
  TPC3DCellID get3DCellID() const { return cellID; }
  TVector3 getRecoPosition() const { return recoPos; }

  void setUseHitInTracking(const bool useInTracking) { registered = useInTracking; }
};



// int readfile(double rate, std::string infilename, std::string outfilename, int evtID, std::vector<hit> &hits, std::map<std::pair<int, int>, std::vector<hit*>> &HitMap, int origindex = -1);
// void checkForInvalidDigits(std::map<std::pair<int, int>, std::vector<hit*>> HitMap);
int readfile(double rate, std::string infilename, std::string outfilename, int evtID, std::vector<hit>& hits,
             std::map<std::pair<int, int>, std::vector<unsigned long long int>>& HitMap, int origindex = -1);
void checkForInvalidDigits(std::vector<hit>& hits, std::map<std::pair<int, int>, std::vector<unsigned long long int>> HitMap);

void cellcreator(std::vector<hit>& hits, std::map<TPC3DCellID, TPC3DCell>& m_CellsWithHits, std::vector<TPC3DCell>& m_TPC3DCells);
void simpleTPCBackgroundRejection(std::vector<TPC3DCell>& cells, std::map<TPC3DCellID, TPC3DCell*> CellMap, std::vector<TPC3DCell*> ConnectedCells,
                                  std::vector<hit>& hits);
void isolatedCellFinder(TPC3DCell* currentCell, std::map<TPC3DCellID, TPC3DCell*>& CellMap, std::vector<TPC3DCell*>& ConnectedCells);
void connectedZCellFinder(TPC3DCell* currentCell, std::map<TPC3DCellID, TPC3DCell*>& CellMap, std::vector<TPC3DCell*>& ConnectedCells);
void invalidateConnectedCells(std::vector<TPC3DCell*>& ConnectedCells, std::vector<hit>& hits);
double calculateMean(const std::vector<double>& values);

TRandom3* myRandom = new TRandom3(0);
UInt_t seed = myRandom->GetSeed();

// Parameters for deadtime implementation
int m_readoutCellSizeXY = 2; // unit cm
double m_readoutPixelPitch = 50.0 / 10000; // unit cm; numerator is pixel pitch in um 
int m_referencePositionXY = -128; // unit cm
double m_driftVelocity = 77.0 / 10000; // unit cm/ns
int m_pixelIntegrationTime = 25; // unit ns
int m_pixelDeadTime = 750; // unit ns

// Parameters for cell creation (see TPC3DCellCreatorModule)
double m_backwardEndplate = -83.12; // unit cm 
double m_readoutCellSizeZT = 250.0;
uint m_occupancyCut = 2000;

// Parameters for simple bkg rejection (see SimpleTPCBackgroundRejectionModule)
bool m_useNeighbourFinder = false;
uint m_MaxConnectedCellsInZ = 5;
uint m_MinConnectedCellsForIsolationCut = 6;


// Index of bkg file to be used
UInt_t bkg_rn = myRandom->Integer(25) + 10000;

// UInt_t bkg_rn = 10000;

void create_overlay()
{
  std::vector<hit> TPCHits;
  // std::map<std::pair<int, int>, std::vector<hit*>> m_currentHitMap;
  std::map<std::pair<int, int>, std::vector<unsigned long long int>> m_currentHitMap;

  // Dummy hitmaps for the beam background, because else the hitmap becomes too large. 
  // Only used for applying effects of pixel dead time which is not a mojor source of hit loss
  std::map<std::pair<int, int>, std::vector<unsigned long long int>> m_currentHitMap1;
  std::map<std::pair<int, int>, std::vector<unsigned long long int>> m_currentHitMap2;
  std::map<std::pair<int, int>, std::vector<unsigned long long int>> m_currentHitMap3;
  std::map<std::pair<int, int>, std::vector<unsigned long long int>> m_currentHitMap4;
  std::map<std::pair<int, int>, std::vector<unsigned long long int>> m_currentHitMap5;
  std::map<std::pair<int, int>, std::vector<unsigned long long int>> m_currentHitMap6;

  std::cout << "SEED = " << seed << std::endl;

  std::string outfile = "/gpfs/group/belle2/users/loeschca/ElectricField/dataSimpleRejection/pitch50/TPCEventOverlay" + std::to_string(
                        seed) + ".root";
  //std::string outfile = "./dataTPC/NewStruct" + std::to_string(seed) + ".root";
  const char* coutfilename = outfile.c_str();

  // factor for easily changing the amount of beam background in overlay. Has to be an integer
  int bkgscaling = 1;

  int original_index = readfile(1.00, "/gpfs/group/belle2/users/loeschca/dataUpsilon4Ss.root",  outfile,  1, TPCHits, m_currentHitMap,
                                -2);
  /*readfile(1.00, "/gpfs/group/belle2/users/loeschca/dataUpsilon4Ss.root",  outfile,  1, TPCHits, m_currentHitMap,
                                -3);
  readfile(1.00, "/gpfs/group/belle2/users/loeschca/dataUpsilon4Ss.root",  outfile,  1, TPCHits, m_currentHitMap,
                                -4);*/

  
  readfile(0.11, "/gpfs/group/belle2/users/loeschca/dataUpsilon4Ss.root",  outfile,  1, TPCHits, m_currentHitMap, original_index);
  readfile(28.8, "/gpfs/group/belle2/users/loeschca/dataBhabhas.root",     outfile,  2, TPCHits, m_currentHitMap);
  readfile(0.48, "/gpfs/group/belle2/users/loeschca/dataGammaGammas.root", outfile,  3, TPCHits, m_currentHitMap);
  readfile(0.11, "/gpfs/group/belle2/users/loeschca/dataMuMus.root",       outfile,  4, TPCHits, m_currentHitMap);
  readfile(0.09, "/gpfs/group/belle2/users/loeschca/dataTauTaus.root",     outfile,  5, TPCHits, m_currentHitMap);
  readfile(0.15, "/gpfs/group/belle2/users/loeschca/dataUUbars.root",      outfile,  6, TPCHits, m_currentHitMap);
  readfile(0.04, "/gpfs/group/belle2/users/loeschca/dataDDbars.root",      outfile,  7, TPCHits, m_currentHitMap);
  readfile(0.04, "/gpfs/group/belle2/users/loeschca/dataSSbars.root",      outfile,  8, TPCHits, m_currentHitMap);
  readfile(0.12, "/gpfs/group/belle2/users/loeschca/dataCCbars.root",      outfile,  9, TPCHits, m_currentHitMap);
  readfile(3.81, "/gpfs/group/belle2/users/loeschca/dataEEEEs.root",       outfile, 10, TPCHits, m_currentHitMap);
  readfile(1.81, "/gpfs/group/belle2/users/loeschca/dataEEMuMus.root",     outfile, 11, TPCHits, m_currentHitMap);
  // Bkg file rates are for a window of 60 us (instead of 120 us) since they are placed only in the central region 
  readfile(bkgscaling * 147.89, "/gpfs/group/belle2/users/loeschca/ElectricField/dataBkgField5/Brems_HER_study-phase3-BG19-" +
           std::to_string(bkg_rn) + ".root",
           outfile, 12, TPCHits, m_currentHitMap1);
  readfile(bkgscaling * 497.57, "/gpfs/group/belle2/users/loeschca/ElectricField/dataBkgField5/Brems_LER_study-phase3-BG19-" +
           std::to_string(bkg_rn) + ".root",
           outfile, 13, TPCHits, m_currentHitMap2);
  readfile(bkgscaling * 980.10, "/gpfs/group/belle2/users/loeschca/ElectricField/dataBkgField5/Coulomb_HER_study-phase3-BG19-" +
           std::to_string(bkg_rn) + ".root",
           outfile, 14, TPCHits, m_currentHitMap3);

  // Coulomb LER from 5 separate files for x5 bkg studies, since one file is to large to create
  for (int b = 0; b < bkgscaling; b++) {
    UInt_t bkg_rn_CLER = myRandom->Integer(25) + 10000;
    readfile(11498.46, "/gpfs/group/belle2/users/loeschca/ElectricField/dataBkgField/Coulomb_LER_study-phase3-BG19-" + std::to_string(
               bkg_rn_CLER) + ".root",  outfile, 15, TPCHits, m_currentHitMap4);
  }

  

  /*
  readfile(bkgscaling * 11498.46, "/gpfs/group/belle2/users/loeschca/ElectricField/dataBkgField5/Coulomb_LER_study-phase3-BG19-" + std::to_string(
             bkg_rn2) + ".root",  outfile, 15, TPCHits, m_currentHitMap);
  readfile(bkgscaling * 11498.46, "/gpfs/group/belle2/users/loeschca/ElectricField/dataBkgField5/Coulomb_LER_study-phase3-BG19-" + std::to_string(
             bkg_rn3) + ".root",  outfile, 15, TPCHits, m_currentHitMap);
  readfile(bkgscaling * 11498.46, "/gpfs/group/belle2/users/loeschca/ElectricField/dataBkgField5/Coulomb_LER_study-phase3-BG19-" + std::to_string(
             bkg_rn4) + ".root",  outfile, 15, TPCHits, m_currentHitMap);
  readfile(bkgscaling * 11498.46, "/gpfs/group/belle2/users/loeschca/ElectricField/dataBkgField5/Coulomb_LER_study-phase3-BG19-" + std::to_string(
             bkg_rn5) + ".root",  outfile, 15, TPCHits, m_currentHitMap);
  */

  
  readfile(bkgscaling * 14.54, "/gpfs/group/belle2/users/loeschca/ElectricField/dataBkgField5/Touschek_HER_study-phase3-BG19-" +
           std::to_string(bkg_rn) + ".root",
           outfile, 16, TPCHits, m_currentHitMap5);
  readfile(bkgscaling * 7164.54,  "/gpfs/group/belle2/users/loeschca/ElectricField/dataBkgField5/Touschek_LER_study-phase3-BG19-" +
           std::to_string(
             bkg_rn) + ".root", outfile, 17, TPCHits, m_currentHitMap6);
  


  // std::cout << TPCHits->x.at(TPCHits->x.size()-1) << std::endl;
  std::cout << "Size of HitMap = " << m_currentHitMap.size() << std::endl;
   
  checkForInvalidDigits(TPCHits, m_currentHitMap);


  std::map<TPC3DCellID, TPC3DCell> m_CellsWithHits;
  std::vector<TPC3DCell> m_TPC3DCells;
  cellcreator(TPCHits, m_CellsWithHits, m_TPC3DCells);

  std::map<TPC3DCellID, TPC3DCell*> m_cellMap;
  std::vector<TPC3DCell*> m_connectedCells;
  simpleTPCBackgroundRejection(m_TPC3DCells, m_cellMap, m_connectedCells, TPCHits);
  



  double entries = TPCHits.size();

  double x, y, z, t;
  double xReco, yReco, zReco, tReco;
  int trackID, pdg;
  Bool_t registered;
  UChar_t eventID;
  UShort_t eventIndex;

  TFile* file = new TFile(coutfilename, "RECREATE");
  TTree* tree = new TTree("tree", "tree");

  tree->Branch("x", &x, "x/D");
  tree->Branch("y", &y, "y/D");
  tree->Branch("z", &z, "z/D");
  tree->Branch("t", &t, "t/D");

  tree->Branch("xReco", &xReco, "xReco/D");
  tree->Branch("yReco", &yReco, "yReco/D");
  tree->Branch("zReco", &zReco, "zReco/D");
  tree->Branch("tReco", &tReco, "tReco/D");

  tree->Branch("eventID", &eventID, "eventID/b");

  tree->Branch("eventIndex", &eventIndex, "eventIndex/s");

  tree->Branch("trackID", &trackID, "trackID/I");

  tree->Branch("pdg", &pdg, "pdg/I");

  tree->Branch("registered", &registered, "registered/O");


  for (int i = 0; i < entries; i++) {
    x = TPCHits.at(i).x;
    y = TPCHits.at(i).y;
    z = TPCHits.at(i).z;
    t = TPCHits.at(i).t;

    xReco = TPCHits.at(i).xReco;
    yReco = TPCHits.at(i).yReco;
    zReco = TPCHits.at(i).zReco;
    tReco = TPCHits.at(i).tReco;

    eventID = TPCHits.at(i).eventID;
    eventIndex = TPCHits.at(i).eventIndex;
    trackID = TPCHits.at(i).trackID;
    pdg = TPCHits.at(i).pdg;

    registered = TPCHits.at(i).registered;

    tree->Fill();
  }



  tree->Write();
  file->Close();

  std::cout << " "  << std::endl;
  std::cout << "File created = " << outfile << std::endl;


  // readfile(2.5*7164.54,  "/gpfs/group/belle2/users/loeschca/ElectricField/dataBkgField/Touschek_LER_study-phase3-BG19-"+std::to_string(bkg_add)+".root", outfile, 17);
  // readfile(2.5*11498.46, "/gpfs/group/belle2/users/loeschca/ElectricField/dataBkgField/Coulomb_LER_study-phase3-BG19-"+std::to_string(bkg_add)+".root",  outfile, 15);

}

//int readfile(double rate, std::string infilename, std::string outfilename, int evtID, std::vector<hit> &hits, std::map<std::pair<int, int>, std::vector<hit*>> &HitMap, int origindex = -1)
int readfile(double rate, std::string infilename, std::string outfilename, int evtID, std::vector<hit>& hits,
             std::map<std::pair<int, int>, std::vector<unsigned long long int>>& HitMap, int origindex = -1)
{
  // Number of events happening in 60 us
  unsigned int eventnumber;
  if (origindex <= -2) {
    eventnumber = 1;
  } else {
    eventnumber = myRandom->Poisson(rate);
  }


  std::cout << "eventnumber = " << eventnumber << " ; evtID = " << evtID << std::endl;

  // Conversion of string to char for root
  const char* cinfilename = infilename.c_str();
  TFile* dataFile = TFile::Open(cinfilename);

  // Creating a TTreeReader
  TTreeReader reader("tree", dataFile);
  //Read the TPCSimHits in the tree entry:
  TTreeReaderArray<double> simHitx(reader, "TPCSimHits.m_x");
  TTreeReaderArray<double> simHity(reader, "TPCSimHits.m_y");
  TTreeReaderArray<double> simHitz(reader, "TPCSimHits.m_z");
  TTreeReaderArray<int> simHitTrackID(reader, "TPCSimHits.m_trackID");
  TTreeReaderArray<int> simHitPDG(reader, "TPCSimHits.m_pdg");

  // Number of events in file
  unsigned int entries = reader.GetEntries();

  // Conversion of string to char for root
  const char* coutfilename = outfilename.c_str();

  // Varaibles for calculating with coordinates
  hit thishit;
  double x, y, z;
  double xReco, yReco, zReco;
  int trackID, pdg;
  Bool_t registered = true;

  UChar_t eventID = evtID;

  UShort_t eventIndex;

  TPC3DCellID cellID;
  TVector3 recoPos;


  // Transverse diffusion coefficient in units cm/sqrt(cm)
  double transdiff = 84.0 / 10000.0;
  // Longitudinal diffusion coefficient in unis cm/sqrt(cm)
  double longdiff = 200.0 / 10000.0;

  // Create unique indices for input file
  unsigned int index;
  std::vector<int> indexvector;
  unsigned int vectorsize = 0;
  unsigned int failsave = 0;
  while (vectorsize < eventnumber) {
    bool unique = true;
    index = myRandom->Integer(entries);
    for (int l = 0; l < vectorsize; ++l) {
      // check for duplicates already created. Origindex is the index of the triggered Y(4S)
      if ((index == indexvector[l]) || index == origindex) {
        unique = false;
      }
    }

    if (failsave > 100) {
      unique = true;
      failsave = 0;
      std::cout << "Failsave catching! EventID = " << evtID << std::endl;
    }

    if (unique) {
      indexvector.push_back(index);
      ++vectorsize;
    } else {
      ++failsave;
    }
    // std::cout << "vectorsize = " << vectorsize << std::endl;
  }

  // Main loop over random indices in which the desired data are calculated
  for (int dice = 0; dice < eventnumber; ++dice) {
    // Getting a random event from the input file
    index = indexvector[dice];
    eventIndex = index;
    reader.SetEntry(index);

    // Setting a random time with respect to t0 (in ns)
    // 7454 is the number of bunch crossings per TPC volume
    int ntimestamp;
    if (origindex == -2) {
      ntimestamp = 0;
    } else if (evtID >= 12) {  // Bkg files are only placed in the central region where overlap with signal event can happen
      ntimestamp = myRandom->Integer(14908) - 7454;
    } else {
      ntimestamp = myRandom->Integer(29816) - 14908;
    }
    
    // Experimental manual event placing
    // if (origindex == -3) ntimestamp = 7454;
    // if (origindex == -4) ntimestamp = -7454;

    double timestamp = ntimestamp * 4.024;
    // Convert time to z coordinate of earliest hits (in cm)
    double zstamp = timestamp * 78.9 / 10000.0;

    // Get number of SimHits for loop and branchsize
    int nsimhits = simHitx.GetSize();

    for (int i_x = 0; i_x < nsimhits; ++i_x) {
      x = simHitx[i_x];
      y = simHity[i_x];
      z = simHitz[i_x] + 83.12;

      trackID = simHitTrackID[i_x];
      pdg = simHitPDG[i_x];

      // std::cout << z << std::endl;

      x += myRandom->Gaus(0.0, transdiff * sqrt(z / 2.0));
      y += myRandom->Gaus(0.0, transdiff * sqrt(z / 2.0));

      // Not saving digits which would diffuse outside of TPC volume
      double radius = sqrt(x * x + y * y);
      if (radius < 44.85 || radius > 109.4) continue;

      z += myRandom->Gaus(0.0, longdiff * sqrt(z));
      z += zstamp;
      double driftTime = z / m_driftVelocity;

      int xCellID = (x - m_referencePositionXY) / m_readoutCellSizeXY;
      int yCellID = (y - m_referencePositionXY) / m_readoutCellSizeXY;

      int xPixel = (x - m_referencePositionXY - m_readoutCellSizeXY * xCellID) / m_readoutPixelPitch;
      int yPixel = (y - m_referencePositionXY - m_readoutCellSizeXY * yCellID) / m_readoutPixelPitch;

      xReco = m_referencePositionXY + m_readoutCellSizeXY * xCellID + m_readoutPixelPitch * xPixel + m_readoutPixelPitch / 2.;
      yReco = m_referencePositionXY + m_readoutCellSizeXY * yCellID + m_readoutPixelPitch * yPixel + m_readoutPixelPitch / 2.;

      // if (origindex == -2) std::cout << "xReco = " << xReco << " ; yReco = " << yReco << std::endl;

      double recoDriftTime = (int)(driftTime / m_pixelIntegrationTime) * m_pixelIntegrationTime + m_pixelIntegrationTime / 2.;
      zReco = recoDriftTime * m_driftVelocity;
      int zCellID = (zReco + 242.69) / (m_driftVelocity * m_readoutCellSizeZT); // 242.69 is length of TPC, add to avoid negative cellIDs from previous TPC volumes

      cellID.setXCell(xCellID);
      cellID.setYCell(yCellID);
      cellID.setZCell(zCellID);

      recoPos.SetXYZ(xReco, yReco, zReco);

      // Fill coordinates to hit
      thishit.x = (x);
      thishit.y = (y);
      thishit.z = (z);
      thishit.t = (driftTime);
      thishit.xReco = (xReco);
      thishit.yReco = (yReco);
      thishit.zReco = (zReco);
      thishit.tReco = (recoDriftTime);
      thishit.trackID = (trackID);
      thishit.pdg = (pdg);
      thishit.registered = (registered);
      thishit.eventID = (eventID);
      thishit.eventIndex = (eventIndex);
      thishit.cellID = (cellID);
      thishit.recoPos = (recoPos);

      hits.push_back(thishit);
      hit* lasthit = &hits.back();


      
      // Fill the hit map (memory problems with x5 Bkg)
      // Reco position in nm to have integers and avoid rounding and precision issues as much as possible
      const int xRecoMap = std::round((m_referencePositionXY + m_readoutCellSizeXY * xCellID + m_readoutPixelPitch * xPixel +
                                       m_readoutPixelPitch / 2.) * 10000000);
      const int yRecoMap = std::round((m_referencePositionXY + m_readoutCellSizeXY * yCellID + m_readoutPixelPitch * yPixel +
                                       m_readoutPixelPitch / 2.) * 10000000);

      auto hitInMap = HitMap.find({xRecoMap, yRecoMap});
      // Reco position xReco, yReco already in the hit map, so append the current digit to the corresponding vector
      if (hitInMap != HitMap.end()) {
        // hitInMap->second.push_back(lasthit);
        hitInMap->second.push_back(hits.size() - 1);
        // std::cout << "New Pixel" << std::endl;
      } else {
        // Reco position xReco, yReco not yet in the hit map, so add it together with the time
        // std::vector<hit*> hitVector = {lasthit};
        std::vector<unsigned long long int> hitVector = {hits.size() - 1};
        HitMap.insert(std::make_pair(std::make_pair(xRecoMap, xRecoMap), hitVector));
        // std::cout << "Same Pixel!!!!!!!" << std::endl;
      }
      
      // branchsize increment here and not in outer loop with branchsize+=nsimhits in order to not count hits which diffuse outside of tpc volume
      // branchsize++;
    }

  }

  if (origindex == -2) return indexvector[0];
  return 0;

}


//void checkForInvalidDigits(std::map<std::pair<int, int>, std::vector<hit*>> HitMap)
void checkForInvalidDigits(std::vector<hit>& hits, std::map<std::pair<int, int>, std::vector<unsigned long long int>> HitMap)
{
  for (/*const*/ auto& hitInMap : HitMap) {
    // only one hit for the given pixel (= reco position xReco, yReco), nothing to do
    if (hitInMap.second.size() == 1) continue;

    // several hit on the same pixel, sort by time first

    std::vector<hit*> hitsOfPixel;
    for (unsigned int i = 0; i < hitInMap.second.size(); i++) {
      hitsOfPixel.push_back(&hits.at(hitInMap.second.at(i)));
    }
    // std::vector<hit*>& hitsOfPixel = &pixelhits; // = hitInMap.second;

    std::sort(hitsOfPixel.begin(), hitsOfPixel.end(), [](const hit * lhs, const hit * rhs)
    { return lhs->getRecoTime() < rhs->getRecoTime(); });

    // for each hit on a given pixel check the times of all subsequent hits (in time) and set them invalid in case they are within the dead time of the pixel
    std::vector<hit*>::iterator currentHit;
    for (currentHit = hitsOfPixel.begin(); currentHit != hitsOfPixel.end() - 1; ++currentHit) {
      for (auto nextHit = currentHit + 1; nextHit != hitsOfPixel.end(); ++nextHit) {
        if (abs((*nextHit)->getRecoTime() - (*currentHit)->getRecoTime()) < m_pixelDeadTime) {
          (*nextHit)->setUseHitInTracking(false);
          //std::cout << " t of hit = " << (*nextHit)->getRecoTime() << std::endl;
          // std::cout << " Used in tracking set to false " << std::endl;
        }
      }
    }
  }
}


void cellcreator(std::vector<hit>& hits, std::map<TPC3DCellID, TPC3DCell>& m_CellsWithHits, std::vector<TPC3DCell>& m_TPC3DCells)
{
  std::vector<double> xCoords, yCoords, zCoords;
  // clear all containers that are filled each event
  m_CellsWithHits.clear();

  uint digitindex = 0;
  // Loop over all digits
  for (auto& tpchit : hits) {
    auto cellPresent = m_CellsWithHits.find(tpchit.get3DCellID());
    if (cellPresent == m_CellsWithHits.end()) {
      TPC3DCell newCell(tpchit.get3DCellID(), m_readoutCellSizeXY, m_readoutCellSizeZT, m_referencePositionXY, m_driftVelocity, m_backwardEndplate);
      newCell.addIndexToTPCDigitIndicesVector(digitindex);
      m_CellsWithHits.insert(make_pair(tpchit.get3DCellID(), newCell));
    } else {
      m_CellsWithHits.at(tpchit.get3DCellID()).addIndexToTPCDigitIndicesVector(digitindex);
    }
    digitindex++;
  }

  for (auto& cellPair : m_CellsWithHits) {
    xCoords.clear();
    yCoords.clear();
    zCoords.clear();

    TPC3DCell cell = cellPair.second;
    const std::vector<uint> tpcDigitIndices = cell.getTPCDigitIndicesInCell();
    const uint nDigits = tpcDigitIndices.size();

    if (nDigits > m_occupancyCut) {
      // B2DEBUG(29, "Too many digits in TPC3DCell (" << nDigits << " which is more than " << m_occupancyCut <<
              // "), mark as bad and continue.");

      cell.setUseCellInTracking(false);
      // cell.getAutomatonCell().setBackgroundFlag(true);
      for (auto& tpcdigitIndex : tpcDigitIndices) {
        hits[tpcdigitIndex].setUseHitInTracking(false);
      }
      m_TPC3DCells.push_back(cell);
      continue;
    }

    for (auto& tpcdigitIndex : tpcDigitIndices) {
      const TVector3& digitRecoPos = hits[tpcdigitIndex].getRecoPosition();
      xCoords.emplace_back(digitRecoPos.X());
      yCoords.emplace_back(digitRecoPos.Y());
      zCoords.emplace_back(digitRecoPos.Z());
    }

    const double xMean = calculateMean(xCoords);
    const double yMean = calculateMean(yCoords);
    const double zMean = calculateMean(zCoords);

    cell.setAverageHitPosition(TVector3(xMean, yMean, zMean));

    m_TPC3DCells.push_back(cell);
  }
}


void simpleTPCBackgroundRejection(std::vector<TPC3DCell>& cells, std::map<TPC3DCellID, TPC3DCell*> CellMap, std::vector<TPC3DCell*> ConnectedCells,
                                  std::vector<hit>& hits)
{
  CellMap.clear();

  for (auto& cell : cells) {
    CellMap.insert(make_pair(cell.getTPC3DCellID(), &cell));
  }

  for (auto& cellPair : CellMap) {
    TPC3DCell* cell = cellPair.second;
    ConnectedCells.clear();
    if (cell->hasTakenFlag() or cell->hasBackgroundFlag()) continue;
    cell->setTakenFlag(true);
    ConnectedCells.emplace_back(cell);
    connectedZCellFinder(cell, CellMap, ConnectedCells);

    // if (m_useNeighbourFinder) {
    //   checkForXYNeighbours();
    // }

    if (ConnectedCells.size() >= m_MaxConnectedCellsInZ) {
      invalidateConnectedCells(ConnectedCells, hits);
    }

    for (auto& unsetCell : ConnectedCells) {
      unsetCell->setTakenFlag(false);
    }

    // find isolated cells with less then m_MinConnectedCellsForIsolationCut connected cells
    // but don't check this cell if it already was invalidated and should not be used to save time
    if (cell->hasTakenFlag() or cell->hasBackgroundFlag()) continue;
    ConnectedCells.clear();
    cell->setTakenFlag(true);
    ConnectedCells.emplace_back(cell);
    isolatedCellFinder(cell, CellMap, ConnectedCells);

    if (ConnectedCells.size() <= m_MinConnectedCellsForIsolationCut) {
      invalidateConnectedCells(ConnectedCells, hits);
    }

    for (auto& unsetCell : ConnectedCells) {
      unsetCell->setTakenFlag(false);
    }

  }
}


void isolatedCellFinder(TPC3DCell* currentCell, std::map<TPC3DCellID, TPC3DCell*>& CellMap, std::vector<TPC3DCell*>& ConnectedCells)
{
  if (ConnectedCells.size() > m_MinConnectedCellsForIsolationCut) return;
  // search for hits
  const vector<TPC3DCellID>& neighbours = currentCell->getTPC3DCellID().getNeighbours();
//   B2INFO("Neighbours size: " << neighbours.size());

  for (auto& neighbour : neighbours) {
    if (ConnectedCells.size() > m_MinConnectedCellsForIsolationCut) return;
    auto neighbourInActiveCells = CellMap.find(neighbour);
    if (neighbourInActiveCells != CellMap.end()) {
      TPC3DCell* nextCell = neighbourInActiveCells->second;

      if (nextCell->hasBackgroundFlag() or nextCell->hasTakenFlag()) continue;

      auto nextCellAlreadyInConnectedCells = std::find(ConnectedCells.begin(), ConnectedCells.end(), nextCell);

      if (nextCellAlreadyInConnectedCells == ConnectedCells.end()) {
        nextCell->setTakenFlag(true);
        ConnectedCells.emplace_back(nextCell);
        if (ConnectedCells.size() > m_MinConnectedCellsForIsolationCut + 1) return;
        isolatedCellFinder(nextCell, CellMap, ConnectedCells);
        if (ConnectedCells.size() > m_MinConnectedCellsForIsolationCut) return;
      }
    }
  }
}


void connectedZCellFinder(TPC3DCell* currentCell, std::map<TPC3DCellID, TPC3DCell*>& CellMap, std::vector<TPC3DCell*>& ConnectedCells)
{
  uint currentZCell = currentCell->getTPC3DCellID().getZCell();
  for (int i = -1; i <= 1; i += 2) {
    TPC3DCellID nextCellID = currentCell->getTPC3DCellID();
    nextCellID.setZCell(currentZCell + i);

    auto cellIDisPresent = CellMap.find(nextCellID);

    if (cellIDisPresent != CellMap.end()) {
      TPC3DCell* cell = cellIDisPresent->second;
      if (cell->hasTakenFlag() or cell->hasBackgroundFlag()) {
        continue;
      } else {
        cell->setTakenFlag(true);
        ConnectedCells.emplace_back(cell);
        connectedZCellFinder(cell, CellMap, ConnectedCells);
      }
    }
  }
}

void invalidateConnectedCells(std::vector<TPC3DCell*>& ConnectedCells, std::vector<hit>& hits)
{
  for (auto& currentCell : ConnectedCells) {
    // a cell is marked as not taken again when it actually has neighours in x,y
    if (not currentCell->hasTakenFlag()) continue;
    currentCell->setUseCellInTracking(false);
    currentCell->setBackgroundFlag(true);
    for (auto& tpcDigitIndex : currentCell->getTPCDigitIndicesInCell()) {
      hits[tpcDigitIndex].setUseHitInTracking(false);
    }
  }
}

double calculateMean(const std::vector<double>& values)
{
  double mean = 0;
  for (auto& value : values) mean += value;
  return mean / values.size();
}