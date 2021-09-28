/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/

#pragma once

//#include <framework/datastore/RelationsObject.h>
//#include <framework/geometry/B2Vector3.h>
//#include <framework/gearbox/Unit.h>

#include "TPC3DCellID.h"

// #include <tracking/trackFindingCDC/ca/AutomatonCell.h>

#include <vector>
#include <string>
#include <ostream>

  /**
   * Class for a 3D cell hit in the TPC volume.
   *
   * - xCell, yCell and zCell all start 0
   * - Reference (0, 0) for x, y is in the (-x, -y) coordinate quadrant (lower left corner in x-y-plane)
   * - zCells start with 0 at the readout plane
   *
   * Internal use of a union gets rid of all the bit shifting which would be
   * neccessary to represent the id as one baseType and get all the
   * components out of it. Disadvantage is that it is not guaranteed to be
   * portable, but neither is bit shifting
   */
  class TPC3DCell {
  public:
    /** Constructor using the unique id */
    TPC3DCell()
    { }

    /** Constructor using layer, ladder and sensor ids */
    TPC3DCell(const TPC3DCellID cellid = 0,
           const double cellSizeXY = 2.,
           const double cellSizeZT = 250.,
           const double referencePositionXY = -128.,
           const double driftVelocity = 77./10000,
           const double readoutPlaneZPos = -83.12) :
      m_cellID(cellid)
    {
      double x = referencePositionXY + cellSizeXY * ((double)cellid.getXCell() + 0.5);
      double y = referencePositionXY + cellSizeXY * ((double)cellid.getYCell() + 0.5);
      double z = readoutPlaneZPos + cellSizeZT * driftVelocity * ((double)cellid.getZCell() + 0.5);
      m_cellCenterPosition.SetXYZ(x, y, z);
    }

    /** Constructor using layer, ladder and sensor ids */
    TPC3DCell(const TPC3DCellID cellid, const TVector3 cellcenter) : m_cellID(cellid), m_cellCenterPosition(cellcenter)
    { }

    /** Get the unique id */
    TPC3DCellID getTPC3DCellID() const            { return m_cellID; }

    /** Get cell position */
    TVector3 getCellCenterPosition() const { return m_cellCenterPosition; }

    /** Get average hit position */
    TVector3 getAverageHitPosition() const { return m_averageHitPosition; }

    /** Get SVD (Singular Value Decomposition) major axis */
    TVector3 getSVDMajorAxis() const { return m_SVDMajorAxis; }

    /** Get SVD (Singular Value Decomposition) semi major axis */
    TVector3 getSVDSemiMajorAxis() const { return m_SVDSemiMajorAxis; }

    /** Get SVD (Singular Value Decomposition) least axis */
    TVector3 getSVDLeastAxis() const { return m_SVDLeastAxis; }

    /** Set cell position */
    void setCellCenterPosition(const TVector3 cellPosition) { m_cellCenterPosition = cellPosition; }

    /** Set average hit position */
    void setAverageHitPosition(const TVector3 hitPosition) { m_averageHitPosition = hitPosition; }

    /** Set SVD major axis */
    void setSVDMajorAxis(const TVector3 major) { m_SVDMajorAxis = major; }

    /** Set SVD semi major axis */
    void setSVDSemiMajorAxis(const TVector3 semimajor) { m_SVDSemiMajorAxis = semimajor; }

    /** Set SVD least axis */
    void setSVDLeastAxis(const TVector3 least) { m_SVDLeastAxis = least; }

    /** Get list of TPCDigit StoreArray indices contained in this cell (instead of setting all relations) */
    std::vector<uint> getTPCDigitIndicesInCell() const { return m_TPCDigitIndices; }

    void addIndexToTPCDigitIndicesVector(const uint hitindex) { m_TPCDigitIndices.push_back(hitindex); }

    void setTPCDigitIndicesInCell(const std::vector<uint>& hitindices) { m_TPCDigitIndices = hitindices; }

    /** Check for equality */
    bool operator==(const TPC3DCell& b) const   { return getTPC3DCellID() == b.getTPC3DCellID(); }

    /** Check for equality */ //TODO I had to out-comment this (02.05.2021), but I don't remember why I included it in the first place.
//     bool operator==(const TPCDigit& b) const   { return getTPC3DCellID() == b.getTPC3DCellID(); }

    /// Getter for the flag indicating whether this digit should be used in tracking
    bool useCellInTracking() const { return m_useCellInTracking; }

    /// set the flag to control whether this digit should be used in tracking
    void setUseCellInTracking(const bool useInTracking) { m_useCellInTracking = useInTracking; }

    /// Getter for the automaton cell.
    // TrackFindingCDC::AutomatonCell& getAutomatonCell() { return m_automatonCell; }

    bool hasTakenFlag() { return c_takenFlag; }
    bool hasBackgroundFlag() { return c_bkgFlag; }

    void setTakenFlag(const bool flag) { c_takenFlag = flag; }
    void setBackgroundFlag(const bool flag) { c_bkgFlag = flag; }


  private:
    /// TPC3DCellID of this cell
    TPC3DCellID m_cellID;

    /// center position of this cell
    TVector3 m_cellCenterPosition;

    /// average position of hits in this cell, not weighted by errors
    TVector3 m_averageHitPosition;

    /// Singular Value Decomposition major axis of hits contained in the cell
    TVector3 m_SVDMajorAxis;

    /// Singular Value Decomposition semi major axis of hits contained in the cell
    TVector3 m_SVDSemiMajorAxis;

    /// Singular Value Decomposition major axis of hits contained in the cell
    TVector3 m_SVDLeastAxis;


    /// vector containing the indices of the TPCDigits that are in this 3D cell
    std::vector<uint> m_TPCDigitIndices;

    /// flag to ignore a cell and the digits contained in it in further steps
    bool m_useCellInTracking = true;

    bool c_takenFlag = false;
    bool c_bkgFlag = false;

    /// CellularAutomatonCell with further status info
    // mutable TrackFindingCDC::AutomatonCell m_automatonCell;  //! ROOT streamer


    //ClassDef(TPC3DCell, 2);    /**< ClassDef */

  };

