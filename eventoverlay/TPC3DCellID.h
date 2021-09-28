/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include <vector>
#include <string>
#include <ostream>


  /**
   * Class to uniquely identify a "3D sensor" in the TPC volume.
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
  class TPC3DCellID {
  public:
    /** The base integer type for TPC3DCells */
    typedef unsigned int baseType;
    enum {
      /** Number of bits available to represent the TCPTPC3DCell in x */
      xBits = 8,
      /** Number of bits available to represent the TCPTPC3DCell in y */
      yBits = 8,
      /** Number of bits available to represent the TCPTPC3DCell in z */
      zBits = 8,
      /** Total bit size of the TPC3DCells */
      Bits  = xBits + yBits + zBits,

      /** Maximum valid x CellID */
      MaxX  = (1 << xBits) - 1,
      /** Maximum valid y CellID */
      MaxY  = (1 << yBits) - 1,
      /** Maximum valid z CellID */
      MaxZ  = (1 << zBits) - 1,
      /** Maximum value for ID */
      MaxID = (1 << Bits) - 1
    };

    /** Constructor using the unique id */
    // cppcheck-suppress noExplicitConstructor
    TPC3DCellID(baseType id = 0)
    {
      m_id.id = id;
    }
    /** Constructor using layer, ladder and sensor ids */
    TPC3DCellID(baseType xCell, baseType yCell, baseType zCell)
    {
      m_id.parts.xCell  = xCell;
      m_id.parts.yCell  = yCell;
      m_id.parts.zCell  = zCell;
    }
    /** Construct ID from string representing the structure */
    explicit TPC3DCellID(const std::string& sensor);
    /** Copy constructor */
    TPC3DCellID(const TPC3DCellID& b): m_id(b.m_id) {}

    /** Assignment operator */
    TPC3DCellID& operator=(const TPC3DCellID& b)        { m_id = b.m_id; return *this; }
    /** Assignment from baseType */
    TPC3DCellID& operator=(baseType id)     { m_id.id = id; return *this; }
    /** Convert to baseType */
    operator baseType() const         { return getID(); }
    /** Convert to string */
    operator std::string() const;
    /** Check for equality */
    bool operator==(const TPC3DCellID& b) const   { return getID() == b.getID(); }
    /** Order by unique id */
    bool operator<(const TPC3DCellID& b) const    { return getID() < b.getID(); }

    /** Get the unique id */
    baseType getID() const            { return m_id.id; }
    /** Get the x cell id */
    baseType getXCell() const         { return m_id.parts.xCell; }
    /** Get the y cell id */
    baseType getYCell() const         { return m_id.parts.yCell; }
    /** Get the z cell id */
    baseType getZCell() const         { return m_id.parts.zCell; }

    /** Set the unique id */
    void setID(baseType id)           { m_id.id = id; }
    /** Set the x cell id */
    void setXCell(baseType xCell)     { m_id.parts.xCell  = xCell;  }
    /** Set the y cell id */
    void setYCell(baseType yCell)     { m_id.parts.yCell  = yCell;  }
    /** Set the z cell id */
    void setZCell(baseType zCell)     { m_id.parts.zCell  = zCell;  }

    /** Get (up to) 26 adjacent cells in all 3 spatial dimensions */
    std::vector<TPC3DCellID> getNeighbours() const
    {

      std::vector<TPC3DCellID> neighbours;
      neighbours.reserve(26); // There are at maximum 26 cells around the current cell

      const int currentX = m_id.parts.xCell;
      const int currentY = m_id.parts.yCell;
      const int currentZ = m_id.parts.zCell;

      for (int x = -1; x <= +1; x++) {
        // don't check out of bounds cells
        if ((x == -1 and currentX == 0) or (x == +1 and currentX == 127)) continue;
        for (int y = -1; y <= +1; y++) {
          // don't check out of bounds cells
          if ((y == -1 and currentY == 0) or (y == +1 and currentY == 127)) continue;
          for (int z = -1; z <= +1; z++) {
            // don't check out of bounds cells
            if ((z == -1 and currentZ == 0) or (z == +1 and currentZ == 127)) continue;
            // cell can't be a neigbour of itself
            if (x == 0 and y == 0 and z == 0) continue;
            neighbours.emplace_back(TPC3DCellID(currentX + x, currentY + y, currentZ + z));
          }
        }
      }
      return neighbours;
    }

    /** Get (up to) 26 adjacent cells in all 3 spatial dimensions */
    std::vector<TPC3DCellID> getNeighbours(const TPC3DCellID cell)
    {

      std::vector<TPC3DCellID> neighbours;
      neighbours.reserve(26); // There are at maximum 26 cells around the current cell

      const int currentX = cell.getXCell();
      const int currentY = cell.getYCell();
      const int currentZ = cell.getZCell();

      for (int x = -1; x <= +1; x++) {
        // don't check out of bounds cells
        if ((x == -1 and currentX == 0) or (x == +1 and currentX == 127)) continue;
        for (int y = -1; y <= +1; y++) {
          // don't check out of bounds cells
          if ((y == -1 and currentY == 0) or (y == +1 and currentY == 127)) continue;
          for (int z = -1; z <= +1; z++) {
            // don't check out of bounds cellsi
            if ((z == -1 and currentZ == 0) or (z == +1 and currentZ == 127)) continue;
            // cell can't be a neigbour of itself
            if (x == 0 and y == 0 and z == 0) continue;
            neighbours.emplace_back(TPC3DCellID(currentX + x, currentY + y, currentZ + z));
          }
        }
      }
      return neighbours;
    }

    /** make this type printable in python with print(TPC3DCellid) */
    std::string __str__() const { return (std::string)(*this); }

  private:

    union {
      /** Unique id */
baseType id: Bits;
      struct {
        /** xCell id */
baseType xCell: xBits;
        /** yCell id */
baseType yCell: yBits;
        /** zCell id */
baseType zCell: zBits;
      } parts /**< Struct to contain all id components */;
    } m_id; /**< Union to store the ID and all components in one go. */
  };

  /** Print id to stream by converting it to string */
  std::ostream& operator<<(std::ostream& out, const TPC3DCellID& id);

