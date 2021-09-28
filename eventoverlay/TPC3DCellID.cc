/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/

#include <TPC3DCellID.h>
#include <sstream>

using namespace std;


  //namespace {
    /**
     * Small helper function to parse TPC3DCellID string representation
     *
     * This function takes an input stream and will return the next component of the TPC3DCellID
     * */
    int getPart(istream& in)
    {
      if (!in.eof()) {
        //Get next char, if it is a dot, ignore it and get the next one
        int next = in.get();
        if (next == '.') next = in.get();
        //If it is a wildcard we return 0 as id, otherwise we put it back in the stream
        if (next == '*' or in.eof()) {
          return 0;
        } else {
          in.unget();
        }
        //If it is the segment separator, we assume the remaining parts to be missing, so return 0
        if (next == '#') return 0;

        //Now get the actual value out of the stream. If this fails something is wrong and it is not
        //a valid id
        int value(0);
        in >> value;
        if (in.fail() && !in.eof()) {
          throw runtime_error("Failed to parse Number");
        }
        return value;
      }
      return 0;
    }
  }

  TPC3DCellID::TPC3DCellID(const std::string& sensor)
  {
    //We parse the Id from string, so set it to 0 first
    m_id.id = 0;
    //create a stream from the string
    istringstream in(sensor);
    try {
      //Get all the parts
      m_id.parts.xCell  = getPart(in);
      m_id.parts.yCell  = getPart(in);
      m_id.parts.zCell  = getPart(in);

    } catch (runtime_error&) {
      //Something went wrong parsing the parts
      m_id.id = 0;
      throw invalid_argument("Could not parse TPC3DCellID: '" + sensor + "'");
    }
    //There is stuff left we also throw an exception as we cannot warn the user
    //without the logging system
    if (!in.eof()) {
      string rest;
      //Get the remainder: get everything in the stream until the next NULL
      //character which should only occur at the end of the string.
      getline(in, rest, '\0');
      throw invalid_argument("Trailing characters after TPC3DCellID " + (string)*this + ": '" + rest + "'");
    }
  }

  TPC3DCellID::operator string() const
  {
    stringstream out;
    if (m_id.parts.xCell) {
      out << m_id.parts.xCell;
    } else {
      out << "*";
    }
    if (m_id.parts.yCell) {
      out << "." << m_id.parts.yCell;
    }
    if (m_id.parts.zCell) {
      out << "." << m_id.parts.zCell;
    }
    return out.str();
  }

  std::ostream& operator<<(std::ostream& out, const TPC3DCellID& id)
  {
    out << ((string)id);
    return out;
  }

  std::vector<TPC3DCellID> TPC3DCellID::getNeighbours() const
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

  std::vector<TPC3DCellID> TPC3DCellID::getNeighbours(const TPC3DCellID cell)
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

  //}


