/******************************************************************************

  This source file is part of the Avogadro project.

  Copyright 2013 Kitware, Inc.

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

******************************************************************************/

#include "pdbformat.h"

#include <avogadro/core/elements.h>
#include <avogadro/core/molecule.h>
#include <avogadro/core/utilities.h>
#include <avogadro/core/vector.h>

#include <istream>
#include <string>

using Avogadro::Core::Atom;
using Avogadro::Core::Bond;
using Avogadro::Core::Elements;
using Avogadro::Core::Molecule;
using Avogadro::Core::lexicalCast;
using Avogadro::Core::startsWith;

using std::string;
using std::istringstream;
using std::getline;
using std::vector;

namespace Avogadro {
namespace Io {

PdbFormat::PdbFormat()
{
}

PdbFormat::~PdbFormat()
{
}

bool PdbFormat::read(std::istream& in, Core::Molecule& mol)
{
  string buffer;

  while (getline(in, buffer)) {
    Matrix4f bioMatrix;           // Stores one BIOMT matrix
    vector<Matrix4f> bioMatrices; // Vector of BIOMT matrices

    Matrix4f symMatrix;           // Stores one SMTRY matrix
    vector<Matrix4f> symMatrices; // Vector of SMTRY matrices

    if(startsWith(buffer, "ENDMDL"))
      break;

    else if(startsWith(buffer, "ATOM") || startsWith(buffer, "HETATM")) {
    // 2 spaces after "ATOM" in other
    // codes, why?
      int serial;
      string name;
      char altLoc;    // Alternate location
      string ResName; // Residue name
      char chainId;   // Chain identifier
      int resSeq;     // Residue sequence number
      char iCode;     // Unique ID for inserted residues
      Vector3 pos;    // Coordinates
      float tempFactor;
      string element; // Element symbol, right justififed
      string charge;

      bool ok(false);
      serial = lexicalCast<int>(buffer.substr(6, 5), ok);
      // PDB columns start numbering at 1
      if (!ok) {
        appendError("Error parsing atom serial number");
        return false;
      }

      name = buffer.substr(12, 4);

      altLoc = lexicalCast<char>(buffer.substr(16, 1), ok);
      if (!ok) {
        appendError("Error parsing alternate location");
        return false;
      }

      ResName = buffer.substr(17, 3);
      
      chainId = lexicalCast<char>(buffer.substr(21, 1), ok);
      if (!ok) {
        appendError("Error parsing chain Identification number");
        return false;
      }

      resSeq = lexicalCast<int>(buffer.substr(22, 4), ok);
      if (!ok) {
        appendError("Error parsing residue sequence number");
        return false;
      }

      iCode = lexicalCast<char>(buffer.substr(26, 1), ok);
      if(!ok) {
        appendError("Error parsing ICode");
        return false;
      }

      pos.x() = lexicalCast<Real>(buffer.substr(30, 8), ok);
      if (!ok) {
      appendError("Failed to parse x coordinate: " + buffer.substr(30, 8));
      return false;
      }

      pos.y() = lexicalCast<Real>(buffer.substr(38, 8), ok);
      if (!ok) {
      appendError("Failed to parse y coordinate: " + buffer.substr(38, 8));
      return false;
      }

      pos.z() = lexicalCast<Real>(buffer.substr(46, 8), ok);
      if (!ok) {
      appendError("Failed to parse z coordinate: " + buffer.substr(46, 8));
      return false;
      }

      tempFactor = lexicalCast<float>(buffer.substr(60, 6), ok);
      if (!ok) {
        appendError("Error parsing Temp factor");
        return false;
      }

      element = buffer.substr(76, 2);
      charge = buffer.substr(78, 2);
    }

    else if (startsWith(buffer, "CONECT")) {
      int a;  // Atom serial number
      int b1; // Serial numbers of subsequent bonded atoms
      int b2;
      int b3;
      int b4;

      bool ok(false);
      a = lexicalCast<int>(buffer.substr(6, 5), ok);
      if (!ok) {
        appendError("Error parsing atom serial number");
        return false;
      }

      b1 = lexicalCast<int>(buffer.substr(11, 5), ok);
      if (!ok) {
        appendError("Error parsing serial number of first bonded atom");
        return false;
      }

      b2 = lexicalCast<int>(buffer.substr(16, 5), ok);
      if (!ok) {
        appendError("Error parsing serial number of second bonded atom");
        return false;
      }

      b3 = lexicalCast<int>(buffer.substr(21, 5), ok);
      if (!ok) {
        appendError("Error parsing serial number of third bonded atom");
        return false;
      }

      b4 = lexicalCast<int>(buffer.substr(26, 5), ok);
      if (!ok) {
        appendError("Error parsing serial number of fourth bonded atom");
        return false;
      }
    }

    else if (startsWith(buffer, "HELIX")) { // Space after HELIX?
      int serial;         // Serial number of helix, starts at 1
      string helixId;     // Alphanumeric helix identifier
      string initResName;
      char initChainId;
      int initSeqNum;
      char initICode;
      string endResName;
      char endChainId;
      int endSeqNum;
      char endICode;
      int helixClass;   // Helix type classification
      int length;       // Length of the helix

      bool ok(false);
      serial = lexicalCast<int>(buffer.substr(7, 3), ok);
      if (!ok) {
        appendError("Error parsing serial number");
        return false;
      }

      helixId = buffer.substr(11, 3);
      initResName = buffer.substr(15, 3);

      initChainId = lexicalCast<char>(buffer.substr(19, 1), ok);
      if (!ok) {
        appendError("Error parsing initial chain ID");
        return false;
      }

      initSeqNum = lexicalCast<int>(buffer.substr(21, 4), ok);
      if (!ok) {
        appendError("Error parsing initial sequence number");
        return false;
      }

      initICode = lexicalCast<char>(buffer.substr(25, 1), ok);
      if (!ok) {
        appendError("Error parsing initial ICode");
        return false;
      }

      endResName = buffer.substr(27, 3);

      endChainId = lexicalCast<char>(buffer.substr(31, 1), ok);
      if (!ok) {
        appendError("Error parsing end chain ID");
        return false;
      }

      endSeqNum = lexicalCast<int>(buffer.substr(33, 4), ok);
      if (!ok) {
        appendError("Error parsing end sequence number");
        return false;
      }

      endICode = lexicalCast<char>(buffer.substr(37, 1), ok);
      if (!ok) {
        appendError("Error parsing end ICode");
        return false;
      }

      helixClass = lexicalCast<int>(buffer.substr(38, 2), ok);
      if (!ok) {
        appendError("Error parsing helix class");
        return false;
      }

      length = lexicalCast<int>(buffer.substr(71, 5), ok);
      if (!ok) {
        appendError("Error parsing helix length");
        return false;
      }
    }

    else if (startsWith(buffer, "SHEET")) { // Space after SHEET?
      int strand;
      string sheetId;
      int numStrands;
      string initResName;
      char initChainId;
      int initSeqNum;
      char initICode;
      string endResName;
      char endChainId;
      int endSeqNum;
      char endICode;
      int sense;        // Sense of strand wrt previous strand, -1, 0 or 1
      string curAtom;   // Cur refers to current strand
      string curResName;
      char curChainId;
      int curResSeq;
      char curICode;
      string prevAtom;
      string prevResName;
      char prevChainId;
      int prevResSeq;
      char prevICode;

      bool ok(false);
      strand = lexicalCast<int>(buffer.substr(7, 3), ok);
      if (!ok) {
        appendError("Error parsing strand number");
        return false;
      }

      sheetId = buffer.substr(11, 3);
      
      numStrands = lexicalCast<int>(buffer.substr(14, 2), ok);
      if (!ok) {
        appendError("Error parsing number of strands");
        return false;
      }

      initResName = buffer.substr(17, 3);
      
      initChainId = lexicalCast<char>(buffer.substr(21, 1), ok);
      if (!ok) {
        appendError("Error parsing end initial chain ID");
        return false;
      }

      initSeqNum = lexicalCast<int>(buffer.substr(22, 4), ok);
      if (!ok) {
        appendError("Error parsing initial sequence number");
        return false;
      }

      initICode = lexicalCast<char>(buffer.substr(26, 1), ok);
      if (!ok) {
        appendError("Error parsing initial ICode");
        return false;
      }

      endResName = buffer.substr(28, 3);

      endChainId = lexicalCast<char>(buffer.substr(32, 1), ok);
      if (!ok) {
        appendError("Error parsing end chain ID");
        return false;
      }

      endSeqNum = lexicalCast<int>(buffer.substr(33, 4), ok);
      if (!ok) {
        appendError("Error parsing end sequence number");
        return false;
      }

      endICode = lexicalCast<char>(buffer.substr(37, 1), ok);
      if (!ok) {
        appendError("Error parsing end ICode");
        return false;
      }

      sense = lexicalCast<int>(buffer.substr(38, 2), ok);
      if (!ok) {
        appendError("Error parsing sense");
        return false;
      }

      curAtom = buffer.substr(41, 4);
      curResName = buffer.substr(45, 3);

      curChainId = lexicalCast<char>(buffer.substr(49, 1), ok);
      if (!ok) {
        appendError("Error parsing chain ID");
        return false;
      }

      curResSeq = lexicalCast<int>(buffer.substr(50, 4), ok);
      if (!ok) {
        appendError("Error parsing current residue sequence");
        return false;
      }

      curICode = lexicalCast<char>(buffer.substr(54, 1), ok);
      if (!ok) {
        appendError("Error parsing current ICode");
        return false;
      }

      prevAtom = buffer.substr(56, 4);
      prevResName = buffer.substr(60, 3);

      prevChainId = lexicalCast<char>(buffer.substr(64, 1), ok);
      if (!ok) {
        appendError("Error parsing previous chain ID");
        return false;
      }

      prevResSeq = lexicalCast<int>(buffer.substr(65, 4), ok);
      if (!ok) {
        appendError("Error parsing previous residue sequence");
        return false;
      }

      prevICode = lexicalCast<char>(buffer.substr(69, 1), ok);
      if (!ok) {
        appendError("Error parsing previous ICode");
        return false;
      }
    }

    else if(startsWith(buffer, "REMARK 350   BIOMT")) {
      // Further possible optimization by reading all matrices together
      bool ok(false);
      int row = lexicalCast<int>(buffer.substr(18, 1), ok) - 1;
      if (!ok) {
        appendError("Error parsing BIOMT matrix row number");
        return false;
      }

      bioMatrix(row, 0) = lexicalCast<float>(buffer.substr(23, 10), ok);
      if (!ok) {
        appendError("Error parsing BIOMT matrix column 0");
        return false;
      }


      bioMatrix(row, 1) = lexicalCast<float>(buffer.substr(33, 10), ok);
      if (!ok) {
        appendError("Error parsing BIOMT matrix column 1");
        return false;
      }


      bioMatrix(row, 2) = lexicalCast<float>(buffer.substr(43, 10), ok);
      if (!ok) {
        appendError("Error parsing BIOMT matrix column 2");
        return false;
      }


      bioMatrix(row, 3) = lexicalCast<float>(buffer.substr(53, 5), ok);
      if (!ok) {
        appendError("Error parsing BIOMT matrix column 3");
        return false;
      }

      if (row == 2)
        bioMatrices.push_back(bioMatrix);
    }

    else if(startsWith(buffer, "REMARK 290   SMTRY")) {
      // Further possible optimization by reading all matrices together
      bool ok(false);
      int row = lexicalCast<int>(buffer.substr(18, 1), ok) - 1;
      if (!ok) {
        appendError("Error parsing SMTRY matrix row number");
        return false;
      }

      symMatrix(row, 0) = lexicalCast<float>(buffer.substr(23, 10), ok);
      if (!ok) {
        appendError("Error parsing SMTRY matrix column 0");
        return false;
      }


      symMatrix(row, 1) = lexicalCast<float>(buffer.substr(33, 10), ok);
      if (!ok) {
        appendError("Error parsing SMTRY matrix column 1");
        return false;
      }


      symMatrix(row, 2) = lexicalCast<float>(buffer.substr(43, 10), ok);
      if (!ok) {
        appendError("Error parsing SMTRY matrix column 2");
        return false;
      }


      symMatrix(row, 3) = lexicalCast<float>(buffer.substr(53, 5), ok);
      if (!ok) {
        appendError("Error parsing SMTRY matrix column 3");
        return false;
      }

      if (row == 2)
        symMatrices.push_back(symMatrix);
    }
  }
  return true;
}

std::vector<std::string> PdbFormat::fileExtensions() const
{
  std::vector<std::string> ext;
  ext.push_back("pdb");
  return ext;
}

std::vector<std::string> PdbFormat::mimeTypes() const
{
  std::vector<std::string> mime;
  mime.push_back("chemical/x-pdb");
  return mime;
}

} // end Io namespace
} // end Avogadro namespace