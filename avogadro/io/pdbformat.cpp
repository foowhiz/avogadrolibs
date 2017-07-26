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
  Matrix4f matrix = Matrix4f::Identity(); // Stores one BIOMT or SMTRY matrix
  vector<Matrix4f> bioMatrices;           // Vector of BIOMT matrices
  vector<Matrix4f> symMatrices;           // Vector of SMTRY matrices

  string buffer;

  while (getline(in, buffer)) { // Read Each line one by one

    if (startsWith(buffer, "ENDMDL"))
      break;

    else if (startsWith(buffer, "ATOM") || startsWith(buffer, "HETATM")) {
      // 2 spaces after "ATOM" in other codes, why?
      int serial;
      bool ok(false);
      serial = lexicalCast<int>(buffer.substr(6, 5), ok);
      // PDB columns start numbering at 1
      if (!ok) {
        appendError("Error parsing atom serial number");
        return false;
      }

      string name;
      name = buffer.substr(12, 4);

      char altLoc; // Alternate location
      if (buffer.substr(16, 1) != " ") {
        altLoc = lexicalCast<char>(buffer.substr(16, 1), ok);
        if (!ok) {
          appendError("Error parsing alternate location");
          return false;
        }
      }

      string ResName; // Residue name
      ResName = buffer.substr(17, 3);

      char chainId; // Chain identifier
      chainId = lexicalCast<char>(buffer.substr(21, 1), ok);
      if (!ok) {
        appendError("Error parsing chain identification number");
        return false;
      }

      int resSeq; // Residue sequence number
      resSeq = lexicalCast<int>(buffer.substr(22, 4), ok);
      if (!ok) {
        appendError("Error parsing residue sequence number");
        return false;
      }

      if (buffer.substr(26, 1) != " ") // Record may not be present
      // TODO: Change to find_first_not_of
      {
        char iCode; // Unique ID for inserted residues
        iCode = lexicalCast<char>(buffer.substr(26, 1), ok);
        if (!ok) {
          appendError("Error parsing ICode");
          return false;
        }
      }

      Vector3 pos; // Coordinates
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

      float tempFactor;
      tempFactor = lexicalCast<float>(buffer.substr(60, 6), ok);
      if (!ok) {
        appendError("Error parsing Temp factor");
        return false;
      }

      string element; // Element symbol, right justififed
      element = buffer.substr(76, 2);

      string charge;
      charge = buffer.substr(78, 2);
    }

    else if (startsWith(buffer, "CONECT")) {

      int a; // Atom serial number
      bool ok(false);
      a = lexicalCast<int>(buffer.substr(6, 5), ok);
      if (!ok) {
        appendError("Error parsing atom serial number");
        return false;
      }

      int b1; // Serial numbers of subsequent bonded atoms
      b1 = lexicalCast<int>(buffer.substr(11, 5), ok);
      if (!ok) {
        appendError("Error parsing serial number of first bonded atom");
        return false;
      }

      int b2;
      b2 = lexicalCast<int>(buffer.substr(16, 5), ok);
      if (!ok &&
          (buffer.substr(16, 5).find_first_not_of(' ') != std::string::npos)) {
        appendError(buffer.substr(16, 5));
        return false;
      }

      int b3;
      b3 = lexicalCast<int>(buffer.substr(21, 5), ok);
      if (!ok &&
          (buffer.substr(21, 5).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing serial number of third bonded atom");
        return false;
      }

      int b4;
      b4 = lexicalCast<int>(buffer.substr(26, 5), ok);
      if (!ok &&
          (buffer.substr(26, 5).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing serial number of fourth bonded atom");
        return false;
      }
      // Up to how many bonds?
    }

    else if (startsWith(buffer, "HELIX")) { // Space after HELIX?

      int serial; // Serial number of helix, starts at 1
      bool ok(false);
      serial = lexicalCast<int>(buffer.substr(7, 3), ok);
      if (!ok) {
        appendError("Error parsing serial number");
        return false;
      }

      string helixId; // Alphanumeric helix identifier
      helixId = buffer.substr(11, 3);

      string initResName;
      initResName = buffer.substr(15, 3);

      char initChainId;
      initChainId = lexicalCast<char>(buffer.substr(19, 1), ok);
      if (!ok) {
        appendError("Error parsing initial chain ID");
        return false;
      }

      int initSeqNum;
      initSeqNum = lexicalCast<int>(buffer.substr(21, 4), ok);
      if (!ok) {
        appendError("Error parsing initial sequence number");
        return false;
      }

      char initICode;
      initICode = lexicalCast<char>(buffer.substr(25, 1), ok);
      if (!ok &&
          (buffer.substr(25, 1).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing initial ICode");
        return false;
      }

      string endResName;
      endResName = buffer.substr(27, 3);

      char endChainId;
      endChainId = lexicalCast<char>(buffer.substr(31, 1), ok);
      if (!ok) {
        appendError("Error parsing end chain ID");
        return false;
      }

      int endSeqNum;
      endSeqNum = lexicalCast<int>(buffer.substr(33, 4), ok);
      if (!ok) {
        appendError("Error parsing end sequence number");
        return false;
      }

      char endICode;
      endICode = lexicalCast<char>(buffer.substr(37, 1), ok);
      if (!ok &&
          (buffer.substr(37, 1).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing end ICode");
        return false;
      }

      int helixClass; // Helix type classification
      helixClass = lexicalCast<int>(buffer.substr(38, 2), ok);
      if (!ok) {
        appendError("Error parsing helix class");
        return false;
      }

      int length; // Length of the helix
      length = lexicalCast<int>(buffer.substr(71, 5), ok);
      if (!ok) {
        appendError("Error parsing helix length");
        return false;
      }
    }

    else if (startsWith(buffer, "SHEET")) { // Space after SHEET?

      int strand;
      bool ok(false);
      strand = lexicalCast<int>(buffer.substr(7, 3), ok);
      if (!ok) {
        appendError("Error parsing strand number");
        return false;
      }

      string sheetId;
      sheetId = buffer.substr(11, 3);

      int numStrands;
      numStrands = lexicalCast<int>(buffer.substr(14, 2), ok);
      if (!ok) {
        appendError("Error parsing number of strands");
        return false;
      }

      string initResName;
      initResName = buffer.substr(17, 3);

      char initChainId;
      initChainId = lexicalCast<char>(buffer.substr(21, 1), ok);
      if (!ok) {
        appendError("Error parsing end initial chain ID");
        return false;
      }

      int initSeqNum;
      initSeqNum = lexicalCast<int>(buffer.substr(22, 4), ok);
      if (!ok) {
        appendError("Error parsing initial sequence number");
        return false;
      }

      char initICode;
      initICode = lexicalCast<char>(buffer.substr(26, 1), ok);
      if (!ok &&
          (buffer.substr(26, 1).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing initial ICode");
        return false;
      }

      string endResName;
      endResName = buffer.substr(28, 3);

      char endChainId;
      endChainId = lexicalCast<char>(buffer.substr(32, 1), ok);
      if (!ok) {
        appendError("Error parsing end chain ID");
        return false;
      }

      int endSeqNum;
      endSeqNum = lexicalCast<int>(buffer.substr(33, 4), ok);
      if (!ok) {
        appendError("Error parsing end sequence number");
        return false;
      }

      char endICode;
      endICode = lexicalCast<char>(buffer.substr(37, 1), ok);
      if (!ok &&
          (buffer.substr(37, 1).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing end ICode");
        return false;
      }

      int sense; // Sense of strand wrt previous strand, -1, 0 or 1
      sense = lexicalCast<int>(buffer.substr(38, 2), ok);
      if (!ok) {
        appendError("Error parsing sense");
        return false;
      }

      string curAtom;
      curAtom = buffer.substr(41, 4);

      string curResName;
      curResName = buffer.substr(45, 3);

      char curChainId;
      curChainId = lexicalCast<char>(buffer.substr(49, 1), ok);
      if (!ok &&
          (buffer.substr(49, 1).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing chain ID");
        return false;
      }

      int curResSeq;
      curResSeq = lexicalCast<int>(buffer.substr(50, 4), ok);
      if (!ok &&
          (buffer.substr(50, 4).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing current residue sequence");
        return false;
      }

      char curICode;
      curICode = lexicalCast<char>(buffer.substr(54, 1), ok);
      if (!ok &&
          (buffer.substr(54, 1).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing current ICode");
        return false;
      }

      string prevAtom;
      prevAtom = buffer.substr(56, 4);

      string prevResName;
      prevResName = buffer.substr(60, 3);

      char prevChainId;
      prevChainId = lexicalCast<char>(buffer.substr(64, 1), ok);
      if (!ok &&
          (buffer.substr(64, 1).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing previous chain ID");
        return false;
      }

      int prevResSeq;
      prevResSeq = lexicalCast<int>(buffer.substr(65, 4), ok);
      if (!ok &&
          (buffer.substr(65, 4).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing previous residue sequence");
        return false;
      }

      char prevICode;
      prevICode = lexicalCast<char>(buffer.substr(69, 1), ok);
      if (!ok &&
          (buffer.substr(69, 1).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing previous ICode");
        return false;
      }
    }

    else if (startsWith(buffer, "REMARK 350   BIOMT")) {
      // Further possible optimization by reading all matrices together
      bool ok(false);
      int row;
      row = lexicalCast<int>(buffer.substr(18, 1), ok) - 1;
      if (!ok) {
        appendError("Error parsing BIOMT matrix row number");
        return false;
      }

      matrix(row, 0) = lexicalCast<float>(buffer.substr(23, 10), ok);
      if (!ok) {
        appendError("Error parsing BIOMT matrix column 0");
        return false;
      }

      matrix(row, 1) = lexicalCast<float>(buffer.substr(33, 10), ok);
      if (!ok) {
        appendError("Error parsing BIOMT matrix column 1");
        return false;
      }

      matrix(row, 2) = lexicalCast<float>(buffer.substr(43, 10), ok);
      if (!ok) {
        appendError("Error parsing BIOMT matrix column 2");
        return false;
      }

      matrix(row, 3) = lexicalCast<float>(buffer.substr(55, 13), ok);
      if (!ok) {
        appendError("Error parsing BIOMT matrix column 3");
        return false;
      }

      if (row == 2)
        bioMatrices.push_back(matrix);
    }

    else if (startsWith(buffer, "REMARK 290   SMTRY")) {
      // Further possible optimization by reading all matrices together
      bool ok(false);
      int row;
      row = lexicalCast<int>(buffer.substr(18, 1), ok) - 1;
      if (!ok) {
        appendError("Error parsing SMTRY matrix row number");
        return false;
      }

      matrix(row, 0) = lexicalCast<float>(buffer.substr(23, 10), ok);
      if (!ok) {
        appendError("Error parsing SMTRY matrix column 0");
        return false;
      }

      matrix(row, 1) = lexicalCast<float>(buffer.substr(33, 10), ok);
      if (!ok) {
        appendError("Error parsing SMTRY matrix column 1");
        return false;
      }

      matrix(row, 2) = lexicalCast<float>(buffer.substr(43, 10), ok);
      if (!ok) {
        appendError("Error parsing SMTRY matrix column 2");
        return false;
      }

      matrix(row, 3) = lexicalCast<float>(buffer.substr(55, 13), ok);
      if (!ok) {
        appendError("Error parsing SMTRY matrix column 3");
        return false;
      }

      if (row == 2)
        symMatrices.push_back(matrix);
    }
  } // End while loop
  return true;
} // End read

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