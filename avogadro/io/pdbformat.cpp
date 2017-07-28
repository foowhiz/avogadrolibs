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
using Avogadro::Core::trimmed;

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
        appendError(std::to_string(serial));
        return false;
      }

      string name;
      name = buffer.substr(12, 4);
      name = trimmed(name);

      /*char altLoc; // Alternate location
      if (buffer.substr(16, 1) != " ") {
        altLoc = lexicalCast<char>(buffer.substr(16, 1), ok);
        if (!ok) {
          appendError("Error parsing alternate location");
          return false;
        }
      }

      appendError(std::to_string(altLoc));

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

      char iCode; // Unique ID for inserted residues
      iCode = lexicalCast<char>(buffer.substr(26, 1), ok);
      if (!ok &&
        (buffer.substr(26, 1).find_first_not_of(' ') != std::string::npos)) {
        appendError("Error parsing ICode");
        return false;
      }*/

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

      /*float tempFactor;
      tempFactor = lexicalCast<float>(buffer.substr(60, 6), ok);
      if (!ok) {
        appendError("Error parsing Temp factor");
        return false;
      }*/

      string element; // Element symbol, right justififed
      element = buffer.substr(76, 2);
      element = trimmed(element);

      /*string charge;
      charge = buffer.substr(78, 2);*/

      unsigned char atomicNum = Elements::atomicNumberFromSymbol(element);
      //appendError(std::to_string(atomicNum));
      Atom newAtom = mol.addAtom(atomicNum);
      newAtom.setPosition3d(pos);
      newAtom.setAtomName(serial, name); // May require casting of serial to Index type
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