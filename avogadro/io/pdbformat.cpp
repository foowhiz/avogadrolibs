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
    if(startsWith(buffer, "ENDMDL"))
      break;

    else if(startsWith(buffer, "ATOM")) {   // 2 spaces after "ATOM" in other
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
      serial = (lexicalCast<int>(buffer.substr(6, 5), ok));
      // PDB columns start numbering at 1
      if (!ok) {
        appendError("Error parsing atom serial number");
        return false;
      }

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