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

      string name;
      name = buffer.substr(12, 4);
      name = trimmed(name);

      if (name == "CA" || name == "O")  // Assuming CA and O occur alternatingly
      {
        Vector3 pos; // Coordinates
        bool ok(false);
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

        string element; // Element symbol, right justififed
        element = buffer.substr(76, 2);
        element = trimmed(element);

        unsigned char atomicNum = Elements::atomicNumberFromSymbol(element);
        //appendError(std::to_string(atomicNum));
        Atom newAtom = mol.addAtom(atomicNum);
        newAtom.setPosition3d(pos);
        newAtom.setAtomName(name); // May require casting of serial to Index type
      }
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