#ifndef AVOGADRO_IO_PDBFORMAT_H
#define AVOGADRO_IO_PDBFORMAT_H

#include "fileformat.h"

namespace Avogadro {
namespace Io {

/**
 * @class PdbFormat pdbformat.h <avogadro/io/pdbformat.h>
 * @brief Parser for the PDB format.
 * @author Tanuj Kumar
 */

class AVOGADROIO_EXPORT PdbFormat : public FileFormat
{
public:
  PdbFormat();
  ~PdbFormat() override;

  Operations supportedOperations() const override
  {
    return Read; //Unsure of what all should be there
  }

  FileFormat* newInstance() const override { return new PdbFormat; }
  std::string identifier() const override { return "Avogadro: PDB"; }
  std::string name() const override { return "PDB"; }
  std::string description() const override
  {
    return "Generic format that contains atoms, bonds, positions."; //To be updated, copied from mdlformat.h
  }

  std::string specificationUrl() const override
  {
    return "http://www.wwpdb.org/documentation/file-format-content/"
           "format33/v3.3.html";
  }

} // end Io namespace
} // end Avogadro namespace
#endif // AVOGADRO_IO_PDBFORMAT_H