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

#include "iotests.h"

#include <gtest/gtest.h>

#include <avogadro/core/molecule.h>

#include <avogadro/io/pdbformat.h>

using Avogadro::Core::Molecule;
using Avogadro::Core::Atom;
using Avogadro::Core::Bond;
using Avogadro::Core::Variant;
using Avogadro::Io::FileFormat;
using Avogadro::Io::PdbFormat;

TEST(PdbTest, readFile)
{
  PdbFormat pdb;
  Molecule molecule;
  bool success = pdb.readFile(std::string(AVOGADRO_DATA) + "/data/5ujw.pdb",
                              molecule);
  EXPECT_TRUE(success);
}