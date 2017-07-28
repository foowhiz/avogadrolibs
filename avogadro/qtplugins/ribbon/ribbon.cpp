/******************************************************************************

  This source file is part of the Avogadro project.

  Copyright 2014 Kitware, Inc.

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

******************************************************************************/

#include "ribbon.h"

#include <avogadro/core/elements.h>
#include <avogadro/core/molecule.h>
#include <avogadro/rendering/geometrynode.h>
#include <avogadro/rendering/groupnode.h>
#include <avogadro/rendering/linestripgeometry.h>

#include <Eigen/Dense>

#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

namespace Avogadro {
namespace QtPlugins {

using Core::Elements;
using Core::Molecule;
using Core::Array;
using Rendering::GeometryNode;
using Rendering::GroupNode;
using Rendering::LineStripGeometry;
using std::vector;
using namespace Eigen;

Ribbon::Ribbon(QObject* p)
  : ScenePlugin(p), m_enabled(true), m_group(nullptr)
{
}
Ribbon::~Ribbon()
{
}

void Ribbon::process(const Molecule& molecule, Rendering::GroupNode& node)
{
  vector<Core::Atom> ca; // Contains all alpha carbons
  vector<Core::Atom> o;   //Contains all carbonyl oxygens

  for(Index i = 0; i < molecule.atomCount(); ++i) {
    Core::Atom atom = molecule.atom(i);
    if (atom.getAtomName() == "CA")
      ca.push_back(atom);

    else if(atom.getAtomName() == "O")
      o.push_back(atom);
    }

  for(size_t t; t < o.size(); ++t)
  {
    Eigen::Matrix<Real, 3, 1> a = ca[t+1].position3d() - ca[t].position3d();
    Eigen::Matrix<Real, 3, 1> b = o[t].position3d() - ca[t].position3d();
    Eigen::Matrix<Real, 3, 1> c;

    //Cross product
    c(0, 0) = a(1, 0)*b(2, 0) - a(2, 0)*b(1, 0);
    c(1, 0) = a(2, 0)*b(0, 0) - a(0, 0)*b(2, 0);
    c(2, 0) = a(0, 0)*b(1, 0) - a(1, 0)*b(0, 0);

    Eigen::Matrix<Real, 3, 1> d;
    d(0, 0) = c(1, 0)*a(2, 0) - c(2, 0)*a(1, 0);
    d(1, 0) = c(2, 0)*a(0, 0) - c(0, 0)*a(2, 0);
    d(2, 0) = c(0, 0)*a(1, 0) - c(1, 0)*a(0, 0);

  }
}

bool Ribbon::isEnabled() const
{
  return m_enabled;
}

void Ribbon::setEnabled(bool enable)
{
  m_enabled = enable;
}

}
}