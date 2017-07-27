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

Ribbon::Ribbon(QObject* p)
  : ScenePlugin(p), m_enabled(true), m_group(nullptr)
{
}
Ribbon::~Ribbon()
{
}

void Ribbon::process(const Molecule& molecule, Rendering::GroupNode& node)
{
  m_group = &node;
  GeometryNode* geometry = new GeometryNode;
  node.addChild(geometry);

  LineStripGeometry* lines = new LineStripGeometry;
  lines->identifier().molecule = &molecule;
  lines->identifier().type = Rendering::BondType;
  geometry->addDrawable(lines);
  
  vector<Core::Atom> ca; // Contains all alpha carbons
  vector<Core::Atom> o;   //Contains all carbonyl oxygens

  for(Index i = 0; i < molecule.bondCount(); ++i) {
    Core::Atom atom = molecule.atom(i);
    if (atom.atomName() == "CA")
      ca.push_back(atom);

    else if(atom.atomName() == "O")
      o.push_back(atom);
  }

  for (size_t s = 0; s < ca.size(); ++s)
  {
    Vector3f pos1 = ca[s].position3d().cast<float>();
    Vector3f pos2 = ca[s+1].position3d().cast<float>();
    Vector3f pos3 = o[s].position3d().cast<float>();

    Vector3ub color1(Elements::color(ca[s].atomicNumber()));
    Vector3ub color2(Elements::color(ca[s+1].atomicNumber()));
    Vector3ub color3(Elements::color(o[s].atomicNumber()));

    Array<Vector3f> points;
    Array<Vector3ub> colors;

    points.push_back(pos1);
    points.push_back(pos2);
    points.push_back(pos3);

    colors.push_back(color1);
    colors.push_back(color2);
    colors.push_back(color3);

    lines->addLineStrip(points, colors, 1.0f);
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