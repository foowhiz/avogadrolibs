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
#include <avogadro/core/mesh.h>
#include <avogadro/core/vector.h>
#include <avogadro/rendering/geometrynode.h>
#include <avogadro/rendering/groupnode.h>
#include <avogadro/rendering/meshgeometry.h>


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
using Core::Mesh;

using Rendering::GeometryNode;
using Rendering::GroupNode;
using Rendering::MeshGeometry;
using std::vector;

Ribbon::Ribbon(QObject* p)
  : ScenePlugin(p), m_enabled(false)
{
}
Ribbon::~Ribbon()
{
}

void Ribbon::process(const Molecule& molecule, Rendering::GroupNode& node)
{
  for(size_t i = 0; i < (molecule.atomCount()-2); i+=2)
  {
    //size_t i = 0;
    Core::Atom ca1 = molecule.atom(i);  // First alpha carbon
    Core::Atom ca2 = molecule.atom(i+2);  // Second alpha carbon
    Core::Atom o = molecule.atom(i+1);  //Carbonyl oxygen between ca1 & ca2
    Vector3f a = ca2.position3d().cast<float>() - ca1.position3d().cast<float>();
    Vector3f b = o.position3d().cast<float>() - ca1.position3d().cast<float>();
    Vector3f c;

    //Cross product
    c(0, 0) = a(1, 0)*b(2, 0) - a(2, 0)*b(1, 0);
    c(1, 0) = a(2, 0)*b(0, 0) - a(0, 0)*b(2, 0);
    c(2, 0) = a(0, 0)*b(1, 0) - a(1, 0)*b(0, 0);
    /*
     *Need to Normalize
    */

    Eigen::Matrix<Real, 3, 1> d;
    d(0, 0) = c(1, 0)*a(2, 0) - c(2, 0)*a(1, 0);
    d(1, 0) = c(2, 0)*a(0, 0) - c(0, 0)*a(2, 0);
    d(2, 0) = c(0, 0)*a(1, 0) - c(1, 0)*a(0, 0);
    /*
     *Need to Normalize
    */

    Core::Array<Vector3f> vertices;
    vertices.push_back(ca1.position3d().cast<float>());
    vertices.push_back(o.position3d().cast<float>());
    vertices.push_back(ca2.position3d().cast<float>());

    Core::Array<Vector3f> normals;  // Should be a good approx
    for (int j = 0; j < 3; ++j)
      normals.push_back(c);

    unsigned char opacity = 255;

    GeometryNode* geometry = new GeometryNode;
    node.addChild(geometry);

    MeshGeometry* mesh1 = new MeshGeometry;
    geometry->addDrawable(mesh1);
    mesh1->setColor(Vector3ub(255, 0, 0));
    mesh1->setOpacity(opacity);
    unsigned int index1 = mesh1->addVertices(vertices, normals);
    mesh1->addTriangle(index1, index1 + 1, index1 + 2);
    //mesh1->setRenderPass(opacity == 255 ? Rendering::OpaquePass
      //                                  : Rendering::TranslucentPass);
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