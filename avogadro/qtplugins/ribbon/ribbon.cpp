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

#include <cmath>

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
  Vector3f temp1;
  temp1 << 0, 0, 0;
  Vector3f temp2;
  temp2 << 0, 0, 0;
  Vector3f tempc;
  tempc<< 0, 0, 0;

  for(size_t i = 0; i < (molecule.atomCount()-2); i+=2)
  {
    Vector3f ca1 = molecule.atom(i).position3d().cast<float>();  // First alpha carbon
    Vector3f ca2 = molecule.atom(i+2).position3d().cast<float>();  // Second alpha carbn
    Vector3f o = molecule.atom(i+1).position3d().cast<float>();  //Carbonyl oxygen between ca1 & ca2
    Vector3f a = ca2 - ca1; // In the direction of the ribbon
    Vector3f b = o - ca1;
    Vector3f c; // Normal to peptide plane

    //Cross product
    c(0, 0) = a(1, 0)*b(2, 0) - a(2, 0)*b(1, 0);
    c(1, 0) = a(2, 0)*b(0, 0) - a(0, 0)*b(2, 0);
    c(2, 0) = a(0, 0)*b(1, 0) - a(1, 0)*b(0, 0);
    
    c = c/(sqrt(pow(c(0, 0), 2) + pow(c(1, 0), 2) + pow(c(2, 0), 2)));  //Normalization
    /*
     *NEED TO USE BETTER FUNCTIONS, PROBABLY FROM EIGEN
    */

    Vector3f d;  // Perpendicular to vector a in peptide plane
    d(0, 0) = c(1, 0)*a(2, 0) - c(2, 0)*a(1, 0);
    d(1, 0) = c(2, 0)*a(0, 0) - c(0, 0)*a(2, 0);
    d(2, 0) = c(0, 0)*a(1, 0) - c(1, 0)*a(0, 0);
    d = d/(sqrt(pow(d(0, 0), 2) + pow(d(1, 0), 2) + pow(d(2, 0), 2)));  //Normalization
    /*
     *NEED TO USE BETTER FUNCTIONS
    */

    Vector3f p = (ca1 + ca2)/2; // Midpoint of ca1 and ca2
    /*
     *Need to Translate by 1.5 Angstroms for reasonable helix diameter
     *in the direction of vector c, if residue in helix
    */

    d = d*3/2;  // d to be multiplied by half of required ribbon width
    // TODO: Use variable instead of 3

    Vector3f p1 = p - d; // p1 and p2 are end points of ribbon width
    Vector3f p2 = p + d;

    if( i != 0){
      Core::Array<Vector3f> vertices;
      vertices.push_back(temp1);
      vertices.push_back(temp2);
      vertices.push_back(p1);
      vertices.push_back(p2);

      Core::Array<Vector3f> normals;  // Should be a good approx of vertex normal
      normals.push_back(tempc);
      normals.push_back(tempc);
      normals.push_back(c);
      normals.push_back(c);

      unsigned char opacity = 255;

      GeometryNode* geometry = new GeometryNode;
      node.addChild(geometry);

      MeshGeometry* mesh1 = new MeshGeometry;
      geometry->addDrawable(mesh1);
      mesh1->setColor(Vector3ub(255, 0, 0));
      mesh1->setOpacity(opacity);
      unsigned int index1 = mesh1->addVertices(vertices, normals);
      Core::Array<unsigned int> indices;
      for(int k = 0; k < 4; ++k)
        indices.push_back(index1 + k);
      mesh1->addTriangles(indices);
      //mesh1->setRenderPass(opacity == 255 ? Rendering::OpaquePass
        //                                  : Rendering::TranslucentPass);
    }

    temp1 = p1;
    temp2 = p2;
    tempc = c;
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