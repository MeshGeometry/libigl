// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Zhongshi Jiang <jiangzs@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_EXACT_GEODESIC_H
#define IGL_EXACT_GEODESIC_H

#include <vector>
#include <variant>

#include "igl_inline.h"
#include <Eigen/Core>

namespace igl 
{
  // Exact geodesic algorithm for triangular mesh with the implementation from https://code.google.com/archive/p/geodesic/, 
  // and the algorithm first described by Mitchell, Mount and Papadimitriou in 1987
  //
  // Inputs:
  //   V  #V by 3 list of 3D vertex positions
  //   F  #F by 3 list of mesh faces
  //   VS #VS by 1 vector specifying indices of source vertices
  //   FS #FS by 1 vector specifying indices of source faces
  //   VT #VT by 1 vector specifying indices of target vertices
  //   FT #FT by 1 vector specifying indices of target faces
  // Output:
  //   D  #VT+#FT by 1 vector of geodesic distances of each target w.r.t. the nearest one in the source set
  //
  // Note: 
  //      Specifying a face as target/source means its center. 
  //
    template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedVS,
    typename DerivedFS,
    typename DerivedVT,
    typename DerivedFT,
    typename DerivedD>
    IGL_INLINE void exact_geodesic(
      const Eigen::MatrixBase<DerivedV> &V,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedVS> &VS,
      const Eigen::MatrixBase<DerivedFS> &FS,
      const Eigen::MatrixBase<DerivedVT> &VT,
      const Eigen::MatrixBase<DerivedFT> &FT,
      Eigen::PlainObjectBase<DerivedD> &D);

    // Inputs:
    //   V  #V by 3 list of 3D vertex positions
    //   F  #F by 3 list of mesh faces
    //   PS 3D vector specifying source point on the surface of the input mesh
    //   PT 3D vector specifying target point on the surface of the input mesh
    // Output:
    //   P  #P by 3 list of points giveing the piecewise linear geodesic connecting PS and PT
    struct VertexPoint {
      int vi;
    };
    struct FacePoint {
      int fi;
    };
    struct EdgePoint {
      int e0, e1;
    };
    using PointType = std::variant<VertexPoint, FacePoint, EdgePoint>;
    template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedPS,
    typename DerivedPT,
    typename DerivedP>
    IGL_INLINE void exact_geodesic_path(
      const Eigen::MatrixBase<DerivedV> &V,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedPS> &PS,
      const PointType Stype,
      const Eigen::MatrixBase<DerivedPT> &PT,
      const PointType Ttype,
      Eigen::PlainObjectBase<DerivedP> &P,
      std::vector<igl::PointType> &Ptypes);
}

#ifndef IGL_STATIC_LIBRARY
#  include "exact_geodesic.cpp"
#endif

#endif