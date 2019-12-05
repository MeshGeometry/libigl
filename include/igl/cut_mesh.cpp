// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2019 Hanxiao Shen <hanxiao@cims.nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <igl/cut_mesh.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/HalfEdgeIterator.h>
#include <igl/is_border_vertex.h>

// wrapper for input/output style
template <typename DerivedV, typename DerivedF, typename DerivedC>
IGL_INLINE void igl::cut_mesh(
  const Eigen::MatrixBase<DerivedV>& V,
  const Eigen::MatrixBase<DerivedF>& F,
  const Eigen::MatrixBase<DerivedC>& C,
  Eigen::PlainObjectBase<DerivedV>& Vn,
  Eigen::PlainObjectBase<DerivedF>& Fn
){
  Vn = V;
  Fn = F;
  typedef typename DerivedF::Scalar Index;
  Eigen::Matrix<Index,Eigen::Dynamic,1> _I;
  cut_mesh(Vn,Fn,C,_I);
}

template <typename DerivedV, typename DerivedF, typename DerivedC, typename DerivedI>
IGL_INLINE void igl::cut_mesh(
  const Eigen::MatrixBase<DerivedV>& V,
  const Eigen::MatrixBase<DerivedF>& F,
  const Eigen::MatrixBase<DerivedC>& C,
  Eigen::PlainObjectBase<DerivedV>& Vn,
  Eigen::PlainObjectBase<DerivedF>& Fn,
  Eigen::PlainObjectBase<DerivedI>& I
){
  Vn = V;
  Fn = F;
  cut_mesh(Vn,Fn,C,I);
}

template <typename DerivedV, typename DerivedF, typename DerivedC, typename DerivedI>
IGL_INLINE void igl::cut_mesh(
  Eigen::PlainObjectBase<DerivedV>& V,
  Eigen::PlainObjectBase<DerivedF>& F,
  const Eigen::MatrixBase<DerivedC>& C,
  Eigen::PlainObjectBase<DerivedI>& I
){
  typedef typename DerivedF::Scalar Index;
  DerivedF FF, FFi;
  igl::triangle_triangle_adjacency(F,FF,FFi);
  // std::cout << "FF.rows(): " << FF.rows() << std::endl;
  // std::cout << "FF.cols(): " << FF.cols() << std::endl;

  igl::cut_mesh(V,F,FF,FFi,C,I);
}

template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedFFi, typename DerivedC, typename DerivedI>
IGL_INLINE void igl::cut_mesh(
  Eigen::PlainObjectBase<DerivedV>& V,
  Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::MatrixBase<DerivedFF>& FF,
  Eigen::MatrixBase<DerivedFFi>& FFi,
  const Eigen::MatrixBase<DerivedC>& C,
  Eigen::PlainObjectBase<DerivedI>& I
){

  typedef typename DerivedF::Scalar Index;
  // std::cout << "x" << std::endl;
  // store current number of occurance of each vertex as the alg proceed
  Eigen::Matrix<Index,Eigen::Dynamic,1> occurence(V.rows());
  occurence.setConstant(1);
  
  // std::cout << "x" << std::endl;
  // set eventual number of occurance of each vertex expected
  Eigen::Matrix<Index,Eigen::Dynamic,1> eventual(V.rows());
  eventual.setZero();
  for(Index i=0;i<F.rows();i++){
    for(Index k=0;k<3;k++){
      Index u = F(i,k);
      Index v = F(i,(k+1)%3);
      if(FF(i,k) == -1){ // add one extra occurance for boundary vertices
        eventual(u) += 1;
      }else if(C(i,k) == 1 && u < v){ // only compute every (undirected) edge ones
        eventual(u) += 1;
        eventual(v) += 1;
      }
    }
  }
  // std::cout << "x" << std::endl;
  
  // original number of vertices
  Index n_v = V.rows(); 
  
  // estimate number of new vertices and resize V
  Index n_new = 0;
  for(Index i=0;i<eventual.rows();i++)
    n_new += ((eventual(i) > 0) ? eventual(i)-1 : 0);
  V.conservativeResize(n_v+n_new,Eigen::NoChange);
  I = DerivedI::LinSpaced(V.rows(),0,V.rows());
  
  // std::cout << "x" << std::endl;
  // pointing to the current bottom of V
  Index pos = n_v;
  for(Index f=0;f<C.rows();f++){
    // std::cout << "f " << f << ": ";
    for(Index k=0;k<3;k++){
      // std::cout << " k " << k;
      Index v0 = F(f,k);
      if(F(f,k) >= n_v) continue; // ignore new vertices
      // std::cout << "y" << std::endl;
      if(C(f,k) == 1 && occurence(v0) != eventual(v0)){
        // std::cout << "y" << std::endl;
        igl::HalfEdgeIterator<DerivedF,DerivedFF,DerivedFFi> he(F,FF,FFi,f,k);
        // std::cout << "first" << std::endl;
        // std::cout << "he.FF.rows(): " << he.FF.rows() << std::endl;
        // std::cout << "he.FF.cols(): " << he.FF.cols() << std::endl;
        Index fi = he.Fi();
        Index ei = he.Ei();

        // rotate clock-wise around v0 until hit another cut
        std::vector<Index> fan;

        // std::cout << "second" << std::endl;
        // std::cout << "he.FF.rows(): " << he.FF.rows() << std::endl;
        // std::cout << "he.FF.cols(): " << he.FF.cols() << std::endl;
        // std::cout << "y" << std::endl;
        do{
          fan.push_back(fi);
          // std::cout << "1" << std::endl;
          he.flipE();
          // std::cout << "1.5" << std::endl;
        // std::cout << "FF.rows(): " << FF.rows() << std::endl;
        // std::cout << "FF.cols(): " << FF.cols() << std::endl;
        // std::cout << "he.FF.rows(): " << he.FF.rows() << std::endl;
        // std::cout << "he.FF.cols(): " << he.FF.cols() << std::endl;
          he.flipF();
          // std::cout << "2" << std::endl;
          fi = he.Fi();
          ei = he.Ei();
          // std::cout << "3" << std::endl;
        }while(C(fi,ei) == 0 && !he.isBorder());
        
        // std::cout << "z" << std::endl;
        // std::cout << "z" << std::endl;
        // make a copy
        V.row(pos) << V.row(v0);
        I(pos) = v0;
        // std::cout << "z" << std::endl;
        // add one occurance to v0
        occurence(v0) += 1;
        // std::cout << "z" << std::endl;
        
        // std::cout << "y" << std::endl;
        // replace old v0
        for(Index f0: fan)
          for(Index j=0;j<3;j++)
            if(F(f0,j) == v0)
              F(f0,j) = pos;
        
        // std::cout << "y" << std::endl;
        // mark cuts as boundary
        FF(f,k) = -1;
        FF(fi,ei) = -1;
        
        pos++;
      }
    }
    // std::cout << '\n';
  }
  // std::cout << "x" << std::endl;
  
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::cut_mesh<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&);
template void igl::cut_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::cut_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::cut_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif