// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Stefan Brugger <stefanbrugger@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "adjacency_list.h"
#include "boundary_loop.h"
#include "oriented_facets.h"
#include "slice.h"
#include "triangle_triangle_adjacency.h"
#include "vertex_triangle_adjacency.h"
#include "is_border_vertex.h"
#include <map>
#include <set>
#include <utility>

#include <iostream>

template <typename DerivedF, typename Index>
IGL_INLINE void igl::boundary_loop(
    const Eigen::MatrixBase<DerivedF> & F,
    std::vector<std::vector<Index> >& L)
{
  using namespace std;
  using namespace Eigen;

  if(F.rows() == 0)
    return;

  VectorXd Vdummy(F.maxCoeff()+1,1);
  DerivedF TT,TTi;
  vector<std::vector<int> > VF, VFi;
  triangle_triangle_adjacency(F,TT,TTi);
  vertex_triangle_adjacency(Vdummy,F,VF,VFi);

  vector<vector<int>> A;
  igl::adjacency_list(F, A);
  MatrixXi E;
  igl::oriented_facets(F, E);
  // Count edges
  map<pair<int,int>,int> edge_repeats;
  auto unordered_make_pair = [](auto a, auto b){ return a < b ? make_pair(a,b) : make_pair(b, a);};
  for (int i = 0; i < E.rows(); i++)
  {
    auto key = unordered_make_pair(E(i, 0), E(i,1));
    if(edge_repeats.count(key) == 0) {
      edge_repeats[key] = 1;
    } else {
      edge_repeats[key]++;
    }
  }
  

  vector<bool> unvisited = is_border_vertex(Vdummy,F);
  set<int> unseen;
  for (size_t i = 0; i < unvisited.size(); ++i)
  {
    if (unvisited[i])
      unseen.insert(unseen.end(),i);
  }
  

  // Stack based approach

  while (!unseen.empty())
  {
    vector<vector<int>> all_paths; // We will take the longest path. This could explode memory we'll see..
    vector<int> path;
    int start = *unseen.begin();
    path.push_back(start);
    unseen.erase(unseen.begin());
    unvisited[start] = false;

    bool finished_loop = false;
    while(true) {
      int v = path.back();
      // if(v == (9528) || v == (9527) || v == (32212) || v == 32213 ) {
      //     std::cout << "cur v: " << v << '\n';
      // }
      //find next
      int next = -1;
      for (int i = 0; i < A[v].size(); i++)
      {
        // if(v == (9528) || v == (9527) || v == (32212) || v == 32213 ) {
        //   std::cout << "A[v][" << i << "]: " << A[v][i] << '\n';
        // }
        // Only consider traversing edges that are on boundary by looking at half-edge count
        if(unvisited[A[v][i]] && edge_repeats[unordered_make_pair(v, A[v][i])] == 1) {
          // if(v == (9528) || v == (9527) || v == (32212) || v == 32213 ) {
          //   std::cout << "found next: " << A[v][i] << '\n';
          // }
          next = A[v][i];
          break;
        } else {
          // if(v == (9528) || v == (9527) || v == (32212) || v == 32213 ) {
          //   std::cout << "unvisited[A[v][i]]: " << unvisited[A[v][i]] <<  ", edge_repeats[unordered_make_pair(v, A[v][i])]: " << edge_repeats[unordered_make_pair(v, A[v][i])] << '\n';
          // }
        }
      }

      if(next == start) { // Found a loop
        break;
      } else if(next != -1) { // Keep following path
        path.push_back(next);
        unvisited[next] = false;
        unseen.erase(next);
      } else { // Back up
        // if(v == (9528) || v == (9527) || v == (32212) || v == 32213 ) {
        //   std::cout << "popping" << '\n';
        // }
        all_paths.push_back(path);
        if(path.size() == 1) {
          break;
        } else {
          path.pop_back();
        }
      }
    }
    auto& longest_path = *std::max_element(all_paths.begin(), all_paths.end(), [](auto& a, auto& b){return a.size() < b.size();});
    // for(auto& p : longest_path) {
    //   std::cout << p << ", ";
    // }
    //   std::cout << "\n";
    L.push_back(longest_path);
  }

  while (!unseen.empty())
  {
    vector<Index> l;

    // Get first vertex of loop
    int start = *unseen.begin();
    unseen.erase(unseen.begin());
    unvisited[start] = false;
    l.push_back(start);

    bool done = false;
    while (!done)
    {
      // Find next vertex
      bool newBndEdge = false;
      int v = l[l.size()-1];
      int next;
      
      // My Way
      // for (int i = 0; i < A[v].size(); i++)
      // {
      //   // Only consider traversing edges that are on boundary by looking at half-edge count
      //   if(unvisited[A[v][i]] && edge_repeats[unordered_make_pair(v, A[v][i])] == 1) {
      //     next = A[v][i];
      //     newBndEdge = true;
      //     break; // Assuming only one incident border vertex
      //   }
      // }
      
      // Old way
      // for (int i = 0; i < (int)VF[v].size() && !newBndEdge; i++)
      // {
      //   int fid = VF[v][i];

      //   if (TT.row(fid).minCoeff() < 0.) // Face contains boundary edge
      //   {
      //     int vLoc = -1;
      //     if (F(fid,0) == v) vLoc = 0;
      //     if (F(fid,1) == v) vLoc = 1;
      //     if (F(fid,2) == v) vLoc = 2;

      //     int vNext = F(fid,(vLoc + 1) % F.cols());

      //     newBndEdge = false;
      //     if (unvisited[vNext] && TT(fid,vLoc) < 0)
      //     {
      //       next = vNext;
      //       newBndEdge = true;
      //     }
      //   }
      // }

      if (newBndEdge)
      {
        l.push_back(next);
        unseen.erase(next);
        unvisited[next] = false;
      }
      else
        done = true;
    }
    L.push_back(l);
  }
}

template <typename DerivedF, typename Index>
IGL_INLINE void igl::boundary_loop(
  const Eigen::MatrixBase<DerivedF>& F,
  std::vector<Index>& L)
{
  using namespace Eigen;
  using namespace std;

  if(F.rows() == 0)
    return;

  vector<vector<int> > Lall;
  boundary_loop(F,Lall);

  int idxMax = -1;
  size_t maxLen = 0;
  for (size_t i = 0; i < Lall.size(); ++i)
  {
    if (Lall[i].size() > maxLen)
    {
      maxLen = Lall[i].size();
      idxMax = i;
    }
  }

  //Check for meshes without boundary
  if (idxMax == -1)
  {
      L.clear();
      return;
  }

  L.resize(Lall[idxMax].size());
  for (size_t i = 0; i < Lall[idxMax].size(); ++i)
  {
    L[i] = Lall[idxMax][i];
  }
}

template <typename DerivedF, typename DerivedL>
IGL_INLINE void igl::boundary_loop(
  const Eigen::MatrixBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedL>& L)
{
  using namespace Eigen;
  using namespace std;

  if(F.rows() == 0)
    return;

  vector<int> Lvec;
  boundary_loop(F,Lvec);

  L.resize(Lvec.size());
  for (size_t i = 0; i < Lvec.size(); ++i)
    L(i) = Lvec[i];
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::boundary_loop<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::boundary_loop<Eigen::Matrix<int, -1, -1, 0, -1, -1>, int>(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >&);
#endif
