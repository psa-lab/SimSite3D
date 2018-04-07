/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the SimSite3D Software project.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

/*
 * $Source: /psa/share/repository/pfizer_proj/src/search/point_and_surf_score.C,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 21:30:14 $
 * 
 * $Log: not supported by cvs2svn $
 * 
 */

#include <iostream>
#include <math_basics.H>
#include <point_and_surf_score.H>

using namespace SimSite3D;

const std::string point_and_surf_score::A_fname = "point_and_surf_score.C";

#if 0
// Rounded to 6 non zero digits
const my_float_t point_and_surf_score::CONSTANT_TERM = -1.07735;
const my_float_t point_and_surf_score::POLAR_SUM_W = 0.0;
const my_float_t point_and_surf_score::POLAR_MISMATCH_W = 0.0;
const my_float_t point_and_surf_score::DD_OR_AA_POLAR_SUM_W = -0.0250271;
const my_float_t point_and_surf_score::DONEPTOR_POLAR_SUM_W = 0.00804633;
const my_float_t point_and_surf_score::HPHOBIC_COUNT_W = 0.0;
const my_float_t point_and_surf_score::UNSAT_POLAR_W = 0.0;
const my_float_t point_and_surf_score::SURFACE_RMSE_W = 0.774077;
#endif

#if 1
// "new" weights -- June 2010
// Rounded to 6 non zero digits
const my_float_t point_and_surf_score::CONSTANT_TERM = -2.62899;
const my_float_t point_and_surf_score::POLAR_SUM_W = 0.0;
const my_float_t point_and_surf_score::POLAR_MISMATCH_W = 0.0;
const my_float_t point_and_surf_score::DD_OR_AA_POLAR_SUM_W = -0.0122689;
const my_float_t point_and_surf_score::DONEPTOR_POLAR_SUM_W = -0.00619202;
const my_float_t point_and_surf_score::HPHOBIC_COUNT_W = 0.0;
const my_float_t point_and_surf_score::UNSAT_POLAR_W = 0.0;
const my_float_t point_and_surf_score::SURFACE_RMSE_W = 2.11849;
#endif

const my_float_t point_and_surf_score::SCALED_CONSTANT_TERM = -2.70337;
const my_float_t point_and_surf_score::SCALED_POLAR_SUM_W = 0.0;
const my_float_t point_and_surf_score::SCALED_POLAR_MISMATCH_W = 0.0;
const my_float_t point_and_surf_score::SCALED_DD_OR_AA_POLAR_SUM_W = -0.237545;
const my_float_t point_and_surf_score::SCALED_DONEPTOR_POLAR_SUM_W = 0.0471545;
const my_float_t point_and_surf_score::SCALED_HPHOBIC_COUNT_W = 0.0;
const my_float_t point_and_surf_score::SCALED_UNSAT_POLAR_W = 0.0;
const my_float_t point_and_surf_score::SCALED_SURFACE_RMSE_W = 3.25161;

const int point_and_surf_score::A_num_terms = 6;

point_and_surf_score::point_and_surf_score()
{
  A_max_surf_pt_dist = 1.5;
  A_max_num_verts = 0;
  A_surf_closest_pts = 0;
  A_dists = 0;
  A_best_normals = 0;
  A_norm_dists = 0;
}

point_and_surf_score::~point_and_surf_score()
{
  if(A_surf_closest_pts) delete [] A_surf_closest_pts;
  A_surf_closest_pts = 0;
  A_dists = 0;
  A_best_normals = 0;
  A_norm_dists = 0;
}

#if 0
my_float_t 
point_and_surf_score::score(const ModelSitemap &model, 
                            const DbaseSitemap &dbase, rigid_align_vi scores)
{
  WeightedSumsScore::score(model, dbase, scores);
  //my_float_t pt_score = WeightedSumsScore::score(model, dbase, scores);

  std::vector<my_float_t> &terms = scores->terms;
  uint start_idx = terms.size();
  terms.resize(terms.size() + A_num_terms);
  std::fill(terms.end() - A_num_terms, terms.end(), 0.0);

  my_float_t &closest_polar_sum = terms[0];
  my_float_t &polar_mismatch_sum = terms[1];
  my_float_t &AA_DD_sum = terms[2];
  my_float_t &doneptor_sum = terms[3];
  my_float_t &phobic_count = terms[4];
  my_float_t &unsat_polar_count = terms[5]; 

  // Surface complementarity
  my_float_t &num_surf_pts = terms[start_idx];
  my_float_t &num_faces = terms[start_idx + 1];
  my_float_t &RMSE = terms[start_idx + 2];
  my_float_t &area = terms[start_idx + 3];
  my_float_t &RMS_norm_err = terms[start_idx + 4];

  const geometry::TransformableTrimesh& model_surf = 
    model.binding_site_mesh_handle();
  const geometry::SimpleTrimeshTwo& dbase_surf = 
    dbase.binding_site_mesh_handle();
  A_max_surf_pt_dist = dbase.max_corr_surf_pt_dist();
  init_storage(model_surf);

  dbase_surf.compare(model_surf.vertices_begin(), 
                     model_surf.number_of_vertices(), A_surf_closest_pts, 
                     A_dists, model_surf.normals_begin(), A_best_normals, 
                     A_norm_dists, A_max_surf_pt_dist);
//  model_surf.compute_features(A_dists, A_max_surf_pt_dist, A_norm_dists, 
  model_surf.compute_features(A_dists, 1.5, A_norm_dists, 
                              &num_surf_pts, &num_faces, &RMSE, &RMS_norm_err,
                              &area);
  if(scale_terms()){
    // The terms before start_idx were scaled in the base class
    std::vector<my_float_t>::const_iterator max_val =
      model.get_max_feature_vals().begin() + start_idx;
    std::vector<my_float_t>::iterator t = terms.begin() + start_idx;
    for( ; t < terms.end(); ++t, ++max_val) *t /= *max_val;

    return SCALED_CONSTANT_TERM + 
           SCALED_POLAR_SUM_W * closest_polar_sum + 
           SCALED_POLAR_MISMATCH_W * polar_mismatch_sum + 
           SCALED_DD_OR_AA_POLAR_SUM_W * AA_DD_sum +
           SCALED_DONEPTOR_POLAR_SUM_W * doneptor_sum +
           SCALED_HPHOBIC_COUNT_W * phobic_count + 
           SCALED_UNSAT_POLAR_W * unsat_polar_count +
           SCALED_SURFACE_RMSE_W * RMSE;

  }else{
    return CONSTANT_TERM + 
           POLAR_SUM_W * closest_polar_sum + 
           POLAR_MISMATCH_W * polar_mismatch_sum + 
           DD_OR_AA_POLAR_SUM_W * AA_DD_sum +
           DONEPTOR_POLAR_SUM_W * doneptor_sum +
           HPHOBIC_COUNT_W * phobic_count + 
           UNSAT_POLAR_W * unsat_polar_count +
           SURFACE_RMSE_W * RMSE;
  }
}
#endif

my_float_t
point_and_surf_score::correspondences(const ModelSitemap &model, 
                                      my_float_t **vertices_ptr, 
                                      my_float_t **closest_pts_ptr, 
                                      size_t *npts)
{
  const geometry::TransformableTrimesh& model_surf = 
    model.binding_site_mesh_handle();
  const int N = model_surf.number_of_vertices();
  *vertices_ptr = new my_float_t[3*N];
  *closest_pts_ptr = new my_float_t[3*N];
  *npts = 0;

  my_float_t SE = 0.0; 
  my_float_t *vert_out = *vertices_ptr;
  my_float_t *cp_out = *closest_pts_ptr;
  const my_float_t *dist = A_dists;
  const my_float_t *saved_cp = A_surf_closest_pts;
  const my_float_t *vert = model_surf.vertices_begin();
  const my_float_t *verts_end = model_surf.vertices_begin() + 3*N; 
  for( ; vert < verts_end; vert += 3, saved_cp += 3, ++dist){
    if(*dist <= A_max_surf_pt_dist){
      std::copy(vert, vert + 3, vert_out);
      std::copy(saved_cp, saved_cp + 3, cp_out);
      vert_out += 3;
      cp_out += 3;
      SE += (*dist)*(*dist);
      *npts += 1;
    }else{
      SE += A_max_surf_pt_dist*A_max_surf_pt_dist;
    }
  }

  return std::sqrt(SE/N);
}

void
point_and_surf_score::correspondences(const ModelSitemap &model, 
                                      std::vector<const my_float_t *> *q_pts,
                                      std::vector<const my_float_t *> *corr_pts)
{
  const geometry::TransformableTrimesh& model_surf = 
    model.binding_site_mesh_handle();
  const int N = model_surf.number_of_vertices();
  q_pts->resize(N);
  corr_pts->resize(N);

  //my_float_t SE = 0.0; 
  //my_float_t *vert_out = *vertices_ptr;
  //my_float_t *cp_out = *closest_pts_ptr;
  const my_float_t *dist = A_dists;
  const my_float_t *saved_cp = A_surf_closest_pts;
  const my_float_t *vert = model_surf.vertices_begin();
  const my_float_t *verts_end = model_surf.vertices_begin() + 3*N; 
  for(int i = 0; vert < verts_end; vert += 3, saved_cp += 3, ++dist, ++i){
    if(*dist <= A_max_surf_pt_dist){
      //std::cout << *dist << " " << vert << " " << saved_cp << std::endl;
      (*q_pts)[i] = vert;
      (*corr_pts)[i] =  saved_cp;
    }else{
      (*q_pts)[i] = 0;
      (*corr_pts)[i] = 0;
    }
  }
}

#if 0
inline void
point_and_surf_score::init_storage(const geometry::TransformableTrimesh& 
                                   model_surf)
{
  const int N = model_surf.number_of_vertices();
  if(N <= A_max_num_verts) return;

  if(A_surf_closest_pts){
    delete[] (A_surf_closest_pts);
   A_surf_closest_pts = 0;
   A_dists = 0;
  }

  A_max_num_verts = N;
  A_surf_closest_pts = new my_float_t[8 * N];
  A_dists = &A_surf_closest_pts[3*N];
  A_best_normals = &A_surf_closest_pts[4*N];
  A_norm_dists = &A_surf_closest_pts[7*N];
}
#endif
