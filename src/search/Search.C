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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <basics.H>
#include <param_tools.H>
#include <mat_ops.H>
#include <Search.H>

using namespace SimSite3D;
const std::string Search::A_fname = "Search.C";

Search::Search(const SearchParameters* args_in)
{
  init_vars();
  if(args_in->status() != BaseParameters::READY) return;
  if(!setup_directories(args_in->proj_output, 0770)) return;
  A_args = args_in;

  std::string path;
  std::string struct_id;
  get_path_and_struct_id(args_in->model_file_name, &path, &struct_id);
 
  if(A_args->score_str == "sample_all_triangles")
    A_model = new AllPairsSitemapTest(
      path, struct_id, *args_in, args_in->normalize, 
      args_in->use_hbond_surfaces, args_in->allow_hphob_triangles);
  else 
    A_model = new ModelSitemap(
      path, struct_id, *args_in, args_in->max_corr_surf_pt_dist, 
      args_in->normalize, args_in->use_hbond_surfaces, 
      args_in->allow_hphob_triangles);
  if(A_model == 0 || A_model->fail()){
    err_msg(A_fname, "Cstr()", "Failed to intialize the model sitemap");
    return;
  }
  my_float_t MAX_ATOM_DIST_TO_CHECK_VDW_OVERLAP = 4.5;
  A_model->bind_site_atoms().bin_coordinates(MAX_ATOM_DIST_TO_CHECK_VDW_OVERLAP);

  // If we are supposed generate random alignments, use clock to seed the 
  // random number generator.  In the future, if we really care about random
  // numbers, the merseinne twister should be used (or similar quality rng).
  if(A_args->num_rand_aligns) std::srand(std::time(0));

  A_fail = false;
}

Search::~Search()
{
  if(A_model) delete A_model;
  init_vars();
}

void
Search::init_vars()
{
  A_fail = true;
  A_model = NULL;
}

bool
Search::run()
{
  if(A_fail){
    warn(A_fname, "run", "Cannot run because an error has occured");
    return false;
  }

  if(A_args->time_process) start_timer();
  align_method_t align_method;
  if(!set_alignment_method(&align_method)) return false;

  if((A_args->do_IK || A_args->save_rigid_scores) &&
     A_args->do_internal_prot_lig_score){
    warn(A_fname, "run", "Currently the saving of rigid scores or ArtSurf is mutually exclusive with protein-ligand scoring.");
    return false;
  }

  if(A_args->do_IK){
    std::vector<align_w_preIK_data_t> alignments;
    if(!run(&alignments, align_method)) return false;
  }else if(A_args->save_rigid_scores){
    std::vector<align_w_rigid_data_t> alignments;
    if(!run(&alignments, align_method)) return false;
  }else if(A_args->do_internal_prot_lig_score){
    std::vector<align_w_pl_score_t> alignments;
    if(!run(&alignments, align_method)) return false;
  }else{
    std::vector<rigid_align_t> alignments;
    if(!run(&alignments, align_method)) return false;
  }

  get_timer_and_write_to_file();
  return true;
}

bool
Search::get_timer_and_write_to_file()
{
  if(!A_args->time_process) return true;

  if(A_args->time_process){ 
    double real, virt, prof;
    A_timer.get(&real, &virt, &prof);
    std::cout << "real time: " << real << "\n";
    std::cout << "virt time: " << virt << "\n";
    std::cout << "prof time: " << prof << "\n";
  }

  std::ifstream results_in;
  if(!open_ifstream(results_in, A_args->ofname)) return false;
  std::vector<std::string> lines;
  int num_param_lines = 12;
  std::string line;
  for(int i = 0; i < num_param_lines && std::getline(results_in, line); ++i)
    lines.push_back(line + "\n");

  double real, virt, prof;
  A_timer.get(&real, &virt, &prof);
  std::ostringstream ostr;
  ostr << std::left << std::fixed;
  ostr.precision(2);
  ostr << "# SimSite3D timing statistics:\n";
  ostr << std::setw(50) << "#   Wall clock time:" << real << " sec.\n";
  ostr << std::setw(50) << "#   CPU time:" << prof << " sec.\n";
  ostr << std::setw(50) << "#   User time:" << virt << " sec.\n";
  ostr << std::setw(50) << "#   Kernel time:" << prof - virt << " sec.\n";
  lines.push_back(ostr.str()); 
  while(std::getline(results_in, line)) lines.push_back(line + "\n");
  results_in.close();

  std::ofstream results_out;
  if(!open_ofstream(results_out, A_args->ofname)) return false;
  for(size_t i = 0; i < lines.size(); ++i)
    results_out << lines[i];
  results_out.close();


  return true;
}

bool
Search::setup_directories(std::string base_dir, mode_t mode)
{
  if(!dir_exists(base_dir, false) && !my_mkdir(base_dir.c_str(), mode)) 
    return false;
  std::string dir = base_dir + "/moved_ligands";
  if(!dir_exists(dir, false) && !my_mkdir(dir.c_str(), mode)) return false;
  dir = base_dir + "/ligand_fragments";
  if(!dir_exists(dir, false) && !my_mkdir(dir.c_str(), mode)) return false;
  return true;
}






#if 0
void
Search::add_identity_alignment()
{
  rigid_align_t my_align;
  std::fill(my_align.R, my_align.R + 9, 0.0);
  my_align.R[0] = my_align.R[4] = my_align.R[8] = 1.0;
  std::fill(my_align.T, my_align.T + 3, 0.0);
  my_float_t Q[] = {1.0, 0.0, 0.0, 0.0};
  my_align.Q = Quaternion(Q, 4);
  alignments.push_back(my_align);
}

void 
Search::add_random_alignments(const uint N, const my_float_t* C)
{
  const my_float_t half_RAND_MAX = 0.5 * RAND_MAX;
  alignments.reserve(6*N);

  std::vector<my_float_t> theta_steps(9);
  theta_steps[0] = 0.0;
  for(size_t i = 1; i < 9; ++i) theta_steps[i] = theta_steps[i-1] + 0.125;
   

  rigid_align_t my_align;
  for(size_t n = 0; n < N; ++n)
    for(size_t t = 1; t < theta_steps.size(); ++t)
      for(size_t d = 1; d < 6; ++d){
        my_align.Q.randomize(theta_steps[t-1]*M_PI, theta_steps[t]*M_PI);
        my_float_t R[9];
        my_align.Q.get_ortho_rot_mat(R);
        for(size_t i = 0; i < 3; ++i)
          for(size_t j = 0; j < 3; ++j)
            my_align.R[3*i + j] = R[3*j + i];
    
        std::copy(C, C + 3, my_align.T);
        my_gemv(3, 3, -1.0, my_align.R, 3, C, 1, 1.0, my_align.T, 1);
        for(size_t i = 0; i < 3; ++i)
          my_align.T[i] += d * ((std::rand() - half_RAND_MAX) / RAND_MAX);
        alignments.push_back(my_align);
      }
}

// This method could easily be incorrect -- verify versus the Python
// implementation 
void 
Search::add_poisson_alignments(const uint N, const my_float_t* C)
{
  const my_float_t half_RAND_MAX = 0.5 * RAND_MAX;
  const size_t nbins = 8;
  const my_float_t q_dist_lim = 1 - std::cos(M_PI/48.0);
  const my_float_t t_dist_lim = 0.05;
  const my_float_t t_dist_lim_squared = t_dist_lim*t_dist_lim;

  alignments.reserve(nbins*N);

  std::vector<my_float_t> q_bins(nbins+1);
  for(size_t i = 0; i <= nbins; ++i) 
    q_bins[i] = 1.0 - std::cos(i * M_PI_2 / nbins);

  Quaternion id_rot;
  rigid_align_t my_align;
  for(size_t bin_num = 1; bin_num <= nbins; ++bin_num)
    for(size_t n = 0; n < N; ++n){
      my_align.Q.randomize();
      my_float_t q_dist = distance(my_align.Q, id_rot);
      //while(q_dist < q_bins[bin_num - 1] || q_bins[bin_num] < q_dist){
      while(q_bins[bin_num] < q_dist){
        my_align.Q.randomize();
        q_dist = distance(my_align.Q, id_rot);
      }

      my_float_t tform[3];
      std::copy(C, C+3, tform);
      my_align.Q.get_ortho_rot_mat(my_align.R);
      my_gemm(1, 3, 3, -1.0, C, 3, my_align.R, 3, tform, 1, 1.0);
/*
      my_float_t R[9];
      my_align.Q.get_ortho_rot_mat(R);
      for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
          my_align.R[3*i + j] = R[3*j + i];
      my_gemv(3, 3, -1.0, my_align.R, 3, C, 1, 1.0, tform, 1);
*/

      for(bool too_close = true; too_close; ){     
        too_close = false;
        std::copy(tform, tform + 3, my_align.T);
        for(size_t i = 0; i < 3; ++i)
          my_align.T[i] += 
            0.125 * (bin_num - 1) * ((std::rand() - half_RAND_MAX) / RAND_MAX);

        std::vector<rigid_align_t>::iterator align_iter;
        if(bin_num == 1) align_iter = alignments.begin();
        else align_iter = alignments.begin() + N * (bin_num - 1);
        for( ; align_iter != alignments.end(); ++align_iter){
          q_dist = distance(my_align.Q, align_iter->Q);
          my_float_t square_t_dist = dist_squared(my_align.T, C);
          if(q_dist < q_dist_lim && square_t_dist < t_dist_lim_squared){
            too_close = true;
            break;
          }
        }
      }
      alignments.push_back(my_align); 
    }
  std::cout << "finished generating alignments\n";
}
#endif
