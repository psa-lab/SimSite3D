/******************************************************************************
 * Copyright (c) 2010, Michigan State University (MSU) Board of Trustees.
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
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/
#include <Search.H>
#include <NoTier1Score.H>
#include <HbondSurfacesScore.H>

using namespace SimSite3D;

template < class tier1_SF, class tier2_SF, typename align_T >
class GridSamplingAndScoring : 
  public ScoreRigidAlignments<tier1_SF, tier2_SF, align_T> {
public:

  typedef ScoreRigidAlignments<tier1_SF, tier2_SF, align_T>
    ScoreBase;
  
  typedef typename ScoreBase::rigid_align_vec           rigid_align_vec;
  typedef typename ScoreBase::rigid_align_vi            rigid_align_vi;
  typedef typename ScoreBase::tier1_score_cmp           tier1_score_cmp;
  typedef typename ScoreBase::tier1_score_cmp           tier2_score_cmp;
  typedef std::pair<my_float_t, align_T>		align_pair;
  typedef typename std::multimap<my_float_t, align_T, tier1_score_cmp>   
                                                        tier1_score_mmap;
  typedef typename std::multimap<my_float_t, align_T, tier2_score_cmp>   
                                                        tier2_score_mmap;
  typedef typename tier1_score_mmap::iterator		tier1_score_mmi;
  typedef typename tier2_score_mmap::iterator		tier2_score_mmi;


  GridSamplingAndScoring(ModelSitemap* model_in, const SearchParameters& params,
                       const my_float_t surf_pt_w, const my_float_t hb_cap_w)
    : ScoreRigidAlignments<tier1_SF, tier2_SF, align_T>(model_in, params)
  {
    //A_model = model_in;
    A_model = ScoreBase::model();
    A_fine_tune_best_tier2_alignment = params.fine_tune_best_tier2_alignment;
    A_fine_tune_tier2_alignments = params.fine_tune_tier2_alignments;
    A_mu = ScoreBase::mu();
    A_sigma = ScoreBase::sigma();
    A_save_rigid_results = params.save_rigid_scores;
    A_surf_pt_w = surf_pt_w;
    A_hb_cap_w = hb_cap_w;
    A_max_num_aligns = params.num_scores_to_keep;
    A_score_cutoff = params.score_cutoff;
    A_max_tier1_aligns = params.max_tier1_aligns;


//    A_rot_grid;
    std::ifstream g_file;
    if(!open_ifstream(g_file, "/home/vanvoor4/code/ISOI/SO3_grid/data.qua_res2")){
//    if(!open_ifstream(g_file, "/home/vanvoor4/code/ISOI/SO3_grid/test_data.txt")){
      std::cout << "unable to open the file /home/vanvoor4/code/ISOI/SO3_grid/data.qua_res2\n";
      return;
    }

#if 1
    for(std::string line; std::getline(g_file, line); ){
      my_float_t q[4];
      g_file >> q[0] >> q[1] >> q[2] >> q[3];
      A_rot_grid.push_back(Quaternion(q,4));
    }
#else
    my_float_t q[] = {1.0, 0.0, 0.0, 0.0};
    A_rot_grid.push_back(Quaternion(q,4));
#endif
    std::cout << "Loaded a quaternion grid with " << A_rot_grid.size() << " grid points\n";

    
  }

  bool
  score_alignments(rigid_align_vec& aligns, DbaseSitemap* dset_site,
                   std::ostream& results_out = std::cout,
                   std::ostream& rigid_results_out = std::cout)
  {
    // This one is a bit different as we will be generating orientations on the 
    // fly as they are based on a sampling grid.  This means that the map
    // will be storing align objects and not iterators to them as we don't want
    // to fill up memory with lots of aligns and we don't want invalid 
    // iterators (each time a vector is modified, existing iterators are 
    // made invalid).

    mol2File *dbase_ligand = 0;
    std::string dbase_struct_id, dbase_mol_id;
    if(!ScoreBase::load_ligand(dset_site->ligand_file_name(),
                               dset_site->atoms_file_name(),
                               &dbase_ligand, &dbase_struct_id, &dbase_mol_id))
      return false;

    // forced flush for debugging
    std::cout << "Comparing " << ScoreBase::model_struct_id() << " to "
              << dbase_struct_id << std::endl;

    const SitemapPointsFile &q_pts = ScoreBase::model()->sitemap_points();
    my_float_t q_centroid[3];
    q_pts.centroid_3D(q_centroid);
    tier1_score_mmap top_tier1_aligns;
    tier2_score_mmap top_tier2_aligns;
   
    atom_vci lig_atom = dbase_ligand->atoms_begin();
    size_t num_hvy = 0;
    for( ; lig_atom < dbase_ligand->atoms_end(); ++lig_atom)
      if(lig_atom->name != H) ++num_hvy;
    size_t cnt = 0; 
    std::cout << "number of dset lig hvy atoms: " << num_hvy << "\n"
              << "number of orientations: " 
              << A_rot_grid.size()*num_hvy << "\n";

//    const my_float_t stride = 1.0 / num_samps;
//    std::cout << "stride: " << stride << "\n";
    
    size_t atom_no = 0;

    
    lig_atom = dbase_ligand->atoms_begin();
    for( ; lig_atom < dbase_ligand->atoms_end(); ++lig_atom){
      if(lig_atom->name == H) continue;
      // could be abit confusing -- we are using "inverse" transform
      // on the query sites -- this implies we want to have the
      // translation from the dset to the query site as it will get reversed

      for(std::vector<Quaternion>::iterator Q = A_rot_grid.begin();
          Q < A_rot_grid.end(); ++Q){

        // get dset "centroid" and rotate it by R, then get translation
        // from the rotated "centroid" to q_centroid
        align_T tmp_align;      
        Q->get_ortho_rot_mat(tmp_align.R);
        std::copy(q_centroid, q_centroid + 3, tmp_align.T);
        my_gemm(1, 3, 3, -1.0, lig_atom->pos, 3, tmp_align.R, 3, tmp_align.T, 3, 1.0);

#if 0
        std::cout << "Dset centroid: " << lig_atom->pos[0] << " "
                  << lig_atom->pos[1] << " " << lig_atom->pos[2] << "\n";
        std::cout << "Q centroid: " << q_centroid[0] << " "
                  << q_centroid[1] << " " << q_centroid[2] << "\n";
        std::cout << "T: " << tmp_align.T[0] << " " << tmp_align.T[1] << " " << tmp_align.T[2] << "\n";

            std::cout << "R: " << tmp_align.R[0] << " " << tmp_align.R[1] << " " << tmp_align.R[2] << "\n"
                      << "   " << tmp_align.R[3] << " " << tmp_align.R[4] << " " << tmp_align.R[5] << "\n"
                      << "   " << tmp_align.R[6] << " " << tmp_align.R[7] << " " << tmp_align.R[8] << "\n";
#endif

        if(ScoreBase::A_tier1_score_class.score_is_noop()){ 
          // call blah with tier2
          blah(&tmp_align, ScoreBase::A_tier2_score_class,
                ScoreBase::A_tier2_score_cmp, &top_tier2_aligns,
                dbase_ligand, true, A_max_num_aligns, 
                A_score_cutoff, A_fine_tune_tier2_alignments,
                dset_site);
        }else{
          // call blah with tier1
          blah(&tmp_align, ScoreBase::A_tier1_score_class,
                ScoreBase::A_tier1_score_cmp, &top_tier1_aligns,
                dbase_ligand, true, A_max_tier1_aligns,
                A_score_cutoff, A_fine_tune_tier2_alignments,
                dset_site);
        }  
        ++cnt;
      }
      atom_no += 1;
    }

    // use tier2 sieve
    if(! ScoreBase::A_tier1_score_class.score_is_noop()){ 
      tier1_score_mmi tier1_iter = top_tier1_aligns.begin();
      for( ; tier1_iter != top_tier1_aligns.end(); ++tier1_iter){
        blah(&(tier1_iter->second), ScoreBase::A_tier2_score_class,
             ScoreBase::A_tier2_score_cmp, &top_tier2_aligns,
             dbase_ligand, true, A_max_num_aligns, 
             A_score_cutoff, A_fine_tune_tier2_alignments,
             dset_site);
      }
    }

    // refine best alignment
    if(top_tier2_aligns.size() && A_fine_tune_best_tier2_alignment){
      refine_best_alignment(ScoreBase::A_tier2_score_class, &top_tier2_aligns,
                            dset_site, dbase_ligand, dbase_struct_id,
                            dbase_mol_id, rigid_results_out);
    }

    if(dbase_ligand){
      std::cout << "prot-lig scoring is not implemented\n";

      tier1_score_mmi tier1_iter = top_tier1_aligns.begin();
      for( ; tier1_iter != top_tier1_aligns.end(); ++tier1_iter){
        if(tier1_iter->second.frag_file) delete tier1_iter->second.frag_file;
        tier1_iter->second.frag_file = 0;
      } 

      tier2_score_mmi tier2_iter = top_tier2_aligns.begin();
      for( ; tier2_iter != top_tier2_aligns.end(); ++tier2_iter){
        if(tier2_iter->second.frag_file) delete tier2_iter->second.frag_file;
        tier2_iter->second.frag_file = 0;
      } 
      // compute affi & orient score
//      if(top_tier2_aligns.begin()->second->compute_prot_lig_score())
//        score_protein_ligand_interactions(&top_tier2_aligns);

      // Write out ligands and reclaim memory
//      handle_ligands(top_tier2_aligns, dbase_ligand, dbase_mol_id);
      delete(dbase_ligand);
      dbase_ligand = 0;
    }

    cnt = 1;
    tier2_score_mmi tier2_iter = top_tier2_aligns.begin();
    for( ; tier2_iter != top_tier2_aligns.end(); ++tier2_iter){
      //std::cout << dbase_struct_id << " " << tier2_iter->second.score << std::endl;
      tier2_iter->second.write_score_fields(results_out, cnt, false, "",
                                            dbase_struct_id, dbase_mol_id);
      results_out << "\n";
    
    }

    return true;
  }

  // somehow pointers are getting messed up with the templates & g++
  // the compiler is getting the wrong signals when using T* and cannot
  // find the function for a pointer to score_method.  If I use T, then
  // it assumes that the template arg is an object and not a pointer or
  // iterator
  template <typename score_T, typename cmp_T> bool
  blah(align_T *align, score_T &score_method, cmp_T score_cmp,
       std::multimap<my_float_t, align_T, cmp_T> *top_aligns,
       mol2File *lig_file, const bool frag_lig_filter,
       const size_t max_num_aligns, const my_float_t score_cutoff,
       const bool fine_tune_all, DbaseSitemap *dset_site)
  {
    typedef typename std::multimap<my_float_t, align_T, cmp_T>::iterator
      my_mmi;
//    std::cout << "Score method uses surface_mesh: " 
//              << score_method.uses_surface_mesh() << std::endl;
    A_model->revert(score_method.uses_surface_mesh(),
                    score_method.uses_hbond_surfaces());
    A_model->inverse_transform(align->R, align->T,
                               score_method.uses_surface_mesh(),
                               score_method.uses_hbond_surfaces());

    // we may want multiple sets of match prints in the future?
    align->match_print.resize(A_model->fit_points_size());

    // Check if we want to fine tune all candidate alignments at this level
    align->score = score_method.score(*A_model, *dset_site, align);
    align->score = (align->score - A_mu) / A_sigma;
    if(fine_tune_all){
//      std::cout << "fine tuning not finished\n";
      fine_tune_align(dset_site, align);
      std::cout << "IK not finished\n";
//      if(A_do_IK)
//        A_my_IK.run(model(), *dset_site, &(*align), "TEST_ID", std::cout);
    }

#if 0
    align->write_score_fields(std::cout, 0, false, "",
                              dset_site->ligand_file_name(), "");
    std::cout << "\n";
#endif


    if(score_cmp(score_cutoff, align->score)) return false;

    // If we have the max number of alignments to keep and this score is
    // worse than all of the previously stored alignemnts, discard it
    my_mmi pos = top_aligns->upper_bound(align->score);
    if(top_aligns->size() >= max_num_aligns && pos == top_aligns->end()){
      return false;
    }

    // If we have a ligand file, determine the ligand fragment in the
    // binding site and check if enough ligand atoms are in the fragment
    align_T new_align = *align;
    new_align.frag_file = 0;
    if(lig_file && (frag_lig_filter || fine_tune_all)){
      int lig_has_enough_atoms =
        get_ligand_fragment(*lig_file, A_model->site_volume_estimate_handle(),
                            A_model->interacting_atoms(), &new_align);
      if(frag_lig_filter && !lig_has_enough_atoms){
        if(new_align.frag_file) delete(new_align.frag_file);
        new_align.frag_file = 0;
        return false;
      }
    }

    // Insert the alignment into the top_aligns, and if we have an extra
    // saved alignment, remove the alignment with the poorest score
    top_aligns->insert(pos, align_pair(new_align.score, new_align));
    if(top_aligns->size() > max_num_aligns){
      my_mmi last = top_aligns->end();
      --last;
      if(last->second.frag_file) delete(last->second.frag_file);
      last->second.frag_file = 0;
      top_aligns->erase(last);
      return true;
    }

    return true;
  }

  template <class score_T> void
  refine_best_alignment(score_T &score_method, tier2_score_mmap *top_aligns,
                        DbaseSitemap* dset_site, mol2File *lig_file,
                        const std::string& db_struct_id,
                        const std::string& db_mol_id,
                        std::ostream& rigid_results_out)
  {
    // Get the best scoring alignment -- need a deep copy
    align_T top_align = top_aligns->begin()->second;
    top_align.frag_file = 0;

    // Refine the alignment
    A_model->revert(score_method.uses_surface_mesh(),
                    score_method.uses_hbond_surfaces());
    A_model->inverse_transform(top_align.R, top_align.T,
                               score_method.uses_surface_mesh(),
                               score_method.uses_hbond_surfaces());
    top_align.score = score_method.score(*A_model, *dset_site, &top_align);
    top_align.score = (top_align.score - A_mu) / A_sigma;
    fine_tune_align(dset_site, &top_align);
    std::cout << "Ik is not implemented\n";
    //if(A_do_IK)
      //A_my_IK.run(A_model, *dset_site, &(*top_align), db_struct_id, std::cout);

    // Clear the given map holding the "unrefined" alignments and insert
    // the refined alignment and its score
    tier2_score_mmi aligns_iter = top_aligns->begin();
    for( ; aligns_iter != top_aligns->end(); ++aligns_iter)
      if(aligns_iter->second.frag_file){
        delete(aligns_iter->second.frag_file);
        aligns_iter->second.frag_file = 0;
      }

    // Update ligand fragmentation if we have a ligand
    if(lig_file){
      get_ligand_fragment(*lig_file, A_model->site_volume_estimate_handle(),
                          A_model->interacting_atoms(), &top_align);
    }
    top_aligns->clear();
    top_aligns->insert(top_aligns->begin(),
                       align_pair(top_align.score, top_align));
  }





private:

  std::vector<Quaternion> A_rot_grid;
  ModelSitemap *A_model;
  bool A_fine_tune_best_tier2_alignment;
  bool A_fine_tune_tier2_alignments;
  bool A_save_rigid_results;
  my_float_t A_surf_pt_w;
  my_float_t A_hb_cap_w;
  my_float_t A_mu;
  my_float_t A_sigma;
  size_t A_max_num_aligns;
  size_t A_max_tier1_aligns;
  my_float_t A_score_cutoff;

};


class test_caps_N_surf_ICP : public Search{
public:
  test_caps_N_surf_ICP(const SearchParameters* args_in) : Search(args_in)
  { ; }

  virtual ~test_caps_N_surf_ICP()
  { ; }

  bool
  run()
  {
    if(fail()){
      warn("test_caps_N_surf_ICP", "run", 
           "Cannot run because an error has occured");
      return false;
    }

    if(args()->time_process) start_timer();
    align_method_t align_method;
    if(!set_alignment_method(&align_method)) return false;
    std::vector<align_w_rigid_data_t> alignments;

    my_float_t surf_pt_w = 1.0, hbond_pt_w = 0.0;
    // Given our that our implementation of ICP uses a closed form to solve
    // for the optimum transformation, we require that the max weight for
    // each be <= 1.0 and one should be 1.0
    if(args()->fine_tune_ratio > 0.0){
      if(args()->fine_tune_ratio > 1.0){
        surf_pt_w = 1.0;
        hbond_pt_w = 1.0 / (args()->fine_tune_ratio);
      }else{
        surf_pt_w = args()->fine_tune_ratio;
        hbond_pt_w = 1.0;
      }
    }if(args()->fine_tune_ratio < 0.0){
      surf_pt_w = 0.0;
      hbond_pt_w = 1.0;
    }
    std::cout << "Surf, hbond weights (resp.): " << surf_pt_w << " " 
              << hbond_pt_w << "\n";
  
    //std::cout << "Using hbond caps+ scoring method -- training method\n";

    if(args()->score_str == "ModelSiteRMSD"){
      std::cout << "Only computing the RMSD of the alignments\n"
                << "NOTE: this requires that the sites are already aligned\n\n";
      score(&alignments, align_method,
              new GridSamplingAndScoring< NoTier1Score, ModelSiteRMSD,
                                         align_w_rigid_data_t >
              (model(), *(args()), surf_pt_w, hbond_pt_w));

    }else if(args()->score_str == "WeightedSumsScore"){
      std::vector<rigid_align_t> blah;
      score(&blah, align_method,
              new GridSamplingAndScoring< NoTier1Score, WeightedSumsScore,
                                         rigid_align_t >
              (model(), *(args()), surf_pt_w, hbond_pt_w));
    
    }else{
      score(&alignments, align_method, 
            new GridSamplingAndScoring< WeightedSumsScore, point_and_surf_score, 
                                      align_w_rigid_data_t >
            (model(), *(args()), surf_pt_w, hbond_pt_w));
    }


    // no if statement here -- we always run the test case
    //get_timer_and_write_to_file();  
    return true;
  }

private:
  
};

int main(const int argc, const char **argv)
{
  std::cout << "\n" << argv[0] << " (" << PACKAGE_NAME << ") " 
            << PACKAGE_VERSION << "\n\n";

  std::cout << "size of rigid: " << sizeof(rigid_align_t) << std::endl;
  std::cout << "size of rigid w: " << sizeof(align_w_rigid_data_t) << std::endl;

  // Do not return -1 here since system, fork, etc return -1 on failure, and we
  // wish to distinguish between system and program failure
  SearchParameters my_params(argc, argv);
  BaseParameters::status_t status = my_params.status();
  if(status == BaseParameters::DISPLAY_HELP_ONLY) return 0;
  else if(status != BaseParameters::READY){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize parameters\n";
    return 1;
  }
  
//  my_params.use_hbond_surfaces = true;


  test_caps_N_surf_ICP my_search(&my_params);
  if(my_search.fail()){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize search\n";
    return 1;
  }

  if(!my_search.run()) return 1;
  return 0; 
}
