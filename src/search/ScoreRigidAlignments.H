/*****************************************************************************e
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the ASCbase Software project.
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

#ifndef _SCORE_RIGID_ALIGNMENTS_HEADER_FILE
#define _SCORE_RIGID_ALIGNMENTS_HEADER_FILE

#include <ScoreMapBase.H>
#include <AlignDataWithPreIKData.H>
#include <AlignDataWithProtLigScore.H>
#include <new_tri_test.H>
#include <DbaseSitemap.H>
#include <ModelSitemap.H>
#include <ProtLigScore.H>
#include <IK_tests.H>


namespace ASCbase{

//! Base class for scoring and storing the rigid alignments of dbase 
//! sitemaps to a given model (query).
/*! At the present we are not searching huge databases, but this is the 
 * future.  We may need to rearrange some of the storage to be more space
 * efficient.  The reason for separating the scoring and generation of
 * alignments was to separate the two ideas and remove the code 
 * dependancies between generation and scoring of alignments.  The thought is
 * to keep the generation of alignments relatively fixed and focus on the
 * scoring of alignments.  If the generation of alignments is drastically
 * changed, it is likely the scoring methods will need to be completely
 * rewritten.
 *
 * notes on template parameters:
 * tier1_SF: must be a derived class of SitesScoreBase
 * tier2_SF: must be a derived class of SitesScoreBase
 * align_T: must be rigid_align_t or derived class
 */
template < class tier1_SF, class tier2_SF, class align_T >
class ScoreRigidAlignments : public ScoreMapBase {
public:
  typedef align_T                                       rigid_align_obj;
  typedef std::vector<align_T>                          rigid_align_vec; 
  typedef typename rigid_align_vec::iterator            rigid_align_vi; 
  typedef typename rigid_align_vec::const_iterator      rigid_align_vci; 
  typedef std::pair<my_float_t, rigid_align_vi>         align_pair;

  typedef typename tier1_SF::score_cmp                  tier1_score_cmp;
  typedef typename tier2_SF::score_cmp                  tier2_score_cmp;
  typedef std::multimap<my_float_t, rigid_align_vi, tier1_score_cmp>
                                                        tier1_score_mmap;
  typedef std::multimap<my_float_t, rigid_align_vi, tier2_score_cmp>
                                                        tier2_score_mmap;
  typedef typename tier1_score_mmap::iterator           tier1_score_mmi;
  typedef typename tier2_score_mmap::iterator           tier2_score_mmi;
  typedef typename tier1_score_mmap::const_iterator     tier1_score_mmci;
  typedef typename tier2_score_mmap::const_iterator     tier2_score_mmci;

  //! Default constructor for scoring methods of rigidly aligned sitemaps 
  /*! NOTE: once an alignment vector is passed to ScoreRigidAlignments DO
   * NOT manipulate it in any way since ScoreRigidAlignments uses iterators
   * pointing to the vector holding the alignments
   *
   * @param model_in The model sitemap (binding pocket, etc.)
   * @param params Reference to the search's runtime parameters
   */
  ScoreRigidAlignments(ModelSitemap* model_in, const SearchParameters& params)
    : ScoreMapBase(params), A_my_IK(params)
  {
    A_orient_num = 0;
    A_model = model_in;
    A_score_cutoff = params.score_cutoff;
    A_max_num_aligns = params.num_scores_to_keep;
    std::cout << "number of scores to keep: " << A_max_num_aligns << std::endl;
    A_max_tier1_aligns = params.max_tier1_aligns;
    A_min_num_atoms = params.min_num_atoms;
    A_do_IK = params.do_IK;
    normalize_scores = params.normalize;
    dbase_ligs = params.dbase_ligs;
    A_proj_output = params.proj_output;
    if(params.score_str.length()) A_score_method_str = params.score_str;
    if(normalize_scores){
      const my_float_t* norms = A_model->get_norm_moments();
      A_mu = norms[0];
      A_sigma = norms[1];
    }else{
      A_mu = 0;
      A_sigma = 1;
    }
    A_struct_id_field = params.add_struct_id_field;
  
    std::string chaff;
    get_path_and_struct_id(A_model->atoms_file_name(), &chaff, 
                           &A_model_struct_id);
    A_model_struct_id = 
      A_model_struct_id.substr(0, A_model_struct_id.length() - 2);
    if(A_prot_lig_score){
      std::cout << "Protein-ligand scoring will be computed for saved "
                << "orientations\n";
      A_query_prot = new PDBStructure(A_model->atoms_file_name());
    }else{ 
      A_query_prot = 0;
      std::cout << "Protein-ligand scoring is OFF\n";
    }

    A_write_ligands = params.write_ligands;
    A_fine_tune_tier2_alignments = params.fine_tune_tier2_alignments;
    A_fine_tune_best_tier2_alignment = params.fine_tune_best_tier2_alignment;
//    if(A_fine_tune_tier2_alignments || A_fine_tune_best_tier2_alignment)
//      std::cout << "need to check for convergence\n";

    if(params.scale_terms){
      std::cout << "Using scaled terms for scoring\n";
      A_tier1_score_class.set_scale_terms(true);      
      A_tier2_score_class.set_scale_terms(true);      
    }
  }

  //! basic destruction
  ~ScoreRigidAlignments() 
  { 
    if(A_query_prot) delete A_query_prot; 
    A_query_prot = 0;
  }

  //! Main scoring routine 
  /*! score the alignments for 1 model, search pair of sitemaps
   *
   * NOTE: once an alignment vector is passed to ScoreRigidAlignments DO
   * NOT manipulate it in any way since ScoreRigidAlignments uses iterators
   * pointing to the vector holding the alignments.  The vector can be 
   * safely modified after this class is no longer needed.
   */
  bool
  score_alignments(rigid_align_vec& aligns, DbaseSitemap* dbase_site,
                   std::ostream& results_out = std::cout,
                   std::ostream& rigid_results_out = std::cout)
  {
    mol2File *dbase_ligand = 0;
    std::string dbase_struct_id, dbase_mol_id;
    if(!load_ligand(dbase_site->ligand_file_name(), 
                    dbase_site->atoms_file_name(), 
                    &dbase_ligand, &dbase_struct_id, &dbase_mol_id)) 
      return false;
    else if(A_prot_lig_score)
      dbase_ligand->calc_charge_sums();

    // forced flush for debugging
    std::cout << "Comparing " << A_model_struct_id << " to " 
              << dbase_struct_id << std::endl;

    tier2_score_mmap top_tier2_aligns;
    if(A_tier1_score_class.score_is_noop()){
      // "Single tier" scoring
//      std::cout << "inside \"single tier\" scoring --- num aligns: " 
//                << aligns.size() << "\n\n";
      if(!score_alignments(aligns, dbase_site, dbase_ligand, dbase_struct_id, 
                           dbase_mol_id, &top_tier2_aligns, rigid_results_out)){
        if(dbase_ligand) delete(dbase_ligand); 
        dbase_ligand = 0; 
        return false; 
      } 
    }else{
      // First tier of scoring (scoring sieve)
      tier1_score_mmap top_tier1_aligns;
      if(!score_tier1_alignments(aligns, dbase_site, dbase_ligand, 
                                 &top_tier1_aligns)){
        if(dbase_ligand) delete(dbase_ligand);
        dbase_ligand = 0;
        return false;
      }

      // Second tier scoring (second scoring sieve)
      if(!score_tier2_alignments(top_tier1_aligns, dbase_site, dbase_ligand, 
                                 dbase_struct_id, dbase_mol_id,
                                 &top_tier2_aligns, rigid_results_out)){
        if(dbase_ligand) delete(dbase_ligand);
        dbase_ligand = 0;
        return false;
      }
    }

    if(dbase_ligand && top_tier2_aligns.size() > 0){
      // compute affi & orient score
      if(top_tier2_aligns.begin()->second->compute_prot_lig_score())
        score_protein_ligand_interactions(&top_tier2_aligns);

      // Write out ligands and reclaim memory
      handle_ligands(top_tier2_aligns, dbase_ligand, dbase_mol_id);
      delete(dbase_ligand);
      dbase_ligand = 0;
    }

    report_alignments(top_tier2_aligns, A_write_ligands, dbase_struct_id, 
                      dbase_mol_id, results_out);
    return true;
  }

  //! Needs to be called before scoring any alignments
  virtual void 
  write_score_header(std::ostream& out)
  {
   if(normalize_scores){
      out << std::setw(50) << "# Mean of raw scores of diverse sitemaps:"
          << A_mu << "\n"
          << std::setw(50) << "# Stdev of raw scores of diverse sitemaps:"
          << A_sigma << "\n";
    }else out << "# ASCbase alignments scores are not normalized\n";
  
    out << "#\n# Fields:\n";
    if(A_struct_id_field) out << "# 0) Database structure id\n";

    std::vector<std::string> fields;

    std::string tmp = 
      "Name of ligand fragment corresponding to the score record (line)";
    fields.push_back(tmp);

    align_T tmp_align; 
    tmp_align.get_score_field_labels(&fields, normalize_scores);
   
    if(ext_SF_id.length() && ext_SF_id != "NONE"){
      std::vector<std::string> sf_names;
      get_ext_SF_names(&sf_names);
      for(uint i = 0; i < sf_names.size(); ++i){
        tmp = sf_names[i];
        tmp += " score of target ligand fragment in query pocket";
        fields.push_back(tmp);
      }    
    }

    if(A_fine_tune_tier2_alignments){
      out << "# NOTE: the matchprint is NOT valid in the RIGID file\n";
    }

    for(uint i = 0; i < fields.size(); ++i){
      out << "# " << std::setw(2) << i+1 << ") " << fields[i] << "\n";
    }

  }

  //! Score the alignments for 1 model, search pair of sitemaps
  /*! NOTE: once an alignment vector is passed to ScoreRigidAlignments DO
   * NOT manipulate it in any way since ScoreRigidAlignments uses iterators
   * pointing to the vector holding the alignments.  The vector can be 
   * safely modified after this class is no longer needed.
   */
  bool 
  score_tier1_alignments(rigid_align_vec& aligns_in, DbaseSitemap* search,
                         mol2File *lig_file, tier1_score_mmap *top_aligns)
  {
    if(aligns_in.empty()) return false;
    if(A_min_num_atoms > 0 && !lig_file) return false;

    // we use an exponent of 30 as we may have float or double and
    // I don't want to worry about overflow because of nonsymmetry of 
    // machine numbers
    my_float_t score_cutoff = (A_tier1_score_cmp(0,1) ? 1E+30 : -1E+30);

    rigid_align_vi align;
    for(align = aligns_in.begin(); align < aligns_in.end(); ++align){
      align->set_number_of_orientations(aligns_in.size());
      blah(align, &A_tier1_score_class, A_tier1_score_cmp, top_aligns, lig_file,
           true, A_max_tier1_aligns, score_cutoff, false, search);
    }
    return true;
  }

  //! Score the alignments that passed the first sieve
  bool 
  score_tier2_alignments(tier1_score_mmap &align_map, DbaseSitemap* dset_site,
                         mol2File *lig_file, const std::string& db_struct_id,
                         const std::string& db_mol_id, 
			 tier2_score_mmap *top_aligns,
                         std::ostream& rigid_results_out)
  {
//    std::cout << "need to adjust ligand orientation after rigid refinement of sites\n";
    if(align_map.empty()) return false;
    if(A_min_num_atoms > 0 && !lig_file) return false;

    // Rename align to something else and get a reference to the align_vi
    tier1_score_mmci align_map_iter = align_map.begin();
    for( ; align_map_iter != align_map.end(); ++align_map_iter){
      blah(align_map_iter->second, &A_tier2_score_class, A_tier2_score_cmp,
           top_aligns, lig_file, false, A_max_num_aligns, A_score_cutoff,
	   A_fine_tune_tier2_alignments, dset_site);
    }

    // If we wanted only to fine tune the best scoring alignment at this level,
    // do it here
    if(top_aligns->size() && A_fine_tune_best_tier2_alignment)
      refine_best_alignment(&A_tier2_score_class, top_aligns, dset_site, 
                            lig_file, db_struct_id, db_mol_id, 
			    rigid_results_out);

    return true;
  }

  //! Score the alignments for 1 model, dset pair of sitemaps
  /*! NOTE: once an alignment vector is passed to ScoreRigidAlignments DO
   * NOT manipulate it in any way since ScoreRigidAlignments uses iterators
   * pointing to the vector holding the alignments.  The vector can be 
   * safely modified after this class is no longer needed.
   *
   * This function is for "single tiered" scoring
   */
  bool 
  score_alignments(rigid_align_vec& aligns_in, DbaseSitemap* dset_site,
                   mol2File *lig_file, const std::string& db_struct_id,
                   const std::string& db_mol_id, tier2_score_mmap *top_aligns, 
                   std::ostream& rigid_results_out)
  {
    if(aligns_in.empty()) return false;
    if(A_min_num_atoms > 0 && !lig_file) return false;

    rigid_align_vi align;
    for(align = aligns_in.begin(); align < aligns_in.end(); ++align){
      align->set_number_of_orientations(aligns_in.size());
      blah(align, &A_tier2_score_class, A_tier2_score_cmp, top_aligns, 
           lig_file, true, A_max_num_aligns, A_score_cutoff, 
           A_fine_tune_tier2_alignments, dset_site);
    }

    // If we wanted only to fine tune the best scoring alignment, do it here
    if(top_aligns->size() && A_fine_tune_best_tier2_alignment)
      refine_best_alignment(&A_tier2_score_class, top_aligns, dset_site, 
                            lig_file, db_struct_id, db_mol_id, 
			    rigid_results_out);
    return true;
  }

  template <class score_T, class cmp_T> bool
  blah(rigid_align_vi align, score_T *score_method, cmp_T score_cmp,
       std::multimap<my_float_t, rigid_align_vi, cmp_T> *top_aligns,
       mol2File *lig_file, const bool frag_lig_filter,
       const size_t max_num_aligns, const my_float_t score_cutoff,
       const bool fine_tune_all, DbaseSitemap *dset_site)
  {
    typedef typename std::multimap<my_float_t, rigid_align_vi, cmp_T>::iterator
      my_mmi;

    A_model->revert(score_method->uses_surface_mesh(),
                    score_method->uses_hbond_surfaces());
    A_model->inverse_transform(align->R, align->T,
                               score_method->uses_surface_mesh(),
                               score_method->uses_hbond_surfaces());

    // we may want multiple sets of match prints in the future?
    align->match_print.resize(A_model->fit_points_size());

    // Check if we want to fine tune all candidate alignments at this level
    align->score = score_method->score(*A_model, *dset_site, &(*align));
    align->score = (align->score - A_mu) / A_sigma;
    if(fine_tune_all){ 
      fine_tune_align(dset_site, &(*align));
      if(A_do_IK)
        A_my_IK.run(A_model, *dset_site, &(*align), "TEST_ID", std::cout);
    }

    if(score_cmp(score_cutoff, align->score)){
//      std::cout << "failed score cut off of " << score_cutoff 
//                << " with a score of " << align->score << "\n";
      // We may need to clean up the memory for the frag file
      if(align->frag_file) delete(align->frag_file);
      align->frag_file = 0;
      return false;
    }

    // If we have the max number of alignments to keep and this score is
    // worse than all of the previously stored alignemnts, discard it
    my_mmi pos = top_aligns->upper_bound(align->score);
    if(top_aligns->size() >= max_num_aligns && pos == top_aligns->end())
      return false;

    // If we have a ligand file, determine the ligand fragment in the
    // binding site and check if enough ligand atoms are in the fragment
    if(lig_file && (frag_lig_filter || fine_tune_all)){
      bool lig_has_enough_atoms =
        get_ligand_fragment(*lig_file, A_model->site_volume_estimate_handle(),
                            A_model->interacting_atoms(), &(*align));
      if(frag_lig_filter && !lig_has_enough_atoms) return false;
    }

    // Insert the alignment into the top_aligns, and if we have an extra
    // saved alignment, remove the alignment with the poorest score
    top_aligns->insert(pos, align_pair(align->score, align));
    if(top_aligns->size() > max_num_aligns){
      my_mmi last = top_aligns->end();
      --last;
      if(last->second->frag_file) delete(last->second->frag_file);
      last->second->frag_file = 0;
      top_aligns->erase(last);
      return true;
    }

    return false;
  }

  //! Refine the best alignment in the map and replace the map with a 
  //! new map containing only the refined alignment
  /*! We have a fine mess here -- The goal is to be consistent
   * Issue: When we fine tune after insertion into the map, we cannot
   * update the key.  
   * Solution: In this case, we assume only 1 score per site pair will be
   * printed.  Dump the entire map and rebuild it with the refined
   * alignment.  Deletion and adding of the refined alignment may not
   * be desirable as another alignment could prove to have a more
   * favorable score if the refinement fails.
   */
  template <class score_T> void
  refine_best_alignment(score_T *score_method, tier2_score_mmap *top_aligns, 
                        DbaseSitemap* dset_site, mol2File *lig_file, 
			const std::string& db_struct_id, 
                        const std::string& db_mol_id, 
                        std::ostream& rigid_results_out)
  {
    // Get the best scoring alignment
    rigid_align_vi top_align = top_aligns->begin()->second;

    // Refine the alignment
    A_model->revert(score_method->uses_surface_mesh(),
                    score_method->uses_hbond_surfaces());
    A_model->inverse_transform(top_align->R, top_align->T,
                               score_method->uses_surface_mesh(),
                               score_method->uses_hbond_surfaces());
    top_align->score = score_method->score(*A_model, *dset_site, &(*top_align));
    top_align->score = (top_align->score - A_mu) / A_sigma;
    fine_tune_align(dset_site, &(*top_align));
    if(A_do_IK) 
      A_my_IK.run(A_model, *dset_site, &(*top_align), db_struct_id, std::cout);
    
    // Clear the given map holding the "unrefined" alignments and insert
    // the refined alignment and its score
    tier2_score_mmi aligns_iter = top_aligns->begin();
    for( ; aligns_iter != top_aligns->end(); ++aligns_iter)
      if(aligns_iter->second->frag_file){
        delete(aligns_iter->second->frag_file);
        aligns_iter->second->frag_file = 0; 
      }
    top_aligns->clear();
    top_aligns->insert(top_aligns->begin(), 
                       align_pair(top_align->score, top_align));

    // Update ligand fragmentation if we have a ligand
    if(lig_file){
      get_ligand_fragment(*lig_file, A_model->site_volume_estimate_handle(),
                          A_model->interacting_atoms(), &(*top_align));
    }
  }

  //! Score protein-ligand interactions for the saved alignments
  void
  score_protein_ligand_interactions(tier2_score_mmap *top_aligns)
  {
    // Protein-ligand scoring 
    tier2_score_mmi align;
    for(align = top_aligns->begin(); align != top_aligns->end(); ++align){
      if(!align->second->frag_file) continue;
       
      // Note: we are potentially scoring lots of fragments and other junk.
      // A decent guess is to keep the same charges as are on the full ligand
      // since the groups are likely to be similarly charged in a 
      // complementary ligand -- i.e. don't compute ligand charge_sums
      ProtLigScore frag_score(*A_query_prot, *(align->second->frag_file),
                              true, false);
      align->second->set_prot_lig_scores(frag_score.affiscore(),
                                         frag_score.orientscore(),
                                         frag_score.ligand_efficiency());
    }
  }

  //! Handle the ligands for the saved alignments
  void
  handle_ligands(tier2_score_mmap &top_aligns, mol2File *lig_file,
                 const std::string &mol_id)
  {
    // write ligands to file
    if(A_write_ligands){
      tier2_score_mmi align = top_aligns.begin();
      for(uint cnt = 1; align != top_aligns.end(); ++align, ++cnt)
        write_ligands(align->second->frag_file, lig_file, align->second->R,
                      align->second->T, cnt, mol_id, A_proj_output);
    }

    // Delete the memory corresponding to the ligand
    tier2_score_mmi align = top_aligns.begin();
    for( ; align != top_aligns.end(); ++align){
      delete align->second->frag_file;
      align->second->frag_file = 0; 
    }
    
    run_external_scoring_method(top_aligns.begin(), top_aligns.end(),
                                A_model->atoms_file_name(), mol_id,
                                A_proj_output);
  }

  //! Loop to print out the saved alignments
  inline void
  report_alignments(const tier2_score_mmap top_aligns, const bool wrote_ligs,
                    const std::string &db_struct_id, 
                    const std::string &db_mol_id,
                    std::ostream& results_out)
  {

    uint cnt = 1;
    tier2_score_mmci align;
    for(align = top_aligns.begin(); align != top_aligns.end(); ++align, ++cnt){
      align->second->write_score_fields(results_out, cnt, wrote_ligs, ext_SF_id,
                                        db_struct_id, db_mol_id);
      results_out << "\n";
    }
  }

  //! Fine tune the transformation -- i.e. apply ICP to the given 
  //! correspondences
  /*! NOTE: the assumption is the two sites are aligned with respect to the
   * current candidate alignment
   */
  virtual void 
  fine_tune_align(DbaseSitemap* search, rigid_align_obj *align)
  {
    std::cout << "\n\n Inside fine tune align\n";
    const my_float_t RMSE_TOL = 0.05;
    const size_t N = 100;

    // Get current correspondences
    my_float_t* query_pts = 0, *db_pts = 0;
    size_t npts;
    A_tier2_score_class.correspondences(*A_model, &query_pts, &db_pts, &npts);
    std::cout << " Number of initial pt correspondences " << npts 
              << "\n" << std::endl;

    // If the the alignment data class supports it, save the initial 
    // transformation, score, and features
    //align->site_atomic_rmsd() ....
    if(A_do_IK)
      A_my_IK.compute_initial_atomic_rmsd(A_model, *search, &(*align));
    align->save_rigid_alignment_vals();

    // We require at least 3 correspondences to fully describe a 3D alignment
    if(npts < 3) return;

    // Save the initial alignment in the event we need to revert to it
    rigid_align_obj initial_align(*align); 

    
#if 0
      std::stringstream ofname_f;
      ofname_f << "blah_" << A_orient_num << "_0";
      A_model->binding_site_mesh_handle().write(ofname_f.str());
      std::cout << "Corresponds: 0\n";
        for(int zz = 0; zz < npts; ++zz){
	  std::cout << query_pts[3*zz + 0] << " "
	            << query_pts[3*zz + 1] << " " << query_pts[3*zz + 2] << " ";
	  std::cout << db_pts[3*zz + 0] << " "
	            << db_pts[3*zz + 1] << " " << db_pts[3*zz + 2] << "\n";
        }
#endif

    //model->binding_site_mesh_handle()->write("blah_0");
    my_float_t prev_RMSE = my_float_max;
    for(size_t n = 0; n < N; ++n){

      // Get most recent transformation 
      Quaternion Q;
      my_float_t T[3];
      lse_3D_fit(query_pts, db_pts, npts, &Q, T);

      // Clean up mem
      if(query_pts) delete [] query_pts;
      if(db_pts) delete [] db_pts;
      query_pts = 0;
      db_pts = 0;

      // Apply the most recent transformation
      my_float_t R[9];
      Q.get_ortho_rot_mat(R);

      // This is a very inelegant solution, but we need to concentrate our
      // resources elsewhere at this point (rather than on numerics).
      // The main issue we have been facing is with a small number of points
      // 3  (the minimum required).  When we have a small number of points, 
      // sometimes too many of the corresponding points are the same;  in such 
      // a cases, the solution is degenerate and the roots of the characteristic
      // equation contain a complex component.
      if(std::isnan(T[0]) || std::isnan(T[1]) || std::isnan(T[2])){
         std::cout << "Fine tuning failed.\n"; 
         std::cout << "Number of pairs of corresponding points: " 
                   << npts << "\n\n";
         *align = initial_align;
         return;
      } 

//      std::cout << "round #: " << n << "\n";
//      std::cout << "Q: " << Q.str() << "\n"
//                << "T: " << T[0] << " " << T[1] << " " << T[2] << "\n"; 
      A_model->inverse_transform(R, T, A_tier2_score_class.uses_surface_mesh(),
                                 A_tier2_score_class.uses_hbond_surfaces());

      // Get current correspondences
      align->score = A_tier2_score_class.score(*A_model, *search, &(*align));
      my_float_t RMSE = 
        A_tier2_score_class.correspondences(*A_model, &query_pts, &db_pts, &npts);


//        std::cout << "RMSE: " << RMSE << "\n" << std::endl;
      if(RMSE <= RMSE_TOL) break;    
      if(prev_RMSE - RMSE < 1E-07) break;
      prev_RMSE = RMSE;

#if 0
        std::stringstream ofname;
        ofname << "blah_" << A_orient_num << "_" << n+1;
        A_model->binding_site_mesh_handle().write(ofname.str());
 
      std::cout << "\nCorresponds: " << n + 1 << "\n";
        for(int zz = 0; zz < npts; ++zz){
	  std::cout << query_pts[3*zz + 0] << " "
	            << query_pts[3*zz + 1] << " " << query_pts[3*zz + 2] << " ";
	  std::cout << db_pts[3*zz + 0] << " "
	            << db_pts[3*zz + 1] << " " << db_pts[3*zz + 2] << "\n";
        }

#endif
    }

    // Get a single transformation that moves the dbase site to the query
    A_model->get_current_inverse_3D_transform(&(align->Q), align->T);
    align->Q.get_ortho_rot_mat(align->R);
    align->score = (align->score - A_mu) / A_sigma;

    // Clean up mem
    if(query_pts) delete [] query_pts;
    if(db_pts) delete [] db_pts;
    query_pts = 0;
    db_pts = 0;
  }

  inline const std::string&
  model_struct_id()
  { return A_model_struct_id; }

    tier1_SF A_tier1_score_class;
    tier2_SF A_tier2_score_class;

  protected:
    
    // protected for now -- testing
    tier1_score_cmp A_tier1_score_cmp;  //!< Get an instance of tier1 cmp
    tier2_score_cmp A_tier2_score_cmp;  //!< Get an instance of tier2 cmp

    //! Load the database ligand file
    bool
    load_ligand(const std::string &lig_fname, const std::string &atoms_fname,
                mol2File **ligand, std::string *struct_id, std::string *mol_id)
    {
      *ligand = 0;
      std::string chaff;
      get_path_and_struct_id(atoms_fname, &chaff, struct_id);
      // Using a brain dead method for "get_path_and_struct_id" and
      // passing in XXXXXX_rad.pdb
      *struct_id = struct_id->substr(0, struct_id->length() - 2);

      // If the minimum fragment size is greater than zero, we require each
      // dbase entry to have a valid ligand file
      if(A_min_num_atoms > 0){
        if(lig_fname == "NONE"){
          message(std::string("Skipping the site ") + *struct_id + 
                  " because it does not have a ligand.\n");
          return false;
        }else if(lig_fname == "NOT_FOUND"){
          message(std::string("Skipping the site ") + *struct_id + 
                  " because the ligand file could not be found");
          return false;
        }else if(!get_path_and_struct_id(lig_fname, &chaff, mol_id)){
          message(std::string("Could not get ligand name and/or directory\n"));
          return false;
        }
        *ligand = new mol2File(lig_fname);
        if(!*ligand || (*ligand)->fail()){
          if(*ligand) delete *ligand;
          *ligand = 0;
          message(std::string("Skipping the site map ") + *struct_id + 
                  " because we could not open the ligand file.\n");
          return false;
        }
      }else{
        bool have_lig = true;
        if(lig_fname == "NOT_FOUND"){
          message(std::string("Ignoring the ligand for ") + *struct_id +
                  " because the ligand file could not be found");
          have_lig = false;
        }else if(lig_fname == "NONE"){
          message(std::string("Ignoring the ligand for ") + *struct_id + 
                  " because it does not have a ligand.\n");
          have_lig = false;
        }else if(!get_path_and_struct_id(lig_fname, &chaff, mol_id)){
          message(std::string("Ignoring the ligand for ") + *struct_id +
                  " because the ligand path could not be parsed");
          have_lig = false;
        }

        if(have_lig){
          *ligand = new mol2File(lig_fname);
          if(!*ligand || (*ligand)->fail()){
            if(*ligand) delete *ligand;
            *ligand = 0;
            message(std::string("Ignoring the ligand for ") + *struct_id + 
                    " because we could not open the ligand file.\n");
          }
        }
      }

      //if(A_prot_lig_score && *ligand) (*ligand)->calc_charge_sums();
      return true;
    }

    bool
    get_ligand_fragment(const mol2File& lig_file, 
                        const BoundingVolume &site_vol, 
                        const CoordFile& interacting_atoms, 
                        rigid_align_t *align)
    {
      mol2File* tmp_frag = new mol2File(lig_file, site_vol, interacting_atoms,
                                        A_min_num_atoms, align->R, align->T);

      // Check to be sure the ligand fragment has a reasonable number of 
      // atoms -- NOTE in the case of empty pocket searches we are 
      // assuming A_min_num_atoms is 0 (currently the definition of an 
      // empty pocket search)
      if(tmp_frag->num_heavy_atoms() < A_min_num_atoms){
        delete tmp_frag;
        tmp_frag = 0;
      }else{
        align->frag_file = tmp_frag;
        // Set the binary string for the ligand fragment atom flags
        align->frag_atoms_flags.resize(lig_file.num_atoms());
        std::vector<bool>::const_iterator f_src =
           tmp_frag->frag_atoms_flags_beg();
        std::vector<bool>::iterator f_dest = 
          align->frag_atoms_flags.begin();
        for( ; f_src < tmp_frag->frag_atoms_flags_end(); ++f_src, ++f_dest)
          *f_dest = *f_src;
        return true;
      }
      return false;
    }

    ModelSitemap* 
    model()
    { return A_model; }

    const my_float_t
    mu() const
    { return A_mu; }

    const my_float_t
    sigma() const
    { return A_sigma; }

    const bool
    fine_tune_tier2_alignments() const
    { return A_fine_tune_tier2_alignments; }

    IK_tests&
    IK_tests_handle()
    { return A_my_IK; }   

    const bool
    do_IK() const
    { return A_do_IK; }

  private:
    //tier1_SF A_tier1_score_class;
    ModelSitemap* A_model; //!< The model sitemap -- this is not a constant pointer since we may desire to transform the model rather than dbase site.
    std::string A_score_method_str;
    IK_tests A_my_IK;
    bool A_struct_id_field;
    std::string A_model_struct_id;
    my_float_t A_mu;      //!< Mean score of the model sitemap 
    my_float_t A_sigma;   //!< Stdev of the model sitemap scores 
    my_float_t A_score_cutoff;  //!< Poorest acceptible score
    size_t A_max_num_aligns;      //!< Number of top scoring hits to keep
    size_t A_max_tier1_aligns;    //!< Max # of alignments to pass to tier2 -- sorted by tier1 cmp
    uint A_min_num_atoms;   //!< num lig atoms required in model pocket
    bool normalize_scores; //!< True => normalize scores, otherwise report raw scores
    std::string dbase_ligs; //!< Path to the database ligands
    std::string A_proj_output;
    PDBStructure *A_query_prot; //!< Pointer to the query protein
    bool A_write_ligands;
    bool A_fine_tune_tier2_alignments;
    bool A_fine_tune_best_tier2_alignment;
    bool A_do_IK;

    // for testing and messing around
    uint A_orient_num;

    //! Write out the saved ligands and ligand fragments
    /*!
     * @param output_dir Directory to place the ligands and ligand fragments
     *        directories
     */
    void 
    write_ligands(mol2File* frag_file, mol2File* lig_file, const my_float_t* R,
                  const my_float_t* T, const uint cnt,
                  const std::string struct_id, const std::string output_dir)
    {
      lig_file->revert();
      lig_file->transform(R, T);
      std::ostringstream fname;
      fname << output_dir << "/moved_ligands/" <<  struct_id << "_"
            << std::setfill('0') << std::setw(5) << cnt << "_l.mol2";
      lig_file->write(fname.str());
      if(frag_file){
        fname.str("");
        fname << output_dir + "/ligand_fragments/" + struct_id << "_"
              << std::setfill('0') << std::setw(5) << cnt << "_f.mol2";
        frag_file->write(fname.str());
      }
    }
 

  };
}
#endif
