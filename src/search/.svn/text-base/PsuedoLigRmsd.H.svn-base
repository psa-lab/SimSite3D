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
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

/*
 * $Source: /psa/share/repository/pfizer_proj/src/search/PsuedoLigRmsd.H,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2008-01-04 15:42:48 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2007/11/01 16:35:44  vanvoor4
 * Now pass in the parameters obj instead of those needed
 *
 * Revision 1.1  2007/10/11 16:21:18  vanvoor4
 * Initial checkin
 *
 * 
 */

#ifndef PSUEDO_LIG_RMSD_HEADER_FILE
#define PSUEDO_LIG_RMSD_HEADER_FILE

#include <ScoreRigidAlignments.H>

namespace SimSite3D{

  //! Do a sort of psuedo lig rmsd -- only useful for testing and not intended
  //! to have a use in the final application
  /*!
   *  NOTE:  this method only works for ligands which are EXACTLY the same.
   *  That is the corresponding atom for the nth query ligand atom is the
   *  nth target ligand atom.
   */
  class PsuedoLigRmsd : public ScoreRigidAlignments<std::less<my_float_t> >{
  public:
    //! Loads the model ligand and calls base cstr
    /*
     * @param model_in Pointer to the model/query sitemap
     * @param params Reference to the search's runtime parameters
     */
    PsuedoLigRmsd(ModelSitemap* model_in, const SearchParameters& params);

    //! basic destruction
    ~PsuedoLigRmsd(); 

  protected:
    //! Given an alignment of search to the query, score said alignment
    /*!
     * @param search Pointer to the sitemap aligned to the model sitemap
     * @return The score of the alignment
     */
    my_float_t score(const DbaseSitemap& search, rigid_align_vi scores);

    my_float_t
    correspondences(my_float_t **query_pts_ptr, my_float_t **db_pts_ptr,
                    size_t *npts)
    {
      *query_pts_ptr = 0;
      *db_pts_ptr = 0;
      *npts = 0;
      return my_float_max;
    }

  private:
    static const std::string _fname; //!< "PsuedoLigRmsd.C"
    mol2File* model_lig;
  };
}
#endif
