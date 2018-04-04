/******************************************************************************
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
#ifndef _SITEMAP_ATOMS_FILE_HEADER_INCLUDED
#define _SITEMAP_ATOMS_FILE_HEADER_INCLUDED

#include <PDBStructure.H>
#include <BoundingVolume.H>

namespace ASCbase{

  //! The structure atoms which gave rise to the sitemap (interaction) points
  class SitemapAtomsFile : public PDBBase {
  public:

    //! Load a sitemap atoms file
    /*! Load the sitemap atoms file given the path stored in filename
     *
     * @param filename Path to the sitemap atoms file to load
     * @param verbosity If is REPORT, send something to stdout
     */
    SitemapAtomsFile(const std::string filename, 
                     const verbose_level_t verbosity = VERBOSE_SILENT);
  
    //! Populate the atoms by appending all the atoms in the input table (map)
    /*! 
     * @param a_tbl Reference to a table (map) holding atoms keyed by atom
     *        number (PDB serial)
     */
    SitemapAtomsFile(const atom_map_t& a_tbl);

    //! Populate the atoms by appending all residues which have at least one
    //! atom within interaction distance of any ligand atom
    SitemapAtomsFile(PDBStructure *prot_atoms, const atom_map_t& a_tbl,
                     const std::string lig_fname, const int min_chain_size,
                     const my_float_t max_interact_dist = A_max_interact_dist);

    //! Populate the atoms by appending all residues which have at least one
    //! atom within the "original" bounding volume plus an additional buffer of
    //! max_interaction_dist
    SitemapAtomsFile(PDBStructure *prot_atoms, BoundingVolume* site_vol,
                     const int min_chain_size, const atom_map_t& a_tbl);

    //! Do nothing dstr
    ~SitemapAtomsFile();

    //! Write all the non-hydrogen atoms to 'out' in xyzr format
    /*! Note: if one want certain types of atoms, etc, they should be 
     * specified through the constructor interface or possibly 
     * use PDBStructure instead
     */
    void
    write_xyzr(std::ostream &out, const atom_vci close_atom, 
               uint* close_atom_idx);

  private:

    static const my_float_t A_max_interact_dist; //!< maximum prot-lig atom interaction distance
    static const std::string A_fname; //!< Name of source file
  };
}

#endif
