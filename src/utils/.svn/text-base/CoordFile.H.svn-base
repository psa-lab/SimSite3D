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

#ifndef COORD_FILE_HEADER_INCLUDED
#define COORD_FILE_HEADER_INCLUDED

#include <basics.H>
#include <atom.H>
#include <PDB_metals.H>

namespace SimSite3D{
  
  class CoordFile;

  inline my_float_t simple_rmsd(const CoordFile& f1, const CoordFile& f2);

  //! Base class to define a coordinate file (PDB, mol2, etc).  The only concept
  //! here is atoms and (for now) connectivity should be handled in the derived
  //! classes.
  /*! NOTE: Once a molecule is read NEVER modify the underlying vectors or 
   * arrays or you will endure the immense pain of invalidated iterators.  At
   * the present this is entirely reasonable because we are not in the 
   * business of growing or shrinking molecules.
   */
  class CoordFile{
  public:
    friend my_float_t simple_rmsd(const CoordFile& f1, const CoordFile& f2);

    //! To be used when a file is to be created; initializes A_fail to false.
    CoordFile();

    //! To be used when a file exists
    /*!
     * @param filename_in Path to the file to open
     * @param verbosity verbose level
     */
    CoordFile(const std::string filename_in, 
              const verbose_level_t verbosity = VERBOSE_SILENT);

    //! Only copy the name of the file since some higher level copy constructors
    //! are used to perform operations on the returned copies.
    CoordFile(const CoordFile& src);

    //! Do nothing virtual dstr
    virtual ~CoordFile();

    //! Each derived class should support some routine to report on the file
    /*!
     * @param out Reference to the ostream to write the report
     */
    virtual void 
    report_stats(std::ostream& out) = 0;
  
    //! Each derived class should support some routine to write a molecule
    /*!
     * @param filename Name of the file to write molecule
     */
    virtual bool 
    write(const std::string filename) = 0;
 
    // Someday may want to look at writing iterators the "correct" way
    // Need to setup the arrays
    atom_vci 
    atoms_begin() const 
    { return atoms.begin(); }

    atom_vci 
    atoms_end() const
    { return atoms.end(); }

    const my_float_t*
    positions_begin() const
    { return atoms.points_begin(); }

    const my_float_t*
    positions_end() const
    { return atoms.points_end(); }

    //! Assumes each atom in a coordinate file has a unique key (atom number)
    atom_vci 
    get_atom(const uint atom_num)
    { 
      if(atoms.size() > atoms_by_num.size()) update_atoms_by_num();
      atom_map_t::const_iterator iter = atoms_by_num.find(atom_num);
      if(iter != atoms_by_num.end()) return iter->second;
      else return atom_t::NULL_ATOM_VCI;
    }


    //! Check if any errors have occured when reading the atoms from file
    bool 
    fail() const 
    { return A_fail; }

    //! Get a string holding the name of the coordinate file
    std::string 
    name() const 
    { return name_of_file; }

    //! Transform the positions
    void
    transform(const my_float_t* R, const my_float_t* T)
    { atoms.transform(R,T); }

    //! Transform the positions by the inverse transformation 
    void
    inverse_transform(const my_float_t* R, const my_float_t* T)
    { atoms.inverse_transform(R,T); }

    //! Revert the current positions to the original positions
    void
    revert() 
    { atoms.revert(); }

    void
    bin_coordinates(const my_float_t max_dist)
    { atoms.bin_points(max_dist); }

    //! Get the closest atom center to the given position
    atom_vci
    closest_atom(const my_float_t* pos, my_float_t* d) const
    { return atoms.closest_point(pos, d); }

    //! Get all atoms within radius of the given positition
    void
    close_atoms(const my_float_t *pos, const my_float_t radius,
                std::multimap<my_float_t, atom_vci>* pts_map) const
                //float_const_iter_map* pts_map) const
    { atoms.close_points(pos, radius, pts_map); }

    //! Determine if any protein atoms significantly overlap/clash with the
    //! atom a
    /*!
     * @param a Reference to the atom that is near the protein 
     * @return True if there is at least one significant overlap, else false.
     */
    const bool
    overlapping_atoms(const atom_t& a) const
    {
      for(atom_vci b = atoms.begin(); b != atoms.end(); ++b)
        if(atoms_overlap(a, *b)) return true;
      return false;
    }
  
    //! Find all protein atoms that overlap atom a
    /*! The purpose of this method is mainly a development/testing method
     * that allows one to see all atoms that have significant overlap with
     * atom a.
     *
     * @param a Reference to the atom that is near the protein 
     * @param overlaps Pointer to a vector to hold the iterators to the atoms
     *        that significantly overlap atom a.
     * @return True if there is at least one significant overlap, else false.
     */
    bool
    overlapping_atoms(const atom_t& a, std::vector<atom_vci>* overlaps)
    {
      overlaps->clear();
      for(atom_vci b = atoms.begin(); b != atoms.end(); ++b)
        if(atoms_overlap(a, *b)) overlaps->push_back(b);
      return (overlaps->size() ? true : false);
    }

    //! Get the number of atoms in the molecule
    const uint
    num_atoms() const
    { return atoms.size(); }

    //! Compute the centroid of the molecule
    const bool
    centroid_3D(my_float_t* C) const
    { return atoms.centroid_3D(C); }

    inline void
    get_atom_bin(atom_vci a, point_storage<atom_t>::bins_vci_type *nbrs, 
                 point_storage<atom_t>::bins_vci_type *moved_nbrs) const
    { return atoms.get_bin(a->pos, nbrs, moved_nbrs); }

    void
    write_xyzr_line(std::ostream& out, const my_float_t* pos,
                    const my_float_t vdw_radius);

    inline void
    write_bin_bounds(std::ostream &out) const
    { atoms.write_bin_bounds(out); }

  protected:
    //! Require all derived classes to implement some sort of read function
    virtual bool read_data(std::ifstream &fin) = 0;

    //! Add atom to the storage
    atom_vci 
    append_atom(const atom_t& a)
    { 
      atoms.push_back(a);
      return atoms.end();
    }

    //! Allow derived classes to modify the atoms via iterators -- beginning 
    //! of the point_storage class
    /*! The name of this function is purposely changed because it using the
     * name atoms_begin() causes g++ to not know which function to choose
     */
    atom_vi 
    atoms_storage_begin()
    { return atoms.begin(); }

    //! Allow derived classes to modify the atoms via iterators -- end 
    //! of the point_storage class
    atom_vi 
    atoms_storage_end()
    { return atoms.end(); }

    //! Object is good to go
    void 
    green_light() 
    { A_fail = false; }
  
    //! An error has occured and object is not safe until error is resolved
    void 
    red_light() 
    { A_fail = true; }

  private:

    void update_atoms_by_num()
    {
      atom_vci a;
      for(a = atoms.begin() + atoms_by_num.size(); a < atoms.end(); ++a)
        atoms_by_num[a->atom_num] = a; 
    }

    //! Check if the atoms have significant and undesireable overlap
    /*! Special rules exist to allow atoms that are participating in hydrogen
     * bonds to be closer than would otherwise be tolerated.  A similar rule
     * exists for acceptor and metal interactions.
     *
     * @param a Reference to the ligand or probe atom
     * @param b Reference to the protein atom (or metal)
     * @return True if there is significant and undesireable overlap between a
     *         and b
     */ 
    const bool 
    atoms_overlap(const atom_t& a, const atom_t& b) const 
    {
      // The assumption is waters and metals shouldn't be in the protein unless
      // they are to be considered as part of the binding site
      my_float_t dd = dist_squared(a.pos, b.pos);

      // Complementary Hbonding atoms
      if((a.act_type == DONOR || a.act_type == DONEPTOR) && 
         (b.act_type == ACCEPTOR || b.act_type == DONEPTOR)){
        if(dd < min_hbond_dist*min_hbond_dist) return true;
      }else if((a.act_type == ACCEPTOR || a.act_type == DONEPTOR) && 
         (b.act_type == DONOR || b.act_type == DONEPTOR)){
        if(dd < min_hbond_dist*min_hbond_dist) return true;
      // Metals
      // Only prot metals for now -- i.e. metals exist in the protein PDB file
      // This was modified as per Matt & Leslie's suggestion -- we could get too
      // close, but in general we are clipping ligand atoms that are rigidly 
      // bonded to the acceptors interacting with the metals
      }else if(b.res == PDB_METAL){    
        const pdb_metal_info_t* m_info = PDB_metals::lookup(b.name);
        if(dd < m_info->min_act_dist * m_info->min_act_dist) return true;
      // Van der Waals contact
      }else{
        my_float_t tmp = (a.vdw_radius + b.vdw_radius - ALLOWED_OVERLAP);
        if(dd < tmp*tmp) return true; 
        return false;
      }
      return false;
    }


    point_storage<atom_t> atoms; //!< Vector of atoms 
    atom_map_t atoms_by_num;    //!< Map (atom num, atom_vci)
    std::string name_of_file;   //!< Name of coordinate file which was (hopefully) read in.
    static const my_float_t ALLOWED_OVERLAP; //!< Tolerance for the sum of two atoms Van der Waals radii to be greater than the distance between them
    static const my_float_t min_hbond_dist;  //!< What is the closest allowable distance for two atoms participating in a hydrogen bond
    bool A_fail;     //!< Did the class fail? -- typically file I/O if it did
    static const std::string A_fname;  //!< Name of the source file
  };

  //! Simple rmsd 
  /*! For some reason g++ complains if this is defined in the source (.C) file
   */
  inline my_float_t 
  simple_rmsd(const CoordFile& f1, const CoordFile& f2)
  {
    atom_vci a1 = f1.atoms.begin();
    atom_vci a2 = f2.atoms.begin();
    my_float_t err = 0;
    my_float_t n = 0;
    for( ; a1 != f1.atoms.end() && a2 != f2.atoms.end(); ++a1, ++a2){
      if(a1->name != a2->name)
        std::cerr << "in RMSD a1 and a2 have different names\n";
      if(a1->name != H || a2->name != H){
        err += dist_squared(a1->pos, a2->pos);
        ++n;
      }
    }
    return std::sqrt(err / n);
  }

}
#endif
