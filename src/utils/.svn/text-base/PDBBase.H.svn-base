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

#ifndef PDBBASE_HEADER_INCLUDED
#define PDBBASE_HEADER_INCLUDED
#include <types.H>
#include <CoordFile.H>

namespace ASCbase{

  //! Type of PDB record/line
  typedef enum{
    PDB_NOT_ATOM_RECORD, 
    PDB_MODRES_RECORD,
    PDB_HET_RECORD,       //!< "HET   " record
    PDB_MODEL_RECORD,
    PDB_HYDROGEN_ATOM,
    PDB_HYDROGEN_HET,
    PDB_HEAVY_ATOM,
    PDB_HEAVY_HET
  }pdb_record_type;

  //! Base class for reading and writing PDB coordinate files
  class PDBBase : public CoordFile{
  public:

    //! Holds the information in one MODRES record (line)
    struct modres_record_type{
      std::string idCode;   //!< PDB ID code -- included for historical reasions?
      std::string resName;  //!< HET ID of the modified residue
      char chainID;         //!< Chain ID of the modified residue
      uint seqNum;          //!< Residue number of the modified residue
      char iCode;           //!< Insertion code 
      residue_type stdRes;  //!< Residue that was modified (1 of the 20 amino acids) 
      std::string comment;  //!< Any pertinant comments -- usually type of modfication

      bool
      contains(atom_vci atom) const
      {
        if(atom->res_num != seqNum) return false;
        if(atom->chainID != chainID) return false;
        if(atom->res_str != resName) return false;
        if(atom->iCode != iCode) return false;
        return true;
      }
    };
    typedef std::vector<modres_record_type> modres_record_vec;
    typedef modres_record_vec::iterator modres_record_vi;
    typedef modres_record_vec::const_iterator modres_record_vci;

    //! Holds the information in one HET record (line)
    struct het_record_type{
      std::string hetID;    //!< 3 character het identifier
      char chainID;         //!< Chain id
      uint seqNum;          //!< Residue number of the het group
      char iCode;           //!< Insertion code
      uint numHetAtoms;     //!< Number of heteroatoms in the het group
      std::string text;     //!< Additional text (column 31+)

      bool
      contains(atom_vci atom) const
      {
        if(seqNum != atom->res_num) return false;
        if(chainID != atom->chainID) return false;
        if(hetID != atom->res_str) return false;
        if(iCode != atom->iCode) return false;
        return true;
      }

      bool
      contains(const atom_t& atom) const
      {
        if(seqNum != atom.res_num) return false;
        if(chainID != atom.chainID) return false;
        if(hetID != atom.res_str) return false;
        if(iCode != atom.iCode) return false;
        return true;
      }
    };
    typedef std::vector<het_record_type> het_record_vec;
    typedef het_record_vec::iterator het_record_vi;
    typedef het_record_vec::const_iterator het_record_vci;

    //! Typically used to create a PDB file
    PDBBase();

    //! Read in a PDB file
    /*! After a successful call to this constructor, the heavy atoms will be
     * accessible by the CoordFile::begin()/end() functions.
     *
     * @param filename Path to the PDB file to read (heavy atoms)
     * @param verbosity Whether to send some comments to stdout or not
     */ 
    PDBBase(const std::string filename,
            const verbose_level_t verbosity = VERBOSE_SILENT);

    //! Do nothing destructor
    ~PDBBase();

    //! Write the atoms to file in PDB format
    /*!
     * @param filename Path to the file to truncate and write PDB atoms lines
     * @return True if file could be opened for writing, else false
     */
    bool 
    write(const std::string filename);

    //! entire existence is for debugging const pointers to pdb files
    const bool
    const_write(const std::string filename) const;

    //! Write as xyzr lines (coordinates + radius)
    /*!
     * @param out reference to a valid ostream
     * @param include_metals If false, ignore all metals in PDB file
     * @param waters Vector of water molecules to include
     * @return True if file could be opened for writing, else false
     */
    bool write_xyzr(std::ostream &out, bool include_metals,
                    const std::vector<std::string> &waters);

    //! Set the header
    /*!
     * @param h Header to write first in PDB file
     */
    void set_header(const std::string h);

    //! does nothing
    void report_stats(std::ostream& out) {;}

    //! True if the pdb atom name is a hydrogen atom, else false
    /*! 
     * Caution: you can't just test column 14 or you will remove
     * metals like Rh, Cd, Pd, Nd, Gd, Md...
     *
     * @param name The 4 columns used to specify the ATOM/HETATM name
     * @return True if is hydrogen, else false
     */
    static bool is_hydrogen(const std::string name);

    //! Determine if the line is an ATOM or HETATM record and if the line is
    //! for a hydrogen atom
    /*! 
     * @param line PDB line to test
     * @return PDB atom line type (not atom line, hetatm, atom or hydrogen)
     */
    static pdb_record_type get_record_type(const std::string line);

    modres_record_vci
    modres_records_begin() const
    { return A_modres_records.begin(); }

    modres_record_vci
    modres_records_end() const
    { return A_modres_records.end(); }

    het_record_vci
    het_records_begin() const
    { return A_het_records.begin(); }

    het_record_vci
    het_records_end() const
    { return A_het_records.end(); }

  protected:
    //! This is cheating, but .... I need a quick and dirty way to do this atm
    bool 
    read_data(std::ifstream &fin)
    {
      return read_data(fin, false);
    }

    //! Read in the heavy atom coordinate lines (HETATM and ATOM) from the path
    //! returned by CoordFile::name()
    /*!
     * @param fin Reference to an input pdb file stream 
     * @param ignore_altLocs If true, use only the first atomic position
     *        for each group of alternate locations
     * @return True if no error occured during reading, else false
     */
    bool read_data(std::ifstream &fin, const bool ignore_altLocs);

    //! Create a new het record for the given atom and add it to the het 
    //! records.
    /*!
     * A number of pdb files do not have HET entries for 1 or more of the 
     * HET groups
     */
    void
    add_het_record(atom_vci atom);

  private:

    //! Given a 4 char atom string and 3 char res_name, load the info into the 
    //! given atom structure
    /*! Given an atom string an residue string, load the Van Der Waals radius,
     * residue enum, atom enum and atom orbit into the given atom structure.
     *
     * @param atom_name 4 char atom string
     * @param res_name 3 char residue string
     * @param atom Pointer to an atom structure 
     * return True if the given residue and atom strings are in table, else 
     *        false
     */
    bool get_res_info(const std::string atom_name, const std::string res_name, 
                      atom_t* atom);

    //! Read a HET record and store it in the het record vector
    void
    read_het_record(std::string& het_line);
 
    //! Read a MODRES record and store it in the modres record vector
    void 
    read_modres_record(std::string& modres_line);

    het_record_vec A_het_records; //!< Vector of "HET   " records
    modres_record_vec A_modres_records;  //!< Vector of MODRES records
    std::string A_header;    //!< String holding the header info
    static const std::string A_fname;  //!< Name of source file
  };
  
}

#endif
