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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/hphob_triads.H,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2008-04-17 18:50:58 $
 *
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2008/03/31 17:49:03  vanvoor4
 * Added a number of sidechain atoms to hphob triads
 *
 * Revision 1.1  2007/09/26 14:54:27  vanvoor4
 * Initial checkin
 *
 *
 *
 */

#ifndef HPHOB_TRIADS_HEADER_INCLUDED
#define HPHOB_TRIADS_HEADER_INCLUDED

namespace SimSite3D{

  typedef struct{
    residue_type residue;  
    atom_type nbr_A;
    atom_type B;
    atom_type nbr_C;
  }hphob_triad_t;
  typedef std::map<residue_type, std::map<atom_type, const hphob_triad_t*> > 
    hphob_triad_map_t;
  typedef std::vector<const hphob_triad_t*> hphob_triad_vec;
  typedef hphob_triad_vec::const_iterator hphob_triad_vci;

  const uint num_hphob_triads = 52;
  const hphob_triad_t hphob_triads[] = {
    {ALA,CA ,CB ,C  },
    {ARG,CG ,CB ,CA },
    {ARG,CD ,CG ,CB },
    {ASP,CG ,CB ,CA },
    {CYS,CA ,CB ,SG },
    {CYS,CB ,SG ,CA },
    {GLN,CG ,CB ,CA },
    {GLN,CD ,CG ,CB },
    {GLU,CG ,CB ,CA },
    {GLU,CD ,CG ,CB },
    {HIS,CA, CB, CG },
    {ILE,CA ,CB ,CG1},
    {ILE,CB ,CG1,CD1},
    {ILE,CB ,CG2,CG1},
    {ILE,CG1,CD1,CB },
    {LEU,CG ,CB ,CA },
    {LEU,CB ,CG ,CA },
    {LEU,CG ,CD1,CD2},
    {LEU,CG ,CD2,CD1},
    {LYS,CG ,CB ,CA },
    {LYS,CD ,CG ,CB },
    {LYS,CE ,CD ,CG },
    {MET,CA ,CB ,CG },
    {MET,CB ,CG ,SD },
    {MET,CE ,SD ,CG },
    {MET,SD ,CE ,CG },
    {PCA,CG ,CB ,CA },
    {PCA,CD ,CG ,CB },
    {PHE,CG ,CB ,CA },
    {PHE,CD1,CG ,CD2},
    {PHE,CE1,CD1,CG },
    {PHE,CE2,CD2,CG },
    {PHE,CZ ,CE1,CD1},
    {PHE,CZ ,CE2,CD2},
    {PHE,CE1,CZ ,CE2},
    {PRO,CG ,CB ,CA },
    {PRO,CD ,CG ,CB },
    {TRP,CG ,CB ,CA },
    {TRP,CD1,CG ,CD2},
    {TRP,CE2,CD2,CE3},
    {TRP,CZ3,CE2,CD2},
    {TRP,CH2,CZ2,CE2},
    {TRP,CH2,CZ3,CE3},
    {TRP,CZ2,CH2,CZ3},
    {TYR,CG ,CB ,CA },
    {TYR,CE1,CD1,CG },
    {TYR,CE2,CD2,CG },
    {TYR,CZ ,CE1,CD1},
    {TYR,CZ ,CE2,CD2},
    {VAL,CA ,CB ,CG1},
    {VAL,CB ,CG1,CG2},
    {VAL,CB ,CG2,CG1}};
}
#endif
