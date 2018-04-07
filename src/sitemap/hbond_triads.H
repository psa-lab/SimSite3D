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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/hbond_triads.H,v $
 * $Revision: 1.2 $
 * $Author: vanvoor4 $
 * $Date: 2007-09-24 15:35:06 $
 *
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2007/08/21 18:13:54  vanvoor4
 * initial checkin
 *
 *
 *
 *
 */

#ifndef HBOND_TRIADS_HEADER_INCLUDED
#define HBOND_TRIADS_HEADER_INCLUDED

#include <types.H>
#include <point.H>
#include <atom.H>

namespace SimSite3D{

  //! Data struct holding the local ideal hbond position parameters
  typedef struct{
    residue_type residue;   //!< Name of the residue
    atom_type nbr_A;        //!< Name of the first neighbor (carbon nbr)
    atom_type B;            //!< Name of the hbond atom
    atom_type nbr_C;        //!< Name of the second neighbor
    my_float_t bond_len;    //!< "Average" hbond length (D-A distance)
    interactionType type;   //!< interaction type to satisify the hbond atom
    int pt_num;             //!< Object number for the hbond atom
    my_float_t alpha;       //!< "In-plane" angle (rotation about Z axis)
    my_float_t beta;        //!< "out-of-plane" angle (rotation about X axis)
  }hbond_triad_t;
  
  const uint num_hbond_triads = 20;
  
  const hbond_triad_t hbond_triads[] = {
  {MAIN_CHAIN,C,O,N,3.0,DONOR, 0, 0.0, 0.0},
  {MAIN_CHAIN,C,N,O,3.0,ACCEPTOR, 0, 0.0, 0.0},
  
  {ARG,CZ,NH1,NH2,3.0,ACCEPTOR, 0, 0.0, 0.0},
  {ARG,CZ,NH2,NH1,3.0,ACCEPTOR, 0, 0.0, 0.0},
  {ARG,CZ,NE,CD,3.0,ACCEPTOR, 0, 0.0, 0.0},
  
  {ASN,CG,OD1,ND2,3.0,DONOR, 0, 0.0, 0.0},
  {ASN,CG,ND2,OD1,3.0,ACCEPTOR, 0, 0.0, 0.0},
  
  {ASP,CG,OD1,OD2,2.9,DONOR, 0, 0.0, 0.0},
  {ASP,CG,OD2,OD1,2.9,DONOR, 0, 0.0, 0.0},
  
  {GLN,CD,OE1,NE2,3.0,DONOR, 0, 0.0, 0.0},
  {GLN,CD,NE2,OE1,3.0,ACCEPTOR, 0, 0.0, 0.0},
  
  {GLU,CD,OE1,OE2,2.9,DONOR, 0, 0.0, 0.0},
  {GLU,CD,OE2,OE1,2.9,DONOR, 0, 0.0, 0.0},
  
  {HIS,CE1,ND1,NE2,3.0,DONEPTOR, 0, 0.0, 0.0},
  // Changed to DONOR -- Discussion with Leslie and Matt based on 
  // Jack Kyte’s Structure in Protein Chemistry pKa diagram for the 
  // high-probability states of His around pH 7
  {HIS,CE1,NE2,ND1,3.0,DONOR, 0, 0.0, 0.0},
  
  {LYS,CE,NZ,CD,2.9,ACCEPTOR, 0, 0.0, 0.0},
  
  {SER,CB,OG,CA,3.0,DONEPTOR, 0, 0.0, 0.0},
  
  {THR,CB,OG1,CA,2.9,DONEPTOR, 0, 0.0, 0.0},
  
  {TRP,CD1,NE1,CE2,3.0,ACCEPTOR, 0, 0.0, 0.0},
  
  {TYR,CZ,OH,CE2,2.9,DONEPTOR, 0, 0.0, 0.0}
  };

  const uint num_new_hbond_triads = 38;
  const hbond_triad_t new_hbond_triads[] = {
  {MAIN_CHAIN,C,O,N,3.0,DONOR, 1, 140.0, 0.0},
  {MAIN_CHAIN,C,O,N,3.0,DONOR, 2, 220.0, 0.0},
  {MAIN_CHAIN,C,N,O,3.0,ACCEPTOR, 1, 120.0, 0.0},
  
  {ARG,CZ,NH1,NH2,3.0,ACCEPTOR, 1, 120.0, 0.0},
  {ARG,CZ,NH1,NH2,3.0,ACCEPTOR, 2, 240.0, 0.0},
  {ARG,CZ,NH2,NH1,3.0,ACCEPTOR, 1, 120.0, 0.0},
  {ARG,CZ,NH2,NH1,3.0,ACCEPTOR, 2, 240.0, 0.0},
  {ARG,CZ,NE,CD,3.0,ACCEPTOR, 1, 120.0, 0.0},
  
  {ASN,CG,OD1,ND2,3.0,DONOR, 1, 140.0, 0.0},
  {ASN,CG,OD1,ND2,3.0,DONOR, 2, 240.0, 0.0},
  {ASN,CG,ND2,OD1,3.0,ACCEPTOR, 1, 120.0, 0.0},
  {ASN,CG,ND2,OD1,3.0,ACCEPTOR, 2, 240.0, 0.0},
  
  {ASP,CG,OD1,OD2,2.9,DONOR, 1, 140.0, 0.0},
  {ASP,CG,OD1,OD2,2.9,DONOR, 2, 240.0, 0.0},
  {ASP,CG,OD2,OD1,2.9,DONOR, 1, 140.0, 0.0},
  {ASP,CG,OD2,OD1,2.9,DONOR, 2, 240.0, 0.0},
  
  {GLN,CD,OE1,NE2,3.0,DONOR, 1, 140.0, 0.0},
  {GLN,CD,OE1,NE2,3.0,DONOR, 2, 240.0, 0.0},
  {GLN,CD,NE2,OE1,3.0,ACCEPTOR, 1, 120.0, 0.0},
  {GLN,CD,NE2,OE1,3.0,ACCEPTOR, 2, 240.0, 0.0},
  
  {GLU,CD,OE1,OE2,2.9,DONOR, 1, 140.0, 0.0},
  {GLU,CD,OE1,OE2,2.9,DONOR, 2, 240.0, 0.0},
  {GLU,CD,OE2,OE1,2.9,DONOR, 1, 120.0, 0.0},
  {GLU,CD,OE2,OE1,2.9,DONOR, 2, 240.0, 0.0},
  
  {HIS,CE1,ND1,NE2,3.0,DONEPTOR, 1, 120.0, 0.0},
  // Changed to DONOR -- Discussion with Leslie and Matt based on 
  // Jack Kyte’s Structure in Protein Chemistry pKa diagram for the 
  // high-probability states of His around pH 7
  {HIS,CE1,NE2,ND1,3.0,DONOR, 1, 120.0, 0.0},
  
  {LYS,CE,NZ,CD,2.9,ACCEPTOR, 1, 110.0, 0.0},
  {LYS,CE,NZ,CD,2.9,ACCEPTOR, 2, 110.0, 120.0},
  {LYS,CE,NZ,CD,2.9,ACCEPTOR, 3, 110.0, 240.0},
  
  {SER,CB,OG,CA,3.0,DONEPTOR, 1, 115.0, 0.0},
  {SER,CB,OG,CA,3.0,DONEPTOR, 2, 115.0, 120.0},
  {SER,CB,OG,CA,3.0,DONEPTOR, 3, 115.0, 240.0},
  
  {THR,CB,OG1,CA,2.9,DONEPTOR, 1, 115.0, 0.0},
  {THR,CB,OG1,CA,2.9,DONEPTOR, 2, 115.0, 120.0},
  {THR,CB,OG1,CA,2.9,DONEPTOR, 3, 115.0, 240.0},
  
  {TRP,CD1,NE1,CE2,3.0,ACCEPTOR, 1, 120.0, 0.0},
  
  {TYR,CZ,OH,CE2,2.9,DONEPTOR, 1, 120.0, 0.0},
  {TYR,CZ,OH,CE2,2.9,DONEPTOR, 2, 240.0, 0.0}
  };

}
#endif
