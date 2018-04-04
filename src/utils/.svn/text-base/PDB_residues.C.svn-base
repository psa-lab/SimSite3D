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

/*
 * $Source: /psa/share/repository/pfizer_proj/src/basics/PDB_residues.C,v $
 * $Revision: 1.5 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 20:57:08 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.4  2008/07/28 15:21:16  vanvoor4
 * Changed to be a static class that is initialized upon call.
 *
 * Revision 1.3  2008/04/17 18:41:34  vanvoor4
 * Added another residue to the enum to distinguish between
 * errors/unknown residues and the NULL residue.
 *
 * Revision 1.2  2008/03/24 15:29:52  vanvoor4
 * Added a check for water residues -- do not want to have main chain
 * atoms for a water residu.
 *
 * Revision 1.1  2008/02/26 19:08:31  vanvoor4
 * Initial checkin -- need a compiled file to setup constant values
 * for the standard PDB residues
 *
 * 
 */

#include <iostream>
#include <PDB_residues.H>

using namespace ASCbase;

atom_table_t PDB_residues::residue_table;
atom_enum_table_t PDB_residues::residue_type_table;
const uint PDB_residues::A_res_array_size = 102;
const pdb_atom_info_t PDB_residues::A_res_array[] = {
  {" N  ", "MAIN_CHAIN",  N  , MAIN_CHAIN, 1.43, 1.70, AMIDE, ALPHA, DONOR, 
    0.0 },
  {" CA ", "MAIN_CHAIN",  CA , MAIN_CHAIN, 1.73, 2.00, SP3, ALPHA, NOTHING, 
    0.0 },
  {" C  ", "MAIN_CHAIN",  C  , MAIN_CHAIN, 1.62, 1.74, SP2, ALPHA, NOTHING, 
    0.0 },
  {" O  ", "MAIN_CHAIN",  O  , MAIN_CHAIN, 1.26, 1.40, SP2, ALPHA, ACCEPTOR, 
    0.0 },
  {" OXT", "MAIN_CHAIN",  OXT, MAIN_CHAIN, 1.24, 1.40, CO2, ALPHA, ACCEPTOR, 
    -1.0 },
  {" CB ", "ALA",  CB , ALA, 1.66, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CB ", "ARG",  CB , ARG, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "ARG",  CG , ARG, 1.68, 2.00, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CD ", "ARG",  CD , ARG, 1.62, 2.00, SP3, DELTA, NOTHING, 0.0 },
  {" NE ", "ARG",  NE , ARG, 1.34, 1.70, PL3, EPSILON, DONOR, 1.0 },
  {" CZ ", "ARG",  CZ , ARG, 1.62, 1.74, CAT, ZETA, NOTHING, 0.0 },
  {" NH1", "ARG",  NH1, ARG, 1.34, 1.80, PL3, ETA, DONOR, 1.0 },
  {" NH2", "ARG",  NH2, ARG, 1.34, 1.80, PL3, ETA, DONOR, 1.0 },
  {" CB ", "ASN",  CB , ASN, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "ASN",  CG , ASN, 1.67, 1.74, SP2, GAMMA, NOTHING, 0.0},
  {" OD1", "ASN",  OD1, ASN, 1.20, 1.40, SP2, DELTA, ACCEPTOR, 0.0 },
  {" ND2", "ASN",  ND2, ASN, 1.31, 1.80, AMIDE, DELTA, DONOR, 0.0 },
  {" AD1", "ASN",  AD1, ASN, 1.25, 1.40, AMBIG, DELTA, DONEPTOR, 0.0 },
  {" AD2", "ASN",  AD2, ASN, 1.25, 1.40, AMBIG, DELTA, DONEPTOR, 0.0 },
  {" CB ", "ASP",  CB , ASP, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "ASP",  CG , ASP, 1.65, 1.74, SP2, GAMMA, NOTHING, 0.0 },
  {" OD1", "ASP",  OD1, ASP, 1.24, 1.40, SP2, DELTA, ACCEPTOR, -1.0 },
  {" OD2", "ASP",  OD2, ASP, 1.24, 1.40, SP2, DELTA, ACCEPTOR, -1.0 },
  {" CB ", "CYS",  CB , CYS, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" SG ", "CYS",  SG , CYS, 1.54, 1.80, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CB ", "GLN",  CB , GLN, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "GLN",  CG , GLN, 1.62, 2.00, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CD ", "GLN",  CD , GLN, 1.67, 1.74, SP2, DELTA, NOTHING, 0.0 },
  {" OE1", "GLN",  OE1, GLN, 1.20, 1.40, SP2, EPSILON, ACCEPTOR, 0.0 },
  {" NE2", "GLN",  NE2, GLN, 1.31, 1.80, AMIDE, EPSILON, DONOR, 0.0 },
  {" AE1", "GLN",  AE1, GLN, 1.25, 1.40, AMBIG, EPSILON, DONEPTOR, 0.0 },
  {" AE2", "GLN",  AE2, GLN, 1.25, 1.40, AMBIG, EPSILON, DONEPTOR, 0.0 },
  {" CB ", "GLU",  CB , GLU, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "GLU",  CG , GLU, 1.62, 2.00, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CD ", "GLU",  CD , GLU, 1.65, 1.74, SP2, DELTA, NOTHING, 0.0 },
  {" OE1", "GLU",  OE1, GLU, 1.24, 1.40, CO2, EPSILON, ACCEPTOR, -1.0 },
  {" OE2", "GLU",  OE2, GLU, 1.24, 1.40, CO2, EPSILON, ACCEPTOR, -1.0 },
  {" CB ", "HIS",  CB , HIS, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "HIS",  CG , HIS, 1.62, 1.74, SP2, GAMMA, NOTHING, 0.0 },
  {" ND1", "HIS",  ND1, HIS, 1.38, 1.54, PL3, DELTA, DONEPTOR, 1.0 },
  {" CD2", "HIS",  CD2, HIS, 1.50, 1.86, SP2, DELTA, NOTHING, 0.0 },
  {" CE1", "HIS",  CE1, HIS, 1.50, 1.86, SP2, EPSILON, NOTHING, 0.0 },
  // Changed to DONOR -- Discussion with Leslie and Matt based on 
  // Jack Kyte’s Structure in Protein Chemistry pKa diagram for the 
  // high-probability states of His around pH 7
  {" NE2", "HIS",  NE2, HIS, 1.38, 1.70, PL3, EPSILON, DONOR, 1.0 },
  {" CB ", "ILE",  CB , ILE, 1.82, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG1", "ILE",  CG1, ILE, 1.68, 2.00, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CG2", "ILE",  CG2, ILE, 1.68, 2.00, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CD1", "ILE",  CD1, ILE, 1.66, 2.00, SP3, DELTA, HYDROPHOB, 0.0 },
  {" CB ", "LEU",  CB , LEU, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "LEU",  CG , LEU, 1.82, 2.00, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CD1", "LEU",  CD1, LEU, 1.66, 2.00, SP3, DELTA, HYDROPHOB, 0.0 },
  {" CD2", "LEU",  CD2, LEU, 1.66, 2.00, SP3, DELTA, HYDROPHOB, 0.0 },
  {" CB ", "LYS",  CB , LYS, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "LYS",  CG , LYS, 1.68, 2.00, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CD ", "LYS",  CD , LYS, 1.68, 2.00, SP3, DELTA, HYDROPHOB, 0.0 },
  {" CE ", "LYS",  CE , LYS, 1.62, 2.00, SP3, EPSILON, NOTHING, 0.0 },
  {" NZ ", "LYS",  NZ , LYS, 1.22, 2.00, SP4, ZETA, DONOR, 1.0 },
  {" CB ", "MET",  CB , MET, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "MET",  CG , MET, 1.68, 2.00, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" SD ", "MET",  SD , MET, 1.67, 1.80, SP3, DELTA, HYDROPHOB, 0.0 },
  {" CE ", "MET",  CE , MET, 1.66, 2.00, SP3, EPSILON, HYDROPHOB, 0.0 },
  {" CB ", "PCA",  CB , PCA, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "PCA",  CG , PCA, 1.69, 2.00, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CD ", "PCA",  CD , PCA, 1.69, 2.00, SP3, DELTA, NOTHING, 0.0 },
  {" OE ", "PCA",  OE , PCA, 1.24, 1.40, CO2, EPSILON, NOTHING, 0.0 },
  {" CB ", "PHE",  CB , PHE, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "PHE",  CG , PHE, 1.62, 1.74, AR , GAMMA, HYDROPHOB, 0.0 },
  {" CD1", "PHE",  CD1, PHE, 1.62, 1.86, AR , DELTA, HYDROPHOB, 0.0 },
  {" CD2", "PHE",  CD2, PHE, 1.62, 1.86, AR , DELTA, HYDROPHOB, 0.0 },
  {" CE1", "PHE",  CE1, PHE, 1.62, 1.86, AR , EPSILON, HYDROPHOB, 0.0 },
  {" CE2", "PHE",  CE2, PHE, 1.62, 1.86, AR , EPSILON, HYDROPHOB, 0.0 },
  {" CZ ", "PHE",  CZ , PHE, 1.62, 1.86, AR , ZETA, HYDROPHOB, 0.0 },
  {" CB ", "PRO",  CB , PRO, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "PRO",  CG , PRO, 1.68, 2.00, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CD ", "PRO",  CD , PRO, 1.68, 2.00, SP3, DELTA, NOTHING, 0.0 },
  {" CB ", "SER",  CB , SER, 1.69, 2.00, SP3, BETA, NOTHING, 0.0 },
  {" OG ", "SER",  OG , SER, 1.30, 1.60, SP3, GAMMA, DONEPTOR, 0.0 },
  {" CB ", "THR",  CB , THR, 1.82, 2.00, SP3, BETA, NOTHING, 0.0 },
  {" OG1", "THR",  OG1, THR, 1.30, 1.60, SP3, GAMMA, DONEPTOR, 0.0 },
  {" CG2", "THR",  CG2, THR, 1.66, 1.74, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CB ", "TRP",  CB , TRP, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "TRP",  CG , TRP, 1.62, 1.74, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CD1", "TRP",  CD1, TRP, 1.62, 1.86, SP2, DELTA, NOTHING, 0.0 },
  {" CD2", "TRP",  CD2, TRP, 1.62, 1.74, AR , DELTA, HYDROPHOB, 0.0 },
  {" NE1", "TRP",  NE1, TRP, 1.35, 1.70, PL3, EPSILON, DONOR, 0.0 },
  {" CE2", "TRP",  CE2, TRP, 1.62, 1.74, AR , EPSILON, HYDROPHOB, 0.0 },
  {" CE3", "TRP",  CE3, TRP, 1.62, 1.86, AR , EPSILON, NOTHING, 0.0 },
  {" CZ2", "TRP",  CZ2, TRP, 1.62, 1.86, AR , ZETA, HYDROPHOB, 0.0 },
  {" CZ3", "TRP",  CZ3, TRP, 1.62, 1.86, AR , ZETA, HYDROPHOB, 0.0 },
  {" CH2", "TRP",  CH2, TRP, 1.62, 1.86, AR , ETA, HYDROPHOB, 0.0 },
  {" CB ", "TYR",  CB , TYR, 1.69, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG ", "TYR",  CG , TYR, 1.62, 1.74, AR , GAMMA, HYDROPHOB, 0.0 },
  {" CD1", "TYR",  CD1, TYR, 1.62, 1.86, AR , DELTA, HYDROPHOB, 0.0 },
  {" CD2", "TYR",  CD2, TYR, 1.62, 1.86, AR , DELTA, HYDROPHOB, 0.0 },
  {" CE1", "TYR",  CE1, TYR, 1.62, 1.86, AR , DELTA, HYDROPHOB, 0.0 },
  {" CE2", "TYR",  CE2, TYR, 1.62, 1.86, AR , DELTA, HYDROPHOB, 0.0 },
  {" CZ ", "TYR",  CZ , TYR, 1.62, 1.74, AR , ZETA, NOTHING, 0.0 },
  {" OH ", "TYR",  OH , TYR, 1.30, 1.60, SP3, ETA, DONEPTOR, 0.0 },
  {" CB ", "VAL",  CB , VAL, 1.82, 2.00, SP3, BETA, HYDROPHOB, 0.0 },
  {" CG1", "VAL",  CG1, VAL, 1.66, 2.00, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" CG2", "VAL",  CG2, VAL, 1.66, 2.00, SP3, GAMMA, HYDROPHOB, 0.0 },
  {" O  ", "HOH",  O  , HOH, 1.36, 1.60, SP3, ALPHA, DONEPTOR, 0.0 },
  {" OH2", "HOH",  O  , HOH, 1.36, 1.60, SP3, ALPHA, DONEPTOR, 0.0 }
};

const pdb_atom_info_t PDB_residues::Pro_main_chain_N =
  {" N  ", "MAIN_CHAIN",  N  , MAIN_CHAIN, 1.43, 1.70, AMIDE, ALPHA, NOTHING, 
    0.0 };

const uint PDB_residues::A_res_conv_array_size = 27;
const residue_conv_type PDB_residues::A_res_conv_array[] = {
  {NULL_RESIDUE, "NULL_RESIDUE"},
  {UNKNOWN_RESIDUE, "UNKNOWN_RESIDUE"},
  {ALA, "ALA"},
  {ARG, "ARG"},
  {ASN, "ASN"},
  {ASP, "ASP"},
  {CYS, "CYS"},
  {GLN, "GLN"},
  {GLU, "GLU"},
  {GLY, "GLY"},
  {HIS, "HIS"},
  {ILE, "ILE"},
  {LEU, "LEU"},
  {LYS, "LYS"},
  {MET, "MET"},
  {PCA, "PCA"},
  {PHE, "PHE"},
  {PRO, "PRO"},
  {SER, "SER"},
  {THR, "THR"},
  {TRP, "TRP"},
  {TYR, "TYR"},
  {VAL, "VAL"},
  {HOH, "HOH"},
//  {TPO, "TPO"}, -- not handled wrt hbonds yet
//  {PTR, "PTR"}, -- not handled wrt hbonds yet
//  {ACE, "ACE"}, -- not handled wrt hbonds yet
  {MAIN_CHAIN, "MAIN_CHAIN"},
  {PDB_METAL, "PDB_METAL"},
  {RESIDUE_ENUM_END, "RESIDUE_ENUM_END"}
};

const uint PDB_residues::A_atom_conv_array_size = 42;
const atom_conv_type PDB_residues::A_atom_conv_array[] = {
  {  AD1, " AD1" },
  {  AD2, " AD2" },
  {  AE1, " AE1" },
  {  AE2, " AE2" },
  {  C  , " C  " },
  {  CA , " CA " },
  {  CB , " CB " },
  {  CD , " CD " },
  {  CD1, " CD1" },
  {  CD2, " CD2" },
  {  CE , " CE " },
  {  CE1, " CE1" },
  {  CE2, " CE2" },
  {  CE3, " CE3" },
  {  CG , " CG " },
  {  CG1, " CG1" },
  {  CG2, " CG2" },
  {  CH2, " CH2" },
  {  CZ , " CZ " },
  {  CZ2, " CZ2" },
  {  CZ3, " CZ3" },
  {  N  , " N  " },
  {  ND1, " ND1" },
  {  ND2, " ND2" },
  {  NE , " NE " },
  {  NE1, " NE1" },
  {  NE2, " NE2" },
  {  NH1, " NH1" },
  {  NH2, " NH2" },
  {  NZ , " NZ " },
  {  O  , " O  " },
  {  OD1, " OD1" },
  {  OD2, " OD2" },
  {  OE , " OE " },
  {  OE1, " OE1" },
  {  OE2, " OE2" },
  {  OG , " OG " },
  {  OG1, " OG1" },
  {  OH , " OH " },
  {  OXT, " OXT" },
  {  SD , " SD " },
  {  SG , " SG " },
};

std::map<residue_type, std::string> PDB_residues::A_residue_to_string;
std::map<std::string, residue_type> PDB_residues::A_string_to_residue;
std::map<atom_type, std::string> PDB_residues::A_atom_to_string;

const size_t PDB_residues::A_hphob_val_array_size = 167;
const pdb_hphob_type PDB_residues::A_hphob_val_array[] = {
  { ALA, CA , -163 },
  { ALA, CB ,   -3 },
  { ALA, C  , -228 },
  { ALA, O  ,  349 },
  { ALA, N  ,  120 },
  { ARG, CA , -196 },
  { ARG, CB , -202 },
  { ARG, C  , -223 },
  { ARG, O  ,  371 },
  { ARG, N  ,  235 },
  { ARG, CD ,  -63 },
  { ARG, CG , -150 },
  { ARG, CZ , -209 },
  { ARG, NE ,    2 },
  { ARG, NH1,  205 },
  { ARG, NH2,  210 },
  { ASN, CA , -214 },
  { ASN, CB , -147 },
  { ASN, C  , -218 },
  { ASN, O  ,  213 },
  { ASN, N  ,  132 },
  { ASN, CG , -195 },
  { ASN, ND2,  229 },
  { ASN, OD1,  224 },
  { ASP, CA , -224 },
  { ASP, CB , -153 },
  { ASP, C  , -229 },
  { ASP, O  ,  268 },
  { ASP, N  ,   66 },
  { ASP, CG , -198 },
  { ASP, OD1,  311 },
  { ASP, OD2,  311 },
  { CYS, CA , -199 },
  { CYS, CB , -135 },
  { CYS, C  , -221 },
  { CYS, O  ,  348 },
  { CYS, N  ,  165 },
  { CYS, SG ,  -96 },
  { GLN, CA , -195 },
  { GLN, CB , -171 },
  { GLN, C  , -226 },
  { GLN, O  ,  267 },
  { GLN, N  ,  180 },
  { GLN, CD , -215 },
  { GLN, CG , -119 },
  { GLN, NE2,  181 },
  { GLN, OE1,  218 },
  { GLU, CA , -195 },
  { GLU, CB , -165 },
  { GLU, C  , -229 },
  { GLU, O  ,  278 },
  { GLU, N  ,  108 },
  { GLU, CD , -221 },
  { GLU, CG , -126 },
  { GLU, OE1,  246 },
  { GLU, OE2,  234 },
  { GLY, CA ,  -96 },
  { GLY, C  , -207 },
  { GLY, O  ,  286 },
  { GLY, N  ,   77 },
  { HIS, CA , -212 },
  { HIS, CB , -136 },
  { HIS, C  , -213 },
  { HIS, O  ,  207 },
  { HIS, N  ,  124 },
  { HIS, CD2,  -95 },
  { HIS, CE1,  -70 },
  { HIS, CG , -175 },
  { HIS, ND1,  206 },
  { HIS, NE2,  118 },
  { ILE, CA , -211 },
  { ILE, CB , -217 },
  { ILE, C  , -235 },
  { ILE, O  ,  310 },
  { ILE, N  ,  145 },
  { ILE, CD1, -148 },
  { ILE, CG1, -204 },
  { ILE, CG2, -111 },
  { LEU, CA , -225 },
  { LEU, CB , -200 },
  { LEU, C  , -235 },
  { LEU, O  ,  313 },
  { LEU, N  ,  149 },
  { LEU, CD1, -110 },
  { LEU, CD2, -118 },
  { LEU, CG , -215 },
  { LYS, CA , -223 },
  { LYS, CB , -175 },
  { LYS, C  , -226 },
  { LYS, O  ,  239 },
  { LYS, N  ,   36 },
  { LYS, CD , -154 },
  { LYS, CE , -106 },
  { LYS, CG , -159 },
  { LYS, NZ ,  241 },
  { MET, CA , -205 },
  { MET, CB , -103 },
  { MET, C  , -205 },
  { MET, O  ,  356 },
  { MET, N  ,  251 },
  { MET, CE ,  -27 },
  { MET, CG , -140 },
  { MET, SD , -220 },
  { PHE, CA , -221 },
  { PHE, CB , -206 },
  { PHE, C  , -235 },
  { PHE, O  ,  268 },
  { PHE, N  ,  171 },
  { PHE, CD1, -168 },
  { PHE, CD2, -159 },
  { PHE, CE1, -136 },
  { PHE, CE2, -187 },
  { PHE, CG , -235 },
  { PHE, CZ , -154 },
  { PRO, CA , -186 },
  { PRO, CB , -103 },
  { PRO, C  , -200 },
  { PRO, O  ,  253 },
  { PRO, N  , -219 },
  { PRO, CD ,  -89 },
  { PRO, CG ,  -71 },
  { SER, CA , -213 },
  { SER, CB ,  -85 },
  { SER, C  , -235 },
  { SER, O  ,  256 },
  { SER, N  ,   10 },
  { SER, OG ,  333 },
  { THR, CA , -214 },
  { THR, CB , -181 },
  { THR, C  , -235 },
  { THR, O  ,  366 },
  { THR, N  ,   57 },
  { THR, CG2,  -42 },
  { THR, OG1,  301 },
  { TRP, CA , -189 },
  { TRP, CB , -186 },
  { TRP, C  , -204 },
  { TRP, O  ,  334 },
  { TRP, N  ,  194 },
  { TRP, CD1,  -86 },
  { TRP, CD2, -168 },
  { TRP, CE2, -135 },
  { TRP, CE3, -208 },
  { TRP, CG , -235 },
  { TRP, CH2, -120 },
  { TRP, CZ2, -163 },
  { TRP, CZ3, -164 },
  { TRP, NE1,  114 },
  { TYR, CA , -229 },
  { TYR, CB , -166 },
  { TYR, C  , -161 },
  { TYR, O  ,  316 },
  { TYR, N  ,  139 },
  { TYR, CD1, -152 },
  { TYR, CD2, -184 },
  { TYR, CE1, -167 },
  { TYR, CE2, -185 },
  { TYR, CG , -213 },
  { TYR, CZ , -199 },
  { TYR, OH ,  400 },
  { VAL, CA , -222 },
  { VAL, CB , -235 },
  { VAL, C  , -235 },
  { VAL, O  ,  299 },
  { VAL, N  ,  247 },
  { VAL, CG1, -150 },
  { VAL, CG2, -129 }
};
hphob_tbl_res_lvl PDB_residues::A_hydrophobic_value;

residue_type
PDB_residues::string_to_residue(const std::string res_name)
{ 
  if(A_string_to_residue.empty()) build_res_conv_maps(); 
  std::map<std::string, residue_type>::const_iterator my_iter;
  my_iter = A_string_to_residue.find(res_name);
  if(my_iter != A_string_to_residue.end()) return my_iter->second;
  return UNKNOWN_RESIDUE;
}

std::string
PDB_residues::residue_to_string(const residue_type res)
{
  if(A_residue_to_string.empty()) build_res_conv_maps();
  std::map<residue_type, std::string>::const_iterator my_iter;
  my_iter = A_residue_to_string.find(res);
  if(my_iter != A_residue_to_string.end()) return my_iter->second;
  return "UNKNOWN_RESIDUE";
}

std::string
PDB_residues::atom_to_string(const atom_type atom)
{
  if(A_atom_to_string.empty()) build_atom_conv_map();
  std::map<atom_type, std::string>::const_iterator my_iter;
  my_iter = A_atom_to_string.find(atom);
  if(my_iter != A_atom_to_string.end()) return my_iter->second; 
  return "UNKNOWN_ATOM";
}

void
PDB_residues::build_res_tables()
{
  std::string prev_res_str = "GLY";
  residue_type prev_res = GLY; 
  for(uint j = 0; j < 5; ++j){
    residue_table[A_res_array[j].pdb_atom_str][prev_res_str] = A_res_array + j;
    residue_type_table[A_res_array[j].atom][prev_res] = A_res_array + j;
  }

  // skip the first 5 as they represent the main chain atoms
  for(uint i = 5; i < A_res_array_size; ++i){
    residue_type cur_res = A_res_array[i].res;
    std::string cur_res_str = A_res_array[i].pdb_res_str;
    // When we get a change in residue, add in the main chain atoms
    // I prefer to have a larger table than check atom types during reading
    // of pdb structures
    if(prev_res != cur_res){
      prev_res = cur_res;
      if(cur_res != HOH)
        for(uint j = 0; j < 5; ++j){
          // Need to check if PCA is cyclic
          if(j == 0 && cur_res == PRO){
            residue_table[" N  "]["PRO"] = &Pro_main_chain_N;
            residue_type_table[N][PRO] = &Pro_main_chain_N;
          }else{
            residue_table[A_res_array[j].pdb_atom_str][cur_res_str] = 
              A_res_array + j;
            residue_type_table[A_res_array[j].atom][cur_res] = A_res_array + j; 
          }
        }
    }

    residue_table[A_res_array[i].pdb_atom_str][cur_res_str] = A_res_array + i;
    residue_type_table[A_res_array[i].atom][cur_res] = A_res_array + i;
  }
}

void 
PDB_residues::build_res_conv_maps()
{
  for(uint i = 0; i < A_res_conv_array_size; ++i){
    A_residue_to_string[A_res_conv_array[i].residue] = A_res_conv_array[i].name;
    A_string_to_residue[A_res_conv_array[i].name] = A_res_conv_array[i].residue;
  }
}

void
PDB_residues::build_atom_conv_map()
{
  for(uint i = 0; i < A_atom_conv_array_size; ++i)
    A_atom_to_string[A_atom_conv_array[i].atom] = A_atom_conv_array[i].name;
}

void
PDB_residues::build_atom_hphob_map()
{
  const pdb_hphob_type *ptr = A_hphob_val_array;
  for(size_t i = 0; i < A_hphob_val_array_size; ++i, ++ptr)
    A_hydrophobic_value[ptr->residue][ptr->atom] = ptr->hydrophobicity;
}
