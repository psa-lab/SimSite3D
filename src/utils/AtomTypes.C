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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/AtomTypes.C,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2008-05-13 17:06:55 $
 * 
 * $Log: not supported by cvs2svn $
 * 
 */

#include <AtomTypes.H>

using namespace SimSite3D;

std::multimap<atom_type, const atom_info_t*> AtomTypes::vdw_table;

const uint AtomTypes::A_atoms_array_size = 81;
const atom_info_t AtomTypes::A_atoms_array[] = {
  // give the null atom a large radius so that problems ensue if it is used
  { NULL_ATOM_TYPE, DEFAULT_ORBIT, 1000.0},   
  { AG, DEFAULT_ORBIT, 1.550 },
  { AU, DEFAULT_ORBIT, 1.900 },
  { AL, DEFAULT_ORBIT, 1.170 },
  { AS, DEFAULT_ORBIT, 1.700 },
  { B , DEFAULT_ORBIT, 1.700 },
  { BA, DEFAULT_ORBIT, 3.100 },
  { BE, DEFAULT_ORBIT, 1.900 },
  { BI, DEFAULT_ORBIT, 1.800 },
  { BR, DEFAULT_ORBIT, 1.950 },
  { C , SP1, 1.7 },
  { C , SP2, 1.7 },
  { C , SP3, 1.7 },
  { C , CAT, 1.7 },
  { C , AR , 1.7 },
  { C , DEFAULT_ORBIT, 1.7 },
  { CA, SP3, 1.7 },
  { CALCIUM, DEFAULT_ORBIT, 1.600 },
  { CADMIUM, DEFAULT_ORBIT, 2.300 }, 
  { CL, DEFAULT_ORBIT, 2.030 },
  { CO, DEFAULT_ORBIT, 1.130 },
  { CR, DEFAULT_ORBIT, 1.130 },
  { CS, DEFAULT_ORBIT, 3.500 },
  { CU, DEFAULT_ORBIT, 1.900 },
  { D , DEFAULT_ORBIT, 1.000 },  // D could well have larger radius than H
  { ER, DEFAULT_ORBIT, 1.590 },  // Erbium
  { F , DEFAULT_ORBIT, 1.550 },
  { FE, DEFAULT_ORBIT, 1.950 },
  { GA, DEFAULT_ORBIT, 1.550 },
  { GE, DEFAULT_ORBIT, 2.000 },
  { H , DEFAULT_ORBIT, 1.000 },
  { HG, DEFAULT_ORBIT, 2.000 },
  { I , DEFAULT_ORBIT, 2.350 },
  { IN, DEFAULT_ORBIT, 1.900 },
  { K , DEFAULT_ORBIT, 1.600 },
  { LA, DEFAULT_ORBIT, 1.830 },
  { LI, DEFAULT_ORBIT, 1.170 },
  { MG, DEFAULT_ORBIT, 1.170 },
  { MN, DEFAULT_ORBIT, 1.190 },
  { MO, DEFAULT_ORBIT, 1.750 },
  { N , SP1, 1.55 },
  { N , SP2, 1.55 },
  { N , SP3, 1.55 },
  { N , SP4, 1.55 },
  { N , PL3, 1.55 },
  { N , AR , 1.55 },
  { N , AMIDE, 1.55 },
  { N , DEFAULT_ORBIT, 1.55 },
  { NA, DEFAULT_ORBIT, 1.600 },
  { NI, DEFAULT_ORBIT, 1.240 },
  { O , SP3, 1.6 },
  { O , SP2, 1.6 },
  { O , CO2, 1.6 },
  { O , DEFAULT_ORBIT, 1.6 },
  { OS, DEFAULT_ORBIT, 2.3 },
  { P , SP3, 2.100 },  
  { P , DEFAULT_ORBIT, 2.100 },
  { PD, DEFAULT_ORBIT, 1.440 },
  { RB, DEFAULT_ORBIT, 3.200 },
  { RE, DEFAULT_ORBIT, 1.300 },
  { RH, DEFAULT_ORBIT, 1.220 },
  { RU, DEFAULT_ORBIT, 1.200 },
  { S , SP3, 1.8 },
  { S , SP2, 1.8 },
  { S , O_ORBITAL, 1.8 },
  { S , O2 , 1.8 },
  { S , DEFAULT_ORBIT, 1.8 },
  { SB, DEFAULT_ORBIT, 1.750 },
  { SE, DEFAULT_ORBIT, 1.850 },
  { SI, DEFAULT_ORBIT, 1.900 },
  { SN, DEFAULT_ORBIT, 2.300 },
  { SR, DEFAULT_ORBIT, 2.700 },
  { TE, DEFAULT_ORBIT, 2.000 },
  { TI, DEFAULT_ORBIT, 1.950 },
  { TL, DEFAULT_ORBIT, 2.200 },
  { U , DEFAULT_ORBIT, 1.750 },
  { V , DEFAULT_ORBIT, 1.060 },
  { Y_, DEFAULT_ORBIT, 1.610 },
  { ZN, DEFAULT_ORBIT, 2.200 },
  { ZR, DEFAULT_ORBIT, 1.420 },
  { UNKNOWN_ATOM, UNKNOWN_ORBIT, 2.000 }
};

void 
AtomTypes::build_vdw_table()
{
  if(vdw_table.size()) return;

  typedef std::pair<atom_type, const atom_info_t*> atm_inf_pair;
  for(uint i = 0; i < A_atoms_array_size; ++i){
    atom_type a = A_atoms_array[i].atom;
    vdw_table.insert(atm_inf_pair(a, &(A_atoms_array[i])));
  }
}
