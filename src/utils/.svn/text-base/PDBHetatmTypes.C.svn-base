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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/PDBHetatmTypes.C,v $
 * $Revision: 1.2 $
 * $Author: vanvoor4 $
 * $Date: 2008-07-29 18:06:41 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2008/05/13 17:13:54  vanvoor4
 * Initial checkin --
 * vdw radii are from table XI in Li & Nussinov
 *
 * 
 */

#include <PDBHetatmTypes.H>

using namespace SimSite3D;

std::map<std::string, atom_type> PDBHetatmTypes::A_string_to_atom_type;

const uint PDBHetatmTypes::A_num_atom_strings = 46;
const atom_conv_type PDBHetatmTypes::A_atom_strings[] = {
  { AG, "AG  " },  // silver
  { AL, "AL  " },  // aluminium
  { AS, "AS  " },  // arsenic
  { BR, "BR  " },  // bromine
  { C , " C  " },  // carbon
  { CA, " CA " },  // alpha carbon
  { CALCIUM, "CA  " }, // calcium
  { CADMIUM, "CD  " }, // cadmium
  { CL, "CL  " }, 
  { CO, "CO  " }, 
  { CU, "CU  " },
  { CR, "CR  " }, 
  { D , " D  " }, 
  { ER, "ER  " }, // eridium
  { F , " F  " }, 
  { FE, "FE  " }, // iron
  { GA, "GA  " },
  { GE, "GE  " },
  { H , " H  " },
  { I , " I  " },
  { K , " K  " },
  { LI, "LI  " }, 
  { MG, "MG  " },
  { MN, "MN  " },
  { MO, "MO  " },
  { N , " N  " },
  { NA, "NA  " },
  { NI, "NI  " },
  { O , " O  " },
  { OS, "OS  " }, // osmium
  { P , " P  " },
  { PD, "PD  " }, // should this be pb instead?
  { RB, "RB  " },
  { RE, "RE  " }, 
  { RH, "RH  " },
  { RU, "RU  " },
  { S , " S  " },
  { SE, "SE  " }, 
  { SI, "SI  " },
  { TI, "TI  " },
  { TL, "TL  " },  // thallium
  { U , " U  " },
  { V , " V  " },
  { Y_, " Y_ " }, 
  { ZN, "ZN  " },
  { ZR, "ZR  " },
  { UNKNOWN_ATOM, "ERR "}
};

void 
PDBHetatmTypes::build_the_map()   
{
  for(uint i = 0; i < A_num_atom_strings; ++i)
    A_string_to_atom_type[A_atom_strings[i].name] = A_atom_strings[i].atom;
}
