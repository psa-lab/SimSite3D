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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/hydro.H,v $
 * $Revision: 1.2 $
 * $Author: vanvoor4 $
 * $Date: 2007-08-21 15:32:21 $
 *
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2007/01/03 16:58:59  vanvoor4
 * Initial checkin.
 *
 * Revision 1.2  2006/08/25 17:29:18  vanvoor4
 * Changed format to be consistent with enum ordering.
 *
 * Revision 1.1  2006/03/30 18:42:22  vanvoor4
 * initial checkin
 *
 *
 */

// 20 + 1 and 37 + 1 is used here to show the enums have 0 reserved to denote
// an unknown value.  37 nonzero atom types are those that are found in the
// 20 residues and are used to determine hydrophobicity.  The rest--of the atom
// types--seem to be (at first glance) nobel gasses, metals, etc.
static int hydro[20+1][37+1] = {
// Placeholder
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
/* ALA : CA CB C O N */
{ 0, 72, 232, 7, 584, 355, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
/* ARG : CA CB C O N CD CG CZ NE NH1 NH2 */
{ 0, 39, 33, 12, 606, 470, 172, 0, 0, 0, 0, 0, 0, 85, 0, 0, 0, 26, 0, 0, 0, 0, 237, 0, 0, 440, 445, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
/* ASN : CA CB C O N CG ND2 OD1 */
{ 0, 21, 88, 17, 448, 367, 0, 0, 0, 0, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 0, 464, 0, 0, 0, 0, 0, 0, 459, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
/* ASP : CA CB C O N CG OD1 OD2 */
{ 0, 11, 82, 6, 503, 301, 0, 0, 0, 0, 0, 0, 0, 37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 546, 546, 0, 0, 0, 0, 0, 0, 0, 0 },
/* CYS : CA CB C O N SG */
{ 0, 36, 100, 14, 583, 400, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 139 },
/* GLN : CA CB C O N CD CG NE2 OE1 */
{ 0, 40, 64, 9, 502, 415, 20, 0, 0, 0, 0, 0, 0, 116, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 416, 0, 0, 0, 0, 0, 453, 0, 0, 0, 0, 0, 0, 0 },
/* GLU : CA CB C O N CD CG OE1 OE2 */
{ 0, 40, 70, 6, 513, 343, 14, 0, 0, 0, 0, 0, 0, 109, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 481, 469, 0, 0, 0, 0, 0, 0 },
/* GLY : CA C O N */
{ 0, 139, 0, 28, 521, 312, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
/* HIS : CA CB C O N CD2 CE1 CG ND1 NE2 */
{ 0, 23, 99, 22, 442, 359, 0, 0, 140, 0, 165, 0, 0, 60, 0, 0, 0, 0, 0, 0, 441, 0, 0, 0, 353, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
/* ILE : CA CB C O N CD1 CG1 CG2 */
{ 0, 24, 18, 0, 545, 380, 0, 87, 0, 0, 0, 0, 0, 0, 31, 124, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
/* LEU : CA CB C O N CD1 CD2 CG */
{ 0, 10, 35, 0, 548, 384, 0, 125, 117, 0, 0, 0, 0, 20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
/* LYS : CA CB C O N CD CE CG NZ */
{ 0, 12, 60, 9, 474, 271, 81, 0, 0, 129, 0, 0, 0, 76, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 476, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
/* MET : CA CB C O N CE CG SD */
{ 0, 30, 132, 30, 591, 486, 0, 0, 0, 208, 0, 0, 0, 95, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15, 0 },
/* PHE : CA CB C O N CD1 CD2 CE1 CE2 CG CZ */
{ 0, 14, 29, 0, 503, 406, 0, 67, 76, 0, 99, 48, 0, 0, 0, 0, 0, 81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
/* PRO : CA CB C O N CD CG */
{ 0, 49, 132, 35, 488, 16, 146, 0, 0, 0, 0, 0, 0, 164, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
/* SER : CA CB C O N OG */
{ 0, 22, 150, 0, 491, 245, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 568, 0, 0, 0, 0, 0 },
/* THR : CA CB C O N CG2 OG1 */
{ 0, 21, 54, 0, 601, 292, 0, 0, 0, 0, 0, 0, 0, 0, 0, 193, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 536, 0, 0, 0, 0 },
/* TRP : CA CB C O N CD1 CD2 CE2 CE3 CG CH2 CZ2 CZ3 NE1 */
{ 0, 46, 49, 31, 569, 429, 0, 149, 67, 0, 0, 100, 27, 0, 0, 0, 115, 0, 72, 71, 0, 0, 0, 349, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
/* TYR : CA CB C O N CD1 CD2 CE1 CE2 CG CZ OH */
{ 0, 6, 69, 74, 551, 374, 0, 83, 51, 0, 68, 50, 0, 22, 0, 0, 0, 36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 635, 0, 0, 0 },
/* VAL : CA CB C O N CG1 CG2 */
{ 0, 13, 0, 0, 534, 482, 0, 0, 0, 0, 0, 0, 0, 0, 85, 106, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
