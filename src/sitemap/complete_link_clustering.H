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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/complete_link_clustering.H,v $
 * $Revision: 1.5 $
 * $Author: vanvoor4 $
 * $Date: 2008-03-31 17:49:39 $
 *
 * $Log: not supported by cvs2svn $
 * Revision 1.4  2007/08/21 18:09:27  vanvoor4
 * Changed the point type from hbond_point to hphob_point
 *
 * Revision 1.3  2007/02/07 15:33:07  vanvoor4
 * point_t type renamed to hbond_point_t and resides in the SimSite3D
 * namespace.
 *
 * Revision 1.2  2006/10/20 13:20:00  vanvoor4
 * vectorified
 *
 * Revision 1.1  2006/03/30 18:39:14  vanvoor4
 * initial checkin
 *
 *
 */

#ifndef  _COMPLETE_LINK_CLUSTERING_H
#define  _COMPLETE_LINK_CLUSTERING_H

#include <HphobPoints.H> 

//! Use complete link clustering to reduce the number of hydrophobic points
void complete_link_clustering(SimSite3D::hphob_point_vec& hphob_pts, 
                              SimSite3D::hphob_point_list* points, 
			      double threshold);

#endif
