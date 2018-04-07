/******************************************************************************
 * search_sitemaps is part of the SimSite3D Software.
 * SimSite3D Software version 1.00, Copyright (C) 2006, Michigan State University
 *      (MSU) Board of Trustees.
 * SimSite3D Software comes with ABSOLUTELY NO WARRANTY.
 *
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

/*
 * $Source: /psa/share/repository/pfizer_proj/src/search/align_sitemaps.C,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2007-01-24 16:18:10 $
 * 
 * $Log: not supported by cvs2svn $
 *
 * 
 * 
 */

#include <SimSite3DMain.H>

using namespace SimSite3D;

int main(const int argc, const char **argv)
{
  SimSite3DMain my_main(argc, argv);
  if(!my_main.fail()) my_main.get_alignments_only();
}
