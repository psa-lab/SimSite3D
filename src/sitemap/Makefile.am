#
# $Source: /psa/share/repository/pfizer_proj/src/gen_points/Makefile.am,v $
# $Revision: 1.15 $
# $Author: vanvoor4 $
# $Date: 2008-05-14 14:25:25 $
#
# $Log: not supported by cvs2svn $
# Revision 1.14  2007/12/17 21:25:31  vanvoor4
# Added a few more source files.
#
# Revision 1.13  2007/11/01 16:25:09  vanvoor4
# Added the comp* binaries to noinst_PROGRAMS
#
# Revision 1.12  2007/11/01 16:20:38  vanvoor4
# Updated to reflect GenPoints -> GenPointsParameters and gen_points.C
# The comp_* binaries should not be installed
#
# Revision 1.11  2007/08/29 20:27:20  vanvoor4
# added a binary to compute simple rmsd
#
# Revision 1.10  2007/08/21 18:34:26  vanvoor4
# Added a library for the Sitemap class and supporting classes so that it
# can also be used by the search routines.
#
# Revision 1.9  2007/03/23 16:42:27  vanvoor4
# Moved PDBFile.* to ../basics
#
# Revision 1.8  2007/02/07 15:26:34  vanvoor4
# Added some source files, the popt lib, and the lib from ../basics.
#
# Revision 1.7  2006/10/20 13:17:19  vanvoor4
# A bit more classified stuff.
#
# Revision 1.6  2006/08/25 17:06:12  vanvoor4
# Redesign of template code.
#
# Revision 1.5  2006/08/25 15:15:14  vanvoor4
# Redesign of template code.
#
#
# Revision 1.4  2006/04/14 20:23:46  vanvoor4
# Changed from template to gen_points.
#
# Revision 1.3  2006/04/11 20:25:13  vanvoor4
# Did not change ligand based generation as yet.
#
# Revision 1.2  2006/04/04 17:36:20  vanvoor4
# Added quicksort.C
#
# Revision 1.1  2006/03/30 18:30:59  vanvoor4
# initial checkin
#
#

INCLUDES = -I. -I../utils

noinst_LIBRARIES = libSimSite3DSoftwareSitemaps.a
noinst_PROGRAMS = comp_simple_rmsd comp_simple_site_rmsd comp_site_rmsd \
	norm_db_ade_site_rmsd
bin_PROGRAMS = gen_points 

libSimSite3DSoftwareSitemaps_a_SOURCES =  GenPointsParameters.C WaterPoints.C \
        HbondPoints.C HphobPoints.C complete_link_clustering.C \
	SitemapPointsFile.C SitemapAtomsFile.C normalize_sitemap.C Sitemap.C \
	HbondBase.C HbondSurfaces.C

gen_points_SOURCES = gen_points.C
gen_points_LDADD = ./libSimSite3DSoftwareSitemaps.a \
	../utils/libSimSite3DSoftwareUtils.a -lm $(POPT)

comp_simple_rmsd_SOURCES = comp_simple_rmsd.C
comp_simple_rmsd_LDADD = ./libSimSite3DSoftwareSitemaps.a \
	../utils/libSimSite3DSoftwareUtils.a -lm

comp_simple_site_rmsd_SOURCES = comp_simple_site_rmsd.C
comp_simple_site_rmsd_LDADD = ./libSimSite3DSoftwareSitemaps.a \
	../utils/libSimSite3DSoftwareUtils.a -lm

comp_site_rmsd_SOURCES = comp_site_rmsd.C
comp_site_rmsd_LDADD = ./libSimSite3DSoftwareSitemaps.a \
	../utils/libSimSite3DSoftwareUtils.a -lm

norm_db_ade_site_rmsd_SOURCES = norm_db_ade_site_rmsd.C
norm_db_ade_site_rmsd_LDADD = ./libSimSite3DSoftwareSitemaps.a \
	../utils/libSimSite3DSoftwareUtils.a -lm
