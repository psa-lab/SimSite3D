
INCLUDES = -I../utils -I../sitemap -I.

noinst_LIBRARIES = libSimSite3DSoftwareScore.a
bin_PROGRAMS = prot_lig_score 

SCORE_SRCS = ProtLigScoreParameters.C ProtLigScore.C ProtLigActTable.C 

libSimSite3DSoftwareScore_a_SOURCES = $(SCORE_SRCS) 

prot_lig_score_SOURCES = $(SCORE_SRCS) prot_lig_score.C

prot_lig_score_LDADD = ../sitemap/libSimSite3DSoftwareSitemaps.a \
	../utils/libSimSite3DSoftwareUtils.a -lm $(POPT)
