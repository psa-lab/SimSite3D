
INCLUDES = -I../utils -I../sitemap -I.

noinst_LIBRARIES = libASCbaseSoftwareScore.a
bin_PROGRAMS = prot_lig_score 

SCORE_SRCS = ProtLigScoreParameters.C ProtLigScore.C ProtLigActTable.C 

libASCbaseSoftwareScore_a_SOURCES = $(SCORE_SRCS) 

prot_lig_score_SOURCES = $(SCORE_SRCS) prot_lig_score.C

prot_lig_score_LDADD = ../sitemap/libASCbaseSoftwareSitemaps.a \
	../utils/libASCbaseSoftwareUtils.a -lm $(POPT)
