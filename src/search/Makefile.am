
INCLUDES = -I../utils -I../sitemap -I../score -I.

noinst_LIBRARIES = libASCbaseSoftwareSearch.a
noinst_PROGRAMS = write_hbond_surfs_meshes new_tri_test test_grid_sampling \
	compute_dataset_matchprints
bin_PROGRAMS = search_sitemaps search_using_ArtSurf

ALIGN_SRCS = ScoreMapBase.C ModelSitemap.C MatchTriangles.C
SCORE_SRCS = WeightedSumsScore.C ExternalScoringFunction.C \
	DrugScoreInterface.C SFCscoreInterface.C ModelHbondSurfaces.C \
	point_and_surf_score.C TriMeshSphere.C HbondSurfacesScore.C \
	SurfDepsJoints.C IK_tests.C
LAPACK_LIBS = /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3

#TESTSCORING.C \
#PsuedoLigRmsd.C EllipsoidalHbondsScore.C SimplifiedMayoHbondsScore.C \
#ModelHbondSurfaces.C 

libASCbaseSoftwareSearch_a_SOURCES = $(ALIGN_SRCS) $(SCORE_SRCS) \
	SearchParameters.C

search_sitemaps_SOURCES = $(ALIGN_SRCS) $(SCORE_SRCS) \
	SearchParameters.C Search.C search_sitemaps.C

search_sitemaps_LDADD = ../score/libASCbaseSoftwareScore.a \
	../sitemap/libASCbaseSoftwareSitemaps.a \
	../utils/libASCbaseSoftwareUtils.a -lm $(POPT) $(LAPACK_LIBS)

#test_moving_stuff_SOURCES = $(ALIGN_SRCS) $(SCORE_SRCS) \
#	SearchParameters.C Search.C IK_tests.C test_moving_stuff.C

#test_moving_stuff_LDADD = ../score/libASCbaseSoftwareScore.a \
#	../sitemap/libASCbaseSoftwareSitemaps.a \
#	../utils/libASCbaseSoftwareUtils.a -lm $(POPT) \
#        /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3


write_hbond_surfs_meshes_SOURCES = $(ALIGN_SRCS) $(SCORE_SRCS) \
	SearchParameters.C write_hbond_surfs_meshes.C

write_hbond_surfs_meshes_LDADD = ../score/libASCbaseSoftwareScore.a \
	../sitemap/libASCbaseSoftwareSitemaps.a \
	../utils/libASCbaseSoftwareUtils.a -lm $(POPT) $(LAPACK_LIBS)


search_using_ArtSurf_SOURCES = $(ALIGN_SRCS) $(SCORE_SRCS) \
	SearchParameters.C Search.C search_using_ArtSurf.C

search_using_ArtSurf_LDADD = ../score/libASCbaseSoftwareScore.a \
	../sitemap/libASCbaseSoftwareSitemaps.a \
	../utils/libASCbaseSoftwareUtils.a -lm $(POPT) $(LAPACK_LIBS)

new_tri_test_SOURCES = $(ALIGN_SRCS) $(SCORE_SRCS) \
	SearchParameters.C Search.C new_tri_test.C

new_tri_test_LDADD = ../score/libASCbaseSoftwareScore.a \
	../sitemap/libASCbaseSoftwareSitemaps.a \
	../utils/libASCbaseSoftwareUtils.a -lm $(POPT) $(LAPACK_LIBS)

test_grid_sampling_SOURCES = $(ALIGN_SRCS) $(SCORE_SRCS) \
	SearchParameters.C Search.C test_grid_sampling.C

test_grid_sampling_LDADD = ../score/libASCbaseSoftwareScore.a \
	../sitemap/libASCbaseSoftwareSitemaps.a \
	../utils/libASCbaseSoftwareUtils.a -lm $(POPT) $(LAPACK_LIBS)

compute_dataset_matchprints_SOURCES = $(ALIGN_SRCS) $(SCORE_SRCS) \
	SearchParameters.C Search.C compute_dataset_matchprints.C

compute_dataset_matchprints_LDADD = ../score/libASCbaseSoftwareScore.a \
	../sitemap/libASCbaseSoftwareSitemaps.a \
	../utils/libASCbaseSoftwareUtils.a -lm $(POPT) $(LAPACK_LIBS)
