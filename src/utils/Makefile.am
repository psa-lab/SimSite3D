INCLUDES = -I.
 
noinst_LIBRARIES = libASCbaseSoftwareUtils.a
noinst_PROGRAMS = compare_surfaces compute_dihedral_angles

libASCbaseSoftwareUtils_a_SOURCES = BaseParameters.C string_basics.C \
	stream_basics.C DiscreteSphere.C AtomTypes.C atom.C PDB_residues.C \
	PDB_metals.C PDBHetatmTypes.C Quaternion.C math_basics.C \
	param_tools.C my_float_array.C mat_ops.C Timer.C CoordFile.C PDBBase.C \
	PDBStructure.C mol2File.C BoundingVolume.C SimpleTrimesh.C \
	SimpleTrimeshTwo.C HbondGeometry.C prot_joint_dep.C \
        VertAttrib.C FaceAttrib.C SimpleTrimeshThree.C ImmovableTrimeshThree.C

compare_surfaces_SOURCES = compare_surfaces.C
compare_surfaces_LDADD = libASCbaseSoftwareUtils.a -lm $(POPT)

compute_dihedral_angles_SOURCES = compute_dihedral_angles.C
compute_dihedral_angles_LDADD = libASCbaseSoftwareUtils.a -lm $(POPT)
