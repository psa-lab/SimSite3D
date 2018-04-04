.. _user_guide:

************************************************
Notes for Using SimSite3D and its Python Scripts
************************************************

.. note:: On the PSA systems, SimSite3D is currently installed from
          SimSite3D_hg, and it has a virtualenv for the python settings.

Using PSA's virtualenv
======================

To use the Python utilities with SimSite3D, I really recommend using the
virtualenv.

To use Python 2.7.2, SimSite3D's python scripts, SimSite3D's Boost::Python
extentions; **and** to setup your paths for SimSite3D::

  source /soft/linux64/.virtualenvs/SimSite3D_4.5/bin/activate.csh

When you need to use other tools or other tools fail, you can deactivate
the virtualenv using::

  deactivate

.. _dataset_setup:

Setting up a Dataset
====================

.. note:: Many of the conventions in this document are not required for the core
          SimSite3D functionality, but are generally required for the Python
          scripts.  Also, if you don't stick to a file naming convention
          (either this one or your own) you may find it difficult to reproduce
          your work or recall what you did when you come back to your results
          at a later date.

SimSite3D generally requires a structure file for each binding site.
Each site may have a ligand file, and will have SimSite3D generated site map
files.  It is recommended to have a specific directory structure and file
naming convention for SimSite3D screening and analysis.

.. _structure_name_template:

Binding Site Naming Convention
------------------------------

A standard naming convention for binding sites is useful for yourself and
others to know which sites you are using in your searches.  You should use
the following template to name sites::

  XXXX_YYY_T.EXT

Where we replace ``XXXX`` with the name of the structure, ``YYY`` with the name
of the ligand, ``T`` with the type of object, and ``EXT`` for the file 
extention.  As an example, suppose we have the PDB structure 1u72 and are 
interested in the methotrexate ligand (MTX).  Then we would have the following::

  1u72_MTX_p.pdb   # PDB structure file
  1u72_MTX_l.mol2  # ligand mol2 file

SimSite3D itself is expected to handle any non-zero length name for both
the protein and ligand, but some of the python scripts still expect
a 4 character protein name and a 3 character ligand name.
Also, do not skip any of the names or '_' characters or most programs will
fail.

.. _directory_structure:

One Directory Structure for SimSite3D Datasets
----------------------------------------------

Let the '~' character denote the root directory where you wish to work with,
manipulate, and save your SimSite3D datasets.  It is good practice to name
the dataset well, and to use the following commands to create the 
proteins, ligands, and dbase directories for that dataset::

  mkdir ~/my_dataset

  mkdir ~/my_dataset/proteins
  mkdir ~/my_dataset/ligands
  mkdir ~/my_dataset/dbase

If you are validating SimSite3D, I recommend aligning the proteins and ligands
before generating the site maps.  Depending on the type of dataset you can use
Dali, SSM, ligand based alignments, or a hybrid approach.

.. _column_order_file:

Column Order File
-----------------

We also need to create a file in the dataset directory called column_order.txt.
This file should contain lines like the following::

  site_id_1|pdb_code_1|plot_name_1|
  site_id_2|pdb_code_2|plot_name_2|

Where ``site_id_i`` is replaced by the name of the site, ``pdb_code_i`` is 
replaced by the pdb code for that structure, and ``plot_name_i`` is replace by 
the name that you would like to use to represent that structure on plots.  
As an example, for a pterin binding site dataset one might have::

  1u72_mtx|1u72|Hs DHFR|
  2fzh_dh1|2fzh|Pc DHFR|
  1qzf_fol|1qzf|Ch DHFR|
  1df7_mtx|1df7|Mt DHFR|
  1j3i_wra|1j3i|Pf DHFR|
  1dr1_hbi|1dr1|Gg DHFR|
  1aoe_gw3|1aoe|Ca DHFR|
  1d1g_mtx|1d1g|Tm DHFR|
  2qx0_ph2|2qx0|Yp HPPK|
  2bmb_pmm2|2bmb|Sc HPPK|
  1q0n_ph2|1q0n|Ec HPPK(t)|
  1rb0_hh2|1rb0|Ec HPPK(b)|
  1cbk_roi|1cbk|Hi HPPK|
  1mmk_h4b|1mmk|Hs PAH|
  1mlw_hbi|1mlw|Hs TPH|
  2toh_hbi|2toh|Rn TH|
  1ltz_hbi|1ltz|Cv PAH|
  2bmb_pmm1|2bmb|Sc DHPS|

.. _simsite3d_scoring_notes:

Notes on SimSite3D Scoring
==========================

The SimSite3D tool **search_sitemaps** writes the search results into a text
file, hereafter referred to as "a/the results file".  
The convention is any line that starts with the '#' character is a comment 
line.  
To help with reproducibility, a number of SimSite3D parameters are written to 
the "comment header" of the results file.
The score lines (if any) will follow the comment header.

SimSite3D Results File Comment Header
-------------------------------------

Shown below is an example of a header from a recent SimSite3D results file.  
We will go through line by line to repeat what is shown.

.. code-block:: python
  :linenos:

  # Parameters for search_sitemaps (SimSite3D) 4.5-rc4 run
  # Local start time:                               Thu Dec 29 11:29:00 2011
  # Working directory (via getcwd()):               /Users/jvanvoorst/data/SimSite3D/results/baseline
  # Searchable sitemaps directory:                  /Users/jvanvoorst/data/SimSite3D/datasets/test/adenines/dbase
  # Corresponding ligands directory:                /Users/jvanvoorst/data/SimSite3D/datasets/test/adenines/ligands
  # Query (Model) Sitemap Name:                     /Users/jvanvoorst/data/SimSite3D/datasets/test/adenines/ade_pockets/1b38_atp_s.csv
  # Max number of scores to keep for each sitemap:  1
  # Score threshold:                                10
  # Minimum number of atoms required in a fragment: 0
  # Average distance metric error tolerance:        0.3
  # Average least squares error tolerance:          0.3
  # Highly hydrophobic query pocket:                Yes
  # SimSite3D timing statistics:
  #   Wall clock time:                              11.64 sec.
  #   CPU time:                                     11.63 sec.
  #   User time:                                    11.58 sec.
  #   Kernel time:                                  0.05 sec.
  # Scoring function terms were scaled:             No
  # Max dist between corresponding surface points:  1.5
  # SimSite3D alignments scores are not normalized
  #
  # Fields:
  # 1 ) Name of ligand fragment corresponding to the score record (line)
  # 2 ) Raw SimSite3D alignment score of target to query
  # 3 ) Rotation matrix to align target to query
  # 4 ) Translation vector to move target to query
  # 5 ) Match print of the query's sitemap points satisfied by sitemap points
  #     in the database hit
  # 6 ) Ligand fragment binary string:  1 or 0 in nth position implies that the
  #     nth mol2 ligand atom is or is not in the mol2 ligand fragment (resp.)
  # 7 ) scoring function terms/alignment features

What questions can I answer with the information can be found in this header?

  #. Which version of SimSite3D (search_sitemaps) was used to search?
  #. When was this search started? (this will depend on the operating system's
     locale)
  #. From which directory did I run the search?
  #. What was the dataset (screening) site maps directory?
  #. In which directory did I instruct search_sitemaps to search for dataset 
     ligand files?
  #. Which file was used to specify the query site map?
  #. What was the maximum number of scores to keep per dataset site map?
  #. What was the score threshold (any alignments that scored numerically 
     higher than this were ignored)?
  #. What was the minimum number of atoms required for dataset ligand 
     fragments?

    * If this value is/was zero, ligands were ignored for the purpose of
      selecting site alignments
    * If this value is/was greater than zero, any listed site alignments will 
      have at least one dataset ligand fragment (with at least that many atoms)
      fully contained (after site alignment) within the query's site map 
      volume

  10. What was the distance matrix error (DME) threshold?
  #. What was the weighted RMSD threshold for triangle matches?
  #. Was the highly hydrophobic ligand flag on/off for search_sitemaps?
  #. Search timing using itimers (note that the timing could be off somewhat,
     and the precision is probably best kept to 1/10ths of seconds).
  #. How much "actual" time elapsed between search_sitemaps starting and
     finishing?
  #. How long was search_sitemaps executing on one CPU?
  #. How much of that time was spent in user space?
  #. How much of the CPU time was spent in the Linux kernel?
  #. Was a scoring function used that required scaled features?
  #. What was the threshold used to determine if two molecular surface points
     (one from query surface and one from dataset surface) could be 
     corresponding points?
  #. Were the scores normalized?
  #. (blank)
  #. What are the fields (columns) listed for hits/matches?
  #. What was the name of the dataset site map?
  #. What was the score for this match?
  #. What rotation matrix may I use to align the dataset site (protein, ligand,
     and/or molecular surface) to the query?
  #. What translation vector corresponds with the rotation given above?
  #. Which points in the query site map had a corresponding point in the
     dataset site?
  #. ...
  #. Which non-hydrogen atoms in the dataset ligand were inside the query
     site map's volume?
  #. ...
  #. What was the scoring function feature vector for this match?

Sundry Comments/Remarks About SimSite3D Scoring
-----------------------------------------------

.. note:: You may skip this section.  It is a number of comments about the 
          scoring process that are primarily of interest to those who
          are developing the tools.

Automatic score normalization is not handled robustly with respect to 
altering of experimental parameters.  
The algorithm and rules were based on the improper assumption that SimSite3D
runs would be mostly static and SimSite3D parameters would be held constant.
This means that to remove all doubt, you should consider running 
search_sitemaps with the --no_normalization flag **and** run search_sitemaps
against a normalization dataset of your choice with **exactly** the same
parameters as the search *except* for the screening dataset and corresponding
ligands directories.

The crux of the issue is the normalization stats are written into the 
query sitemap file, and once there they are not computed again.  
Also, SimSite3D now has too many parameters to reliably and reasonably list 
them with a particular set of normalization stats.  
We can still keep this step automated, but, in my opinion, we should recompute 
the normalization stats for each search (rather than writing them into the 
query's sitemap files).
In fact, although it puts more burden on users, it is much safer to 
compute the normalization as an explicit and separate step (however, adding
additional steps will likely chase away most users).

Finally, score normalization was deemed a very important part of SimSite3D.
Thus, score normalization is enabled by default at site map generation time
and at search time.  
For now, if you want to do things your own way, you will need to
use the --no_normalization flag (at least when searching both your 
dataset and your normalization dataset, and to save time, when you 
create your dataset site maps).


.. _setup_simsite_and_plotting:

Howto use SimSite3D and Plotting Tools for Validation and Testing
==================================================================

.. note:: If you have not setup your dataset directory(ies) please do so first
          before continuing with this section.

Directory Names
---------------

As a generally good practice, I recommend creating a results directory for
SimSite3D and, if you prefer, a subdirectory for your specific datasets.  
One example of this is::

  mkdir ~/SimSite3D_results  # main results directory
  mkdir ~/SimSite3D_results/test_some_thing # my current results directory
  mkdir ~/SimSite3D_results/test_some_thing/pterins # my current pterins results directory if I want results segregated by dataset

SimSite3D Searches
------------------

Run SimSite3D to create a results file that starts with the same prefix as the
query site.  For example if I am using 2toh_hbi as the query site I would use::

  search_sitemaps --proj_output /path/to/my/results/2toh_hbi_results

If we run such a command for each of the query sites in the pterins dataset,
we can plot an NxN matrix that shows the relative score of each query 
versus each dataset site.  The plotting is quite easy if you followed 
a consistent naming convention and can be miserable if you have files
all over the place.

Computing RMSD of Alignment
---------------------------

.. note:: To compute RMSD of alignment you must have aligned the sites before
          doing the SimSite3D searches.  I recommend doing this by having 
          aligned structures and ligands, and then creating the site maps.
          Following such a pattern generally results in less work and issues
          of "was this aligned or not".

The general ideas and methods of "computing RMSD of site alignment" as was done 
for Jeff's dissertation is presented here (of course, you are not required
to use this method, but if you want to reproduce similar results, it is 
likely necessary). 
You will need to figure out how to best align your sites, and follow the 
following steps:

  #. Align the sites/proteins/ligands using whatever tool you like to some 
     reference 
     (Jeff used one "reference" structure for each dataset and aligned all
     others to that structure)
  #. Save the transformations -- you will wish you had later if you don't
  #. Align the proteins (and ligands if there are any used to generate sites)
  #. Generate the site maps
  #. Do a search
  #. Use a program to transform the pocket for each dataset site using the
     saved transformation (in the results .out) file
  #. Compute the RMSD of alignment as the RMSD of the transformed pocket with
     the coordinates its initial alignment to the reference




