# $Source: /psa/share/repository/pfizer_proj/acinclude.m4,v $
# $Revision: 1.1 $
# $Author: vanvoor4 $
# $Date: 2007-02-07 16:08:35 $
#
# $Log: not supported by cvs2svn $
#
#
# The format of this file is similar to that used by Dr. Karlis Kaugars of
# Western Michigan University.

dnl
dnl AC_DOXYGEN
dnl
AC_DEFUN([AC_DOXYGEN],[
  AC_ARG_WITH(doxygen, [  --with-doxygen=PATH           path to doxygen executable],
                         doxygen="$withval",
                         doxygen="/usr/bin/doxygen")

  AC_MSG_CHECKING(for doxygen)
  if test -x ${doxygen} ; then
    have_doxygen=yes
    DOXYGEN="${doxygen}"
    AC_MSG_RESULT(yes)
  else
    have_doxygen=no
    DOXYGEN=""
    AC_MSG_RESULT(no)
  fi
  AC_SUBST(DOXYGEN)
])
