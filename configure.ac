AC_INIT(SimSite3D, 4.6-rc1)

AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_MAKE_SET
# Use libtool for generating static and/or dynamic libs instead of ranlib
AC_PROG_LIBTOOL
# For make rules to auto update libtool conf files
AC_SUBST(LIBTOOLS_DEPS)
AC_ISC_POSIX
AC_HEADER_STDC
AC_PROG_CXX

AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)

AC_CHECK_LIB([popt], [poptGetContext],
	[POPT="-lpopt"], AC_MSG_WARN(Did not find a popt library))
AC_SUBST(POPT)

# These next lines handle the configuration and execution of documentation

AC_ARG_ENABLE(refman,
[  --enable-refman         generate reference manual [default=yes]],
[case "${enableval}" in
  yes) usedoxy=true ;;
  no)  usedoxy=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-doxygen) ;;
esac],[usedoxy=true])

# depends on having doxygen around
if test x${usedoxy} = xtrue; then
  AC_DOXYGEN
  if test x${have_doxygen} = xno ; then
    usedoxy = false
  fi
fi

AM_CONDITIONAL(DOXY, test x$usedoxy = xtrue)

echo "Build target is $target"
case "$target" in
  sparc-sun-solaris*)
    if test x${GCC} = xyes; then
      CFLAGS="-g -O -Wall -Wundef -Wshadow -Wmissing-declarations";
    else
      CFLAGS="-g -xCC";
    fi;;
# Next line for debugging since gdb doesn't like optimized code.
#  *) CFLAGS="-g -Wall -Wundef -Wshadow"
#  *) CFLAGS="-O2 -pg"
#  *) CFLAGS=" -g -fPIC -Wall -Wundef -Wshadow"
#  *) CFLAGS="-O2 -g -fPIC -Wall -Wundef -Wshadow -fprofile-generate"
  *) CFLAGS="-O2 -g -fPIC -Wall -Wundef -Wshadow"
#  *) CFLAGS="-O2 -fPIC"
esac
CXXFLAGS="$CFLAGS"

AC_CONFIG_FILES([
	Makefile
	src/Makefile
	src/utils/Makefile
	src/sitemap/Makefile
        src/score/Makefile
	src/search/Makefile
	src/doc/Makefile
	src/doc/Doxyfile
	])
AC_OUTPUT
