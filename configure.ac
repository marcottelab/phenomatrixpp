#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.64])
AC_LANG([C++])

AC_INIT([phenomatrix], [0.1.0], [john@gunfreeschools.org])
AC_CONFIG_AUX_DIR([/usr/share/automake-1.7])
AC_CONFIG_SRCDIR([src/oracle.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
have_boost_system=no
have_boost_program_options=no
have_boost_filesystem=no
have_boost=no

AX_BOOST_BASE([1.39.0])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_PROGRAM_OPTIONS

#AC_ARG_WITH(boost,
#    [ --with-boost=<path> prefix of boost installation. e.g. /usr/local or /home/user/usr/local],
#    [BOOST_PREFIX=$with_boost],
#    AC_MSG_ERROR([You must call configure with the --with-boost option.
#    This tells configure where to find the Boost library and headers.
#    e.g. --with-boost=/usr/local or --with-boost=/home/user/usr/local])
#)

AC_SEARCH_LIBS([get_system_category],[boost_system],[have_boost_system=yes])

# Checks for header files.
if test "x${have_boost_system}" = xyes; then
  AC_CHECK_HEADERS([boost_system.h], [], [have_boost=no])
fi

#if test "x${have_boost}" = xno; then
#  AC_MSG_ERROR([Could not find some or all of the boost libraries. These are required for building.])
#  (exit 1); exit 1;
#fi

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AM_INIT_AUTOMAKE

AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/version.h
])
AC_OUTPUT

