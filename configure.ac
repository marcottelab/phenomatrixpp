#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.61])

#AC_LANG(C++)

AC_INIT([phenomatrix], [0.1.0], [john.woods@marcottelab.org])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([/usr/share/automake-1.7])
AC_CONFIG_SRCDIR([src/oracle.h])
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_HOST
AC_DEFINE_UNQUOTED([CHOST], ["$host"], [Canonical host])
AC_CANONICAL_BUILD
AC_DEFINE_UNQUOTED([CBUILD], ["$build"], [Canonical build])
AC_CANONICAL_TARGET
AC_DEFINE_UNQUOTED([CTARGET], ["$target"], [Canonical target])

# Checks for programs.
AC_PROG_CXX
AC_PROG_LIBTOOL

# Checks for libraries.
#AX_BOOST_BASE([1.41.0])
#AX_BOOST_SYSTEM
#AX_BOOST_FILESYSTEM
#AX_BOOST_PROGRAM_OPTIONS
BOOST_REQUIRE([1.41.0])
BOOST_FILESYSTEM
BOOST_SYSTEM
BOOST_PROGRAM_OPTIONS
# AC_CHECK_HEADERS(boost/filesystem.hpp, /usr/local/include)
#
#have_boost_system=no
#AC_CHECK_FUNC(get_system_category, AC_DEFINE(HAVE_BOOST_SYSTEM),)
#AC_SEARCH_LIBS([get_system_category],[boost_system],[have_boost_system=yes],[have_boost_system=no])

# Checks for header files.
#AC_CHECK_HEADERS([boost_system.h], [], [have_boost=no])

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
])
AC_OUTPUT

