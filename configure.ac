#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(splicegraphmerge, 0.1.0, shaomingfu@gmail.com)
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Check home directories for external libraries
AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost], [home directory for boost]), BOOST_HOME=$withval, BOOST_HOME=)

# Assign CXXFLAGS
AS_IF([test "x$BOOST_HOME" != "x"], [AC_SUBST([CXXFLAGS], ["-I$BOOST_HOME $CXXFLAGS"])], [])

# Checks for header files.
AC_CHECK_HEADERS([stdint.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT32_T

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
