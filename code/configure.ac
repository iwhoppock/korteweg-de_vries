# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([3_KdV_Hoppock], [0.01])
AM_INIT_AUTOMAKE([-Wall foreign])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])


# Checks for programs.

AC_PROG_CC_C99
if test -z "$CC" || test "$CC" = "no"; then
  AC_MSG_ERROR([No suitable C99 compiler found])
fi
AM_PROG_CC_C_O

# add -Wall if possible
AX_CHECK_COMPILER_FLAGS([-Wall], [CFLAGS="$CFLAGS -Wall"])

# find flags for using OpenMP
AC_OPENMP
CFLAGS="$CFLAGS $OPENMP_CFLAGS"

# link with the math library if needed
AC_CHECK_LIB([m], [sin])


# for library support
AC_PROG_RANLIB
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

echo ""
echo "                       Configure Summary"
echo "                       ================="
echo ""
echo "Compiling Options:"
echo "------------------"
echo "                     C Compiler: $CC"
echo "                       CPPFLAGS: $CPPFLAGS"
echo "                         CFLAGS: $CFLAGS"
echo "                        LDFLAGS: $LDFLAGS"
                                                                   