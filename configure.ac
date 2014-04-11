AC_INIT(radx, version-1.0.0)

# set language to cpp
AC_PROG_CXX
AC_LANG(C++)

# check for MacOS and append the architecture flag
if test `uname` == "Darwin"
then
  CXXFLAGS="$CXXFLAGS -arch x86_64"
  echo "checking whether we're compiling on Darwin/MacOS... yes"
else
  echo "checking whether we're compiling on Darwin/MacOS... no"
fi

# check for gsl
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
ac_cv_lib_gsl_matrix=no
SAVED_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS -lgsl"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include "gsl/gsl_matrix_double.h"],
    [gsl_matrix_calloc(5,5)])],
  [ac_cv_lib_gsl_matrix=yes],
  [])
LDFLAGS=$SAVED_LDFLAGS
if test "$ac_cv_lib_gsl_matrix" != yes
then 
  echo "checking for linking against gsl matrix... no"
else
  echo "checking for linking against gsl matrix... yes"
fi

GSLERROR="could not find or link to GSL -- radX requires GSL and \
cannot be built without it. GSL can be downloaded from \
http://www.gnu.org/software/gsl/ -- please ensure GSL is \
correctly installed on your system."
if test "$ac_cv_lib_m_cos" != yes 
then
  AC_MSG_ERROR([$GSLERROR])
fi
if test "$ac_cv_lib_gsl_gsl_blas_dgemm" != yes
then
  AC_MSG_ERROR([$GSLERROR])
fi 
if test "$ac_cv_lib_gslcblas_cblas_dgemm" != yes
then
  AC_MSG_ERROR([$GSLERROR])
fi
if test "$ac_cv_lib_gsl_matrix" != yes
then
  AC_MSG_ERROR([$GSLERROR])
fi

AC_OUTPUT(Makefile)
AC_OUTPUT(src/Makefile)
#AC_OUTPUT(src/test/Makefile)
AC_OUTPUT(src/progs/Makefile)
AC_OUTPUT(src/common/Makefile)
AC_OUTPUT(src/common/config.hpp)
#AC_OUTPUT(src/test/regressionTests.xml)
