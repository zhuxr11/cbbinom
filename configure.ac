AC_INIT([cbbinom],[0.2.0],[zxr6@163.com])

dnl Find pkg-config program
PKG_PROG_PKG_CONFIG

dnl Select optional include paths, from a configure option
dnl or from an environment variable.
gmp_include_path=""
AC_ARG_WITH([gmp-include],
            AS_HELP_STRING([--with-gmp-include=GMP_INCLUDE],[specify the GMP include path]),
            [gmp_include_path=$withval])
CBBINOM_CPPFLAGS="-I."
if test [ -n "$gmp_include_path" ] ; then
   CBBINOM_CPPFLAGS="${CBBINOM_CPPFLAGS} -I${gmp_include_path}"
else
  if test [ -n "${CBBINOM_GMP_INCLUDE}" ] ; then
     CBBINOM_CPPFLAGS="${CBBINOM_CPPFLAGS} -I${CBBINOM_GMP_INCLUDE}"
  fi
fi

mpfr_include_path=""
AC_ARG_WITH([mpfr-include],
            AS_HELP_STRING([--with-mpfr-include=MPFR_INCLUDE],[specify the MPFR include path]),
            [mpfr_include_path=$withval])
if test [ -n "$mpfr_include_path" ] ; then
   CBBINOM_CPPFLAGS="${CBBINOM_CPPFLAGS} -I${mpfr_include_path}"
else
  if test [ -n "${CBBINOM_MPFR_INCLUDE}" ] ; then
     CBBINOM_CPPFLAGS="${CBBINOM_CPPFLAGS} -I${CBBINOM_MPFR_INCLUDE}"
  else
    if test -n "${PKG_CONFIG}"; then
      MPFR_CFLAGS_PKG_CONFIG=`"${PKG_CONFIG}" --cflags mpfr`
      CBBINOM_CPPFLAGS="${CBBINOM_CPPFLAGS} ${MPFR_CFLAGS_PKG_CONFIG}"
    fi
  fi
fi

dnl Select optional library paths
AC_ARG_WITH([gmp-lib],
            AS_HELP_STRING([--with-gmp-lib=GMP_LIB],[the location of GMP library]),
            [gmp_lib_path=$withval])
if test [ -n "gmp_lib_path" ] ; then
   LIBS="-L${gmp_lib_path} -lgmp ${LIBS}"
else
  if test [ -n "${CBBINOM_GMP_LIBS}" ] ; then
     LIBS="-L${CBBINOM_GMP_LIBS} -lgmp ${LIBS}"
  fi
fi

AC_ARG_WITH([mpfr-lib],
            AS_HELP_STRING([--with-mpfr-lib=MPFR_LIB],[the location of MPFR library]),
            [mpfr_lib_path=$withval])
if test [ -n "mpfr_lib_path" ] ; then
   LIBS="-L${mpfr_lib_path} -lmpfr ${LIBS}"
else
  if test [ -n "${CBBINOM_MPFR_LIBS}" ] ; then
     LIBS="-L${CBBINOM_MPFR_LIBS} -lmpfr ${LIBS}"
  else
    if test -n "${PKG_CONFIG}"; then
      MPFR_LIBS_PKG_CONFIG=`"${PKG_CONFIG}" --libs mpfr`
      LIBS="${MPFR_LIBS_PKG_CONFIG} ${LIBS}"
    fi
  fi
fi

dnl Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

CPPFLAGS="${CPPFLAGS} ${CBBINOM_CPPFLAGS}"

dnl Check the headers can be found
AC_CHECK_HEADERS(mpfr.h gmp.h)
if test "${ac_cv_header_gmp_h}" = no; then
  AC_MSG_FAILURE([
  ---------------------------------------------
   'gmp' and its header files are required.

   Please install:

   * gmp (Linux/MacOS)
   * Rtools40 or above, which includes MPFR and GMP (Windows)

   and try again.

   If you believe this library is installed on your system but
   this script is simply unable to find it, you can specify the
   include and lib paths manually:

   R CMD INSTALL ${PACKAGE_NAME} \\
     --configure-vars='LIBS=-L/path/to/libs CPPFLAGS=-I/path/to/headers'
  ---------------------------------------------])
fi
if test "${ac_cv_header_mpfr_h}" = no; then
  AC_MSG_FAILURE([
  ---------------------------------------------
   'mpfr' and its header files are required.

   Please install:

   * mpfr (Linux/MacOS)
   * Rtools40 or above, which includes MPFR and GMP (Windows)

   and try again.

   If you believe this library is installed on your system but
   this script is simply unable to find it, you can specify the
   include and lib paths manually:

   R CMD INSTALL ${PACKAGE_NAME} \\
     --configure-vars='LIBS=-L/path/to/libs CPPFLAGS=-I/path/to/headers'
  ---------------------------------------------])
fi

dnl substitute CBBINOM_CPPFLAGS and LIBS
AC_SUBST(CBBINOM_CPPFLAGS)
AC_SUBST(LIBS)
AC_CONFIG_HEADERS([src/config.h])
dnl and do subsitution in the src/Makevars.in and src/config.h
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

echo "
  --------------------------------------------------
  Configuration for ${PACKAGE_NAME} ${PACKAGE_VERSION}

    cppflags: ${CPPFLAGS} ${CBBINOM_CPPFLAGS}
    libs:     ${LIBS}

  --------------------------------------------------
"
