dnl @synopsis TAC_ARG_CONFIG_CHACO
dnl
dnl Test a variety of Chaco options:
dnl --enable-chaco       - Turns Chaco on
dnl --with-chaco-libs    - specify Chaco libraries
dnl --with-chaco-libdir  - specify location of Chaco libraries
dnl
dnl If any of these options are set, HAVE_CHACO will be defined for both
dnl Autoconf and Automake, and HAVE_CHACO will be defined in the
dnl generated config.h file
dnl
dnl
dnl @author Kevin Long <krlong@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_CONFIG_CHACO],
[

AC_ARG_ENABLE(chaco,
[AC_HELP_STRING([--enable-chaco],[Chaco mesh partitioning])],
[HAVE_PKG_CHACO=$enableval],
[HAVE_PKG_CHACO=no]
)


AC_ARG_WITH(chaco-libs,
[AC_HELP_STRING([--with-chaco-libs="LIBS"],[Chaco libraries @<:@"-lchaco"@:>@])],
[
  CHACO_LIB=${withval}
  AC_MSG_CHECKING(user-defined Chaco libraries)
  AC_MSG_RESULT([${CHACO_LIB}])
]
)


AC_ARG_WITH(chaco-libdir,
[AC_HELP_STRING([--with-chaco-libdir=DIR],[Chaco library directory. Do not use -L])],
[
  CHACO_LIBDIR=${withval}
  AC_MSG_CHECKING(user-defined Chaco library directory)
  AC_MSG_RESULT([${CHACO_LIBDIR}])
]
)

AC_MSG_CHECKING(whether we are using Chaco)
AC_MSG_RESULT([${HAVE_PKG_CHACO}])

if test "X${HAVE_PKG_CHACO}" = "Xyes"; then
  AC_DEFINE(HAVE_CHACO,,[define if we want to use Chaco])
fi

dnl Define Automake version of HAVE_CHACO if appropriate

AM_CONDITIONAL(HAVE_CHACO, [test "X${HAVE_PKG_CHACO}" = "Xyes"])

#AC_SUBST(CHACO_LIB)   


])