# aclocal.m4 generated automatically by aclocal 1.6.3 -*- Autoconf -*-

# Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002
# Free Software Foundation, Inc.
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

# Like AC_CONFIG_HEADER, but automatically create stamp file. -*- Autoconf -*-

# Copyright 1996, 1997, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

AC_PREREQ([2.52])

# serial 6

# When config.status generates a header, we must update the stamp-h file.
# This file resides in the same directory as the config header
# that is generated.  We must strip everything past the first ":",
# and everything past the last "/".

# _AM_DIRNAME(PATH)
# -----------------
# Like AS_DIRNAME, only do it during macro expansion
AC_DEFUN([_AM_DIRNAME],
       [m4_if(regexp([$1], [^.*[^/]//*[^/][^/]*/*$]), -1,
	      m4_if(regexp([$1], [^//\([^/]\|$\)]), -1,
		    m4_if(regexp([$1], [^/.*]), -1,
			  [.],
			  patsubst([$1], [^\(/\).*], [\1])),
		    patsubst([$1], [^\(//\)\([^/].*\|$\)], [\1])),
	      patsubst([$1], [^\(.*[^/]\)//*[^/][^/]*/*$], [\1]))[]dnl
])# _AM_DIRNAME


# The stamp files are numbered to have different names.
# We could number them on a directory basis, but that's additional
# complications, let's have a unique counter.
m4_define([_AM_STAMP_Count], [0])


# _AM_STAMP(HEADER)
# -----------------
# The name of the stamp file for HEADER.
AC_DEFUN([_AM_STAMP],
[m4_define([_AM_STAMP_Count], m4_incr(_AM_STAMP_Count))dnl
AS_ESCAPE(_AM_DIRNAME(patsubst([$1],
                               [:.*])))/stamp-h[]_AM_STAMP_Count])


# _AM_CONFIG_HEADER(HEADER[:SOURCES], COMMANDS, INIT-COMMANDS)
# ------------------------------------------------------------
# We used to try to get a real timestamp in stamp-h.  But the fear is that
# that will cause unnecessary cvs conflicts.
AC_DEFUN([_AM_CONFIG_HEADER],
[# Add the stamp file to the list of files AC keeps track of,
# along with our hook.
AC_CONFIG_HEADERS([$1],
                  [# update the timestamp
echo 'timestamp for $1' >"_AM_STAMP([$1])"
$2],
                  [$3])
])# _AM_CONFIG_HEADER


# AM_CONFIG_HEADER(HEADER[:SOURCES]..., COMMANDS, INIT-COMMANDS)
# --------------------------------------------------------------
AC_DEFUN([AM_CONFIG_HEADER],
[AC_FOREACH([_AM_File], [$1], [_AM_CONFIG_HEADER(_AM_File, [$2], [$3])])
])# AM_CONFIG_HEADER

# Add --enable-maintainer-mode option to configure.
# From Jim Meyering

# Copyright 1996, 1998, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 1

AC_DEFUN([AM_MAINTAINER_MODE],
[AC_MSG_CHECKING([whether to enable maintainer-specific portions of Makefiles])
  dnl maintainer-mode is disabled by default
  AC_ARG_ENABLE(maintainer-mode,
[  --enable-maintainer-mode enable make rules and dependencies not useful
                          (and sometimes confusing) to the casual installer],
      USE_MAINTAINER_MODE=$enableval,
      USE_MAINTAINER_MODE=no)
  AC_MSG_RESULT([$USE_MAINTAINER_MODE])
  AM_CONDITIONAL(MAINTAINER_MODE, [test $USE_MAINTAINER_MODE = yes])
  MAINT=$MAINTAINER_MODE_TRUE
  AC_SUBST(MAINT)dnl
]
)

# AM_CONDITIONAL                                              -*- Autoconf -*-

# Copyright 1997, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 5

AC_PREREQ(2.52)

# AM_CONDITIONAL(NAME, SHELL-CONDITION)
# -------------------------------------
# Define a conditional.
AC_DEFUN([AM_CONDITIONAL],
[ifelse([$1], [TRUE],  [AC_FATAL([$0: invalid condition: $1])],
        [$1], [FALSE], [AC_FATAL([$0: invalid condition: $1])])dnl
AC_SUBST([$1_TRUE])
AC_SUBST([$1_FALSE])
if $2; then
  $1_TRUE=
  $1_FALSE='#'
else
  $1_TRUE='#'
  $1_FALSE=
fi
AC_CONFIG_COMMANDS_PRE(
[if test -z "${$1_TRUE}" && test -z "${$1_FALSE}"; then
  AC_MSG_ERROR([conditional \"$1\" was never defined.
Usually this means the macro was only invoked conditionally.])
fi])])

# Do all the work for Automake.                            -*- Autoconf -*-

# This macro actually does too much some checks are only needed if
# your package does certain things.  But this isn't really a big deal.

# Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002
# Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 8

# There are a few dirty hacks below to avoid letting `AC_PROG_CC' be
# written in clear, in which case automake, when reading aclocal.m4,
# will think it sees a *use*, and therefore will trigger all it's
# C support machinery.  Also note that it means that autoscan, seeing
# CC etc. in the Makefile, will ask for an AC_PROG_CC use...


AC_PREREQ([2.52])

# Autoconf 2.50 wants to disallow AM_ names.  We explicitly allow
# the ones we care about.
m4_pattern_allow([^AM_[A-Z]+FLAGS$])dnl

# AM_INIT_AUTOMAKE(PACKAGE, VERSION, [NO-DEFINE])
# AM_INIT_AUTOMAKE([OPTIONS])
# -----------------------------------------------
# The call with PACKAGE and VERSION arguments is the old style
# call (pre autoconf-2.50), which is being phased out.  PACKAGE
# and VERSION should now be passed to AC_INIT and removed from
# the call to AM_INIT_AUTOMAKE.
# We support both call styles for the transition.  After
# the next Automake release, Autoconf can make the AC_INIT
# arguments mandatory, and then we can depend on a new Autoconf
# release and drop the old call support.
AC_DEFUN([AM_INIT_AUTOMAKE],
[AC_REQUIRE([AM_SET_CURRENT_AUTOMAKE_VERSION])dnl
 AC_REQUIRE([AC_PROG_INSTALL])dnl
# test to see if srcdir already configured
if test "`cd $srcdir && pwd`" != "`pwd`" &&
   test -f $srcdir/config.status; then
  AC_MSG_ERROR([source directory already configured; run "make distclean" there first])
fi

# Define the identity of the package.
dnl Distinguish between old-style and new-style calls.
m4_ifval([$2],
[m4_ifval([$3], [_AM_SET_OPTION([no-define])])dnl
 AC_SUBST([PACKAGE], [$1])dnl
 AC_SUBST([VERSION], [$2])],
[_AM_SET_OPTIONS([$1])dnl
 AC_SUBST([PACKAGE], [AC_PACKAGE_TARNAME])dnl
 AC_SUBST([VERSION], [AC_PACKAGE_VERSION])])dnl

_AM_IF_OPTION([no-define],,
[AC_DEFINE_UNQUOTED(PACKAGE, "$PACKAGE", [Name of package])
 AC_DEFINE_UNQUOTED(VERSION, "$VERSION", [Version number of package])])dnl

# Some tools Automake needs.
AC_REQUIRE([AM_SANITY_CHECK])dnl
AC_REQUIRE([AC_ARG_PROGRAM])dnl
AM_MISSING_PROG(ACLOCAL, aclocal-${am__api_version})
AM_MISSING_PROG(AUTOCONF, autoconf)
AM_MISSING_PROG(AUTOMAKE, automake-${am__api_version})
AM_MISSING_PROG(AUTOHEADER, autoheader)
AM_MISSING_PROG(MAKEINFO, makeinfo)
AM_MISSING_PROG(AMTAR, tar)
AM_PROG_INSTALL_SH
AM_PROG_INSTALL_STRIP
# We need awk for the "check" target.  The system "awk" is bad on
# some platforms.
AC_REQUIRE([AC_PROG_AWK])dnl
AC_REQUIRE([AC_PROG_MAKE_SET])dnl

_AM_IF_OPTION([no-dependencies],,
[AC_PROVIDE_IFELSE([AC_PROG_][CC],
                  [_AM_DEPENDENCIES(CC)],
                  [define([AC_PROG_][CC],
                          defn([AC_PROG_][CC])[_AM_DEPENDENCIES(CC)])])dnl
AC_PROVIDE_IFELSE([AC_PROG_][CXX],
                  [_AM_DEPENDENCIES(CXX)],
                  [define([AC_PROG_][CXX],
                          defn([AC_PROG_][CXX])[_AM_DEPENDENCIES(CXX)])])dnl
])
])

# Copyright 2002  Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA

# AM_AUTOMAKE_VERSION(VERSION)
# ----------------------------
# Automake X.Y traces this macro to ensure aclocal.m4 has been
# generated from the m4 files accompanying Automake X.Y.
AC_DEFUN([AM_AUTOMAKE_VERSION],[am__api_version="1.6"])

# AM_SET_CURRENT_AUTOMAKE_VERSION
# -------------------------------
# Call AM_AUTOMAKE_VERSION so it can be traced.
# This function is AC_REQUIREd by AC_INIT_AUTOMAKE.
AC_DEFUN([AM_SET_CURRENT_AUTOMAKE_VERSION],
	 [AM_AUTOMAKE_VERSION([1.6.3])])

# Helper functions for option handling.                    -*- Autoconf -*-

# Copyright 2001, 2002  Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 2

# _AM_MANGLE_OPTION(NAME)
# -----------------------
AC_DEFUN([_AM_MANGLE_OPTION],
[[_AM_OPTION_]m4_bpatsubst($1, [[^a-zA-Z0-9_]], [_])])

# _AM_SET_OPTION(NAME)
# ------------------------------
# Set option NAME.  Presently that only means defining a flag for this option.
AC_DEFUN([_AM_SET_OPTION],
[m4_define(_AM_MANGLE_OPTION([$1]), 1)])

# _AM_SET_OPTIONS(OPTIONS)
# ----------------------------------
# OPTIONS is a space-separated list of Automake options.
AC_DEFUN([_AM_SET_OPTIONS],
[AC_FOREACH([_AM_Option], [$1], [_AM_SET_OPTION(_AM_Option)])])

# _AM_IF_OPTION(OPTION, IF-SET, [IF-NOT-SET])
# -------------------------------------------
# Execute IF-SET if OPTION is set, IF-NOT-SET otherwise.
AC_DEFUN([_AM_IF_OPTION],
[m4_ifset(_AM_MANGLE_OPTION([$1]), [$2], [$3])])

#
# Check to make sure that the build environment is sane.
#

# Copyright 1996, 1997, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 3

# AM_SANITY_CHECK
# ---------------
AC_DEFUN([AM_SANITY_CHECK],
[AC_MSG_CHECKING([whether build environment is sane])
# Just in case
sleep 1
echo timestamp > conftest.file
# Do `set' in a subshell so we don't clobber the current shell's
# arguments.  Must try -L first in case configure is actually a
# symlink; some systems play weird games with the mod time of symlinks
# (eg FreeBSD returns the mod time of the symlink's containing
# directory).
if (
   set X `ls -Lt $srcdir/configure conftest.file 2> /dev/null`
   if test "$[*]" = "X"; then
      # -L didn't work.
      set X `ls -t $srcdir/configure conftest.file`
   fi
   rm -f conftest.file
   if test "$[*]" != "X $srcdir/configure conftest.file" \
      && test "$[*]" != "X conftest.file $srcdir/configure"; then

      # If neither matched, then we have a broken ls.  This can happen
      # if, for instance, CONFIG_SHELL is bash and it inherits a
      # broken ls alias from the environment.  This has actually
      # happened.  Such a system could not be considered "sane".
      AC_MSG_ERROR([ls -t appears to fail.  Make sure there is not a broken
alias in your environment])
   fi

   test "$[2]" = conftest.file
   )
then
   # Ok.
   :
else
   AC_MSG_ERROR([newly created file is older than distributed files!
Check your system clock])
fi
AC_MSG_RESULT(yes)])

#  -*- Autoconf -*-


# Copyright 1997, 1999, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 3

# AM_MISSING_PROG(NAME, PROGRAM)
# ------------------------------
AC_DEFUN([AM_MISSING_PROG],
[AC_REQUIRE([AM_MISSING_HAS_RUN])
$1=${$1-"${am_missing_run}$2"}
AC_SUBST($1)])


# AM_MISSING_HAS_RUN
# ------------------
# Define MISSING if not defined so far and test if it supports --run.
# If it does, set am_missing_run to use it, otherwise, to nothing.
AC_DEFUN([AM_MISSING_HAS_RUN],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
test x"${MISSING+set}" = xset || MISSING="\${SHELL} $am_aux_dir/missing"
# Use eval to expand $SHELL
if eval "$MISSING --run true"; then
  am_missing_run="$MISSING --run "
else
  am_missing_run=
  AC_MSG_WARN([`missing' script is too old or missing])
fi
])

# AM_AUX_DIR_EXPAND

# Copyright 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# For projects using AC_CONFIG_AUX_DIR([foo]), Autoconf sets
# $ac_aux_dir to `$srcdir/foo'.  In other projects, it is set to
# `$srcdir', `$srcdir/..', or `$srcdir/../..'.
#
# Of course, Automake must honor this variable whenever it calls a
# tool from the auxiliary directory.  The problem is that $srcdir (and
# therefore $ac_aux_dir as well) can be either absolute or relative,
# depending on how configure is run.  This is pretty annoying, since
# it makes $ac_aux_dir quite unusable in subdirectories: in the top
# source directory, any form will work fine, but in subdirectories a
# relative path needs to be adjusted first.
#
# $ac_aux_dir/missing
#    fails when called from a subdirectory if $ac_aux_dir is relative
# $top_srcdir/$ac_aux_dir/missing
#    fails if $ac_aux_dir is absolute,
#    fails when called from a subdirectory in a VPATH build with
#          a relative $ac_aux_dir
#
# The reason of the latter failure is that $top_srcdir and $ac_aux_dir
# are both prefixed by $srcdir.  In an in-source build this is usually
# harmless because $srcdir is `.', but things will broke when you
# start a VPATH build or use an absolute $srcdir.
#
# So we could use something similar to $top_srcdir/$ac_aux_dir/missing,
# iff we strip the leading $srcdir from $ac_aux_dir.  That would be:
#   am_aux_dir='\$(top_srcdir)/'`expr "$ac_aux_dir" : "$srcdir//*\(.*\)"`
# and then we would define $MISSING as
#   MISSING="\${SHELL} $am_aux_dir/missing"
# This will work as long as MISSING is not called from configure, because
# unfortunately $(top_srcdir) has no meaning in configure.
# However there are other variables, like CC, which are often used in
# configure, and could therefore not use this "fixed" $ac_aux_dir.
#
# Another solution, used here, is to always expand $ac_aux_dir to an
# absolute PATH.  The drawback is that using absolute paths prevent a
# configured tree to be moved without reconfiguration.

# Rely on autoconf to set up CDPATH properly.
AC_PREREQ([2.50])

AC_DEFUN([AM_AUX_DIR_EXPAND], [
# expand $ac_aux_dir to an absolute path
am_aux_dir=`cd $ac_aux_dir && pwd`
])

# AM_PROG_INSTALL_SH
# ------------------
# Define $install_sh.

# Copyright 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

AC_DEFUN([AM_PROG_INSTALL_SH],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
install_sh=${install_sh-"$am_aux_dir/install-sh"}
AC_SUBST(install_sh)])

# AM_PROG_INSTALL_STRIP

# Copyright 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# One issue with vendor `install' (even GNU) is that you can't
# specify the program used to strip binaries.  This is especially
# annoying in cross-compiling environments, where the build's strip
# is unlikely to handle the host's binaries.
# Fortunately install-sh will honor a STRIPPROG variable, so we
# always use install-sh in `make install-strip', and initialize
# STRIPPROG with the value of the STRIP variable (set by the user).
AC_DEFUN([AM_PROG_INSTALL_STRIP],
[AC_REQUIRE([AM_PROG_INSTALL_SH])dnl
# Installed binaries are usually stripped using `strip' when the user
# run `make install-strip'.  However `strip' might not be the right
# tool to use in cross-compilation environments, therefore Automake
# will honor the `STRIP' environment variable to overrule this program.
dnl Don't test for $cross_compiling = yes, because it might be `maybe'.
if test "$cross_compiling" != no; then
  AC_CHECK_TOOL([STRIP], [strip], :)
fi
INSTALL_STRIP_PROGRAM="\${SHELL} \$(install_sh) -c -s"
AC_SUBST([INSTALL_STRIP_PROGRAM])])

# serial 4						-*- Autoconf -*-

# Copyright 1999, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.


# There are a few dirty hacks below to avoid letting `AC_PROG_CC' be
# written in clear, in which case automake, when reading aclocal.m4,
# will think it sees a *use*, and therefore will trigger all it's
# C support machinery.  Also note that it means that autoscan, seeing
# CC etc. in the Makefile, will ask for an AC_PROG_CC use...



# _AM_DEPENDENCIES(NAME)
# ----------------------
# See how the compiler implements dependency checking.
# NAME is "CC", "CXX", "GCJ", or "OBJC".
# We try a few techniques and use that to set a single cache variable.
#
# We don't AC_REQUIRE the corresponding AC_PROG_CC since the latter was
# modified to invoke _AM_DEPENDENCIES(CC); we would have a circular
# dependency, and given that the user is not expected to run this macro,
# just rely on AC_PROG_CC.
AC_DEFUN([_AM_DEPENDENCIES],
[AC_REQUIRE([AM_SET_DEPDIR])dnl
AC_REQUIRE([AM_OUTPUT_DEPENDENCY_COMMANDS])dnl
AC_REQUIRE([AM_MAKE_INCLUDE])dnl
AC_REQUIRE([AM_DEP_TRACK])dnl

ifelse([$1], CC,   [depcc="$CC"   am_compiler_list=],
       [$1], CXX,  [depcc="$CXX"  am_compiler_list=],
       [$1], OBJC, [depcc="$OBJC" am_compiler_list='gcc3 gcc'],
       [$1], GCJ,  [depcc="$GCJ"  am_compiler_list='gcc3 gcc'],
                   [depcc="$$1"   am_compiler_list=])

AC_CACHE_CHECK([dependency style of $depcc],
               [am_cv_$1_dependencies_compiler_type],
[if test -z "$AMDEP_TRUE" && test -f "$am_depcomp"; then
  # We make a subdir and do the tests there.  Otherwise we can end up
  # making bogus files that we don't know about and never remove.  For
  # instance it was reported that on HP-UX the gcc test will end up
  # making a dummy file named `D' -- because `-MD' means `put the output
  # in D'.
  mkdir conftest.dir
  # Copy depcomp to subdir because otherwise we won't find it if we're
  # using a relative directory.
  cp "$am_depcomp" conftest.dir
  cd conftest.dir

  am_cv_$1_dependencies_compiler_type=none
  if test "$am_compiler_list" = ""; then
     am_compiler_list=`sed -n ['s/^#*\([a-zA-Z0-9]*\))$/\1/p'] < ./depcomp`
  fi
  for depmode in $am_compiler_list; do
    # We need to recreate these files for each test, as the compiler may
    # overwrite some of them when testing with obscure command lines.
    # This happens at least with the AIX C compiler.
    echo '#include "conftest.h"' > conftest.c
    echo 'int i;' > conftest.h
    echo "${am__include} ${am__quote}conftest.Po${am__quote}" > confmf

    case $depmode in
    nosideeffect)
      # after this tag, mechanisms are not by side-effect, so they'll
      # only be used when explicitly requested
      if test "x$enable_dependency_tracking" = xyes; then
	continue
      else
	break
      fi
      ;;
    none) break ;;
    esac
    # We check with `-c' and `-o' for the sake of the "dashmstdout"
    # mode.  It turns out that the SunPro C++ compiler does not properly
    # handle `-M -o', and we need to detect this.
    if depmode=$depmode \
       source=conftest.c object=conftest.o \
       depfile=conftest.Po tmpdepfile=conftest.TPo \
       $SHELL ./depcomp $depcc -c conftest.c -o conftest.o >/dev/null 2>&1 &&
       grep conftest.h conftest.Po > /dev/null 2>&1 &&
       ${MAKE-make} -s -f confmf > /dev/null 2>&1; then
      am_cv_$1_dependencies_compiler_type=$depmode
      break
    fi
  done

  cd ..
  rm -rf conftest.dir
else
  am_cv_$1_dependencies_compiler_type=none
fi
])
AC_SUBST([$1DEPMODE], [depmode=$am_cv_$1_dependencies_compiler_type])
])


# AM_SET_DEPDIR
# -------------
# Choose a directory name for dependency files.
# This macro is AC_REQUIREd in _AM_DEPENDENCIES
AC_DEFUN([AM_SET_DEPDIR],
[rm -f .deps 2>/dev/null
mkdir .deps 2>/dev/null
if test -d .deps; then
  DEPDIR=.deps
else
  # MS-DOS does not allow filenames that begin with a dot.
  DEPDIR=_deps
fi
rmdir .deps 2>/dev/null
AC_SUBST([DEPDIR])
])


# AM_DEP_TRACK
# ------------
AC_DEFUN([AM_DEP_TRACK],
[AC_ARG_ENABLE(dependency-tracking,
[  --disable-dependency-tracking Speeds up one-time builds
  --enable-dependency-tracking  Do not reject slow dependency extractors])
if test "x$enable_dependency_tracking" != xno; then
  am_depcomp="$ac_aux_dir/depcomp"
  AMDEPBACKSLASH='\'
fi
AM_CONDITIONAL([AMDEP], [test "x$enable_dependency_tracking" != xno])
AC_SUBST([AMDEPBACKSLASH])
])

# Generate code to set up dependency tracking.   -*- Autoconf -*-

# Copyright 1999, 2000, 2001, 2002 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

#serial 2

# _AM_OUTPUT_DEPENDENCY_COMMANDS
# ------------------------------
AC_DEFUN([_AM_OUTPUT_DEPENDENCY_COMMANDS],
[for mf in $CONFIG_FILES; do
  # Strip MF so we end up with the name of the file.
  mf=`echo "$mf" | sed -e 's/:.*$//'`
  # Check whether this is an Automake generated Makefile or not.
  # We used to match only the files named `Makefile.in', but
  # some people rename them; so instead we look at the file content.
  # Grep'ing the first line is not enough: some people post-process
  # each Makefile.in and add a new line on top of each file to say so.
  # So let's grep whole file.
  if grep '^#.*generated by automake' $mf > /dev/null 2>&1; then
    dirpart=`AS_DIRNAME("$mf")`
  else
    continue
  fi
  grep '^DEP_FILES *= *[[^ @%:@]]' < "$mf" > /dev/null || continue
  # Extract the definition of DEP_FILES from the Makefile without
  # running `make'.
  DEPDIR=`sed -n -e '/^DEPDIR = / s///p' < "$mf"`
  test -z "$DEPDIR" && continue
  # When using ansi2knr, U may be empty or an underscore; expand it
  U=`sed -n -e '/^U = / s///p' < "$mf"`
  test -d "$dirpart/$DEPDIR" || mkdir "$dirpart/$DEPDIR"
  # We invoke sed twice because it is the simplest approach to
  # changing $(DEPDIR) to its actual value in the expansion.
  for file in `sed -n -e '
    /^DEP_FILES = .*\\\\$/ {
      s/^DEP_FILES = //
      :loop
	s/\\\\$//
	p
	n
	/\\\\$/ b loop
      p
    }
    /^DEP_FILES = / s/^DEP_FILES = //p' < "$mf" | \
       sed -e 's/\$(DEPDIR)/'"$DEPDIR"'/g' -e 's/\$U/'"$U"'/g'`; do
    # Make sure the directory exists.
    test -f "$dirpart/$file" && continue
    fdir=`AS_DIRNAME(["$file"])`
    AS_MKDIR_P([$dirpart/$fdir])
    # echo "creating $dirpart/$file"
    echo '# dummy' > "$dirpart/$file"
  done
done
])# _AM_OUTPUT_DEPENDENCY_COMMANDS


# AM_OUTPUT_DEPENDENCY_COMMANDS
# -----------------------------
# This macro should only be invoked once -- use via AC_REQUIRE.
#
# This code is only required when automatic dependency tracking
# is enabled.  FIXME.  This creates each `.P' file that we will
# need in order to bootstrap the dependency handling code.
AC_DEFUN([AM_OUTPUT_DEPENDENCY_COMMANDS],
[AC_CONFIG_COMMANDS([depfiles],
     [test x"$AMDEP_TRUE" != x"" || _AM_OUTPUT_DEPENDENCY_COMMANDS],
     [AMDEP_TRUE="$AMDEP_TRUE" ac_aux_dir="$ac_aux_dir"])
])

# Copyright 2001 Free Software Foundation, Inc.             -*- Autoconf -*-

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 2

# AM_MAKE_INCLUDE()
# -----------------
# Check to see how make treats includes.
AC_DEFUN([AM_MAKE_INCLUDE],
[am_make=${MAKE-make}
cat > confinc << 'END'
doit:
	@echo done
END
# If we don't find an include directive, just comment out the code.
AC_MSG_CHECKING([for style of include used by $am_make])
am__include="#"
am__quote=
_am_result=none
# First try GNU make style include.
echo "include confinc" > confmf
# We grep out `Entering directory' and `Leaving directory'
# messages which can occur if `w' ends up in MAKEFLAGS.
# In particular we don't look at `^make:' because GNU make might
# be invoked under some other name (usually "gmake"), in which
# case it prints its new name instead of `make'.
if test "`$am_make -s -f confmf 2> /dev/null | fgrep -v 'ing directory'`" = "done"; then
   am__include=include
   am__quote=
   _am_result=GNU
fi
# Now try BSD make style include.
if test "$am__include" = "#"; then
   echo '.include "confinc"' > confmf
   if test "`$am_make -s -f confmf 2> /dev/null`" = "done"; then
      am__include=.include
      am__quote="\""
      _am_result=BSD
   fi
fi
AC_SUBST(am__include)
AC_SUBST(am__quote)
AC_MSG_RESULT($_am_result)
rm -f confinc confmf
])

dnl @synopsis TAC_ARG_CONFIG_MPI
dnl
dnl Test a variety of MPI options:
dnl --enable-mpi       - Turns MPI compiling mode on
dnl --with-mpi         - specify root directory of MPI
dnl --with-mpi-compilers - Turns on MPI compiling mode and sets the MPI C++
dnl                       compiler = mpicxx or mpiCC (if mpicxx not available),
dnl                       the MPI C compiler = mpicc and 
dnl                       the MPI Fortran compiler = mpif77
dnl --with-mpi-incdir - specify include directory for MPI 
dnl --with-mpi-libs    - specify MPI libraries
dnl --with-mpi-libdir  - specify location of MPI libraries
dnl
dnl If any of these options are set, HAVE_MPI will be defined for both
dnl Autoconf and Automake, and HAVE_MPI will be defined in the
dnl generated config.h file
dnl
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_CONFIG_MPI],
[

AC_ARG_ENABLE(mpi,
[AC_HELP_STRING([--enable-mpi],[MPI support])],
[HAVE_PKG_MPI=$enableval],
[HAVE_PKG_MPI=no]
)

AC_ARG_WITH(mpi-compilers,
[AC_HELP_STRING([--with-mpi-compilers=PATH],
[use MPI compilers mpicc, mpif77, and mpicxx (or mpiCC) in the specified path or in the default path if no path is specified. Enables MPI])],
[
  if test X${withval} != Xno; then
    HAVE_PKG_MPI=yes
    if test X${withval} = Xyes; then
      # Check for mpicxx, if it does not exist, use mpiCC instead.
      AC_CHECK_PROG(MPI_CXX, mpicxx, mpicxx, mpiCC)
      MPI_CC=mpicc
      MPI_F77=mpif77
    else
      MPI_TEMP_CXX=${withval}/mpicxx
      if test -f ${MPI_TEMP_CXX}; then
        MPI_CXX=${MPI_TEMP_CXX}
      else
        MPI_CXX=${withval}/mpiCC
      fi
      MPI_CC=${withval}/mpicc
      MPI_F77=${withval}/mpif77
    fi
  fi
]
)

AC_ARG_WITH(mpi,
[AC_HELP_STRING([--with-mpi=MPIROOT],[use MPI root directory (enables MPI)])],
[
  HAVE_PKG_MPI=yes
  MPI_DIR=${withval}
  AC_MSG_CHECKING(MPI directory)
  AC_MSG_RESULT([${MPI_DIR}])
]
)

#AC_ARG_WITH(mpi-include,
#[AC_HELP_STRING([--with-mpi-include],[Obsolete.  Use --with-mpi-incdir=DIR instead.  Do not prefix DIR with '-I'.])],
#[AC_MSG_ERROR([--with-mpi-include is an obsolte option.   Use --with-mpi-incdir=DIR instead.  Do not prefix DIR with '-I'.  For example '--with-mpi-incdir=/usr/lam_path/include'.])]
#)

AC_ARG_WITH(mpi-libs,
[AC_HELP_STRING([--with-mpi-libs="LIBS"],[MPI libraries @<:@"-lmpi"@:>@])],
[
  MPI_LIBS=${withval}
  AC_MSG_CHECKING(user-defined MPI libraries)
  AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_ARG_WITH(mpi-incdir,
[AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@  Do not use -I])],
[
  MPI_INC=${withval}
  AC_MSG_CHECKING(user-defined MPI includes)
  AC_MSG_RESULT([${MPI_INC}])
]
)

AC_ARG_WITH(mpi-libdir,
[AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@  Do not use -L])],
[
  MPI_LIBDIR=${withval}
  AC_MSG_CHECKING(user-defined MPI library directory)
  AC_MSG_RESULT([${MPI_LIBDIR}])
]
)

AC_MSG_CHECKING(whether we are using MPI)
AC_MSG_RESULT([${HAVE_PKG_MPI}])

if test "X${HAVE_PKG_MPI}" = "Xyes"; then
   AC_DEFINE(HAVE_MPI,,[define if we want to use MPI])
fi

dnl Define Automake version of HAVE_MPI if appropriate

AM_CONDITIONAL(HAVE_MPI, [test "X${HAVE_PKG_MPI}" = "Xyes"])


dnl
dnl --------------------------------------------------------------------
dnl Check for MPI compilers (must be done *before* AC_PROG_CXX,
dnl AC_PROG_CC and AC_PROG_F77)
dnl 
dnl --------------------------------------------------------------------

if test -n "${MPI_CXX}"; then
  if test -f ${MPI_CXX}; then
    MPI_CXX_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_CXX_EXISTS, ${MPI_CXX}, yes, no)
  fi

  if test "X${MPI_CXX_EXISTS}" = "Xyes"; then
    CXX=${MPI_CXX}
  else
    echo "-----"
    echo "Cannot find MPI C++ compiler ${MPI_CXX}."
    echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH"
    echo "or specify a C++ compiler using CXX=<compiler>"
    echo "Do not use --with-mpi-compilers if using CXX=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI C++ compiler (${MPI_CXX}) not found.])
  fi
fi

if test -n "${MPI_CC}"; then
  if test -f ${MPI_CC}; then
    MPI_CC_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_CC_EXISTS, ${MPI_CC}, yes, no)
  fi

  if test "X${MPI_CC_EXISTS}" = "Xyes"; then
    CC=${MPI_CC}
  else
    echo "-----"
    echo "Cannot find MPI C compiler ${MPI_CC}."
    echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH"
    echo "or specify a C compiler using CC=<compiler>"
    echo "Do not use --with-mpi-compilers if using CC=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI C compiler (${MPI_CC}) not found.])
  fi
fi

if test -n "${MPI_F77}"; then
  if test -f ${MPI_F77}; then
    MPI_F77_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_F77_EXISTS, ${MPI_F77}, yes, no)
  fi

  if test "X${MPI_F77_EXISTS}" = "Xyes"; then
    F77=${MPI_F77}
  else
    echo "-----"
    echo "Cannot find MPI Fortran compiler ${MPI_F77}."
    echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH"
    echo "or specify a Fortran 77 compiler using F77=<compiler>"
    echo "Do not use --with-mpi-compilers if using F77=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI Fortran 77 compiler (${MPI_F77}) not found.])
  fi
fi
])

dnl @synopsis TAC_ARG_WITH_AR
dnl
dnl Test for --with-ar="ar_program ar_flags".
dnl Default is "ar cru"
dnl 
dnl Generates an Automake conditional USE_ALTERNATE_AR that can be tested.  
dnl Generates the user-specified archiver command in @ALTERNATE_AR@.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_AR],
[
AC_ARG_WITH(ar,
AC_HELP_STRING([--with-ar], [override archiver command (default is "ar cru")]),
[
AC_MSG_CHECKING(user-defined archiver)
AC_MSG_RESULT([${withval}])
USE_ALTERNATE_AR=yes
ALTERNATE_AR="${withval}"
]
)

if test -n "${SPECIAL_AR}" && test "X${USE_ALTERNATE_AR}" != "Xyes";
then
  USE_ALTERNATE_AR=yes
  ALTERNATE_AR="${SPECIAL_AR}"
fi

AC_MSG_CHECKING(for special archiver command)
if test "X${USE_ALTERNATE_AR}" = "Xyes"; then
   AC_MSG_RESULT([${ALTERNATE_AR}])
   AM_CONDITIONAL(USE_ALTERNATE_AR, true)
else
   AC_MSG_RESULT([none])
   AM_CONDITIONAL(USE_ALTERNATE_AR, false)
fi
AC_SUBST(ALTERNATE_AR)
])


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
dnl @synopsis TAC_ARG_ENABLE_FEATURE(FEATURE_NAME, FEATURE_DESCRIPTION, HAVE_NAME, DEFAULT_VAL)
dnl
dnl Test for --enable-${FEATURE_NAME} and set to DEFAULT_VAL value if feature not specified.
dnl Also calls AC_DEFINE to define HAVE_${HAVE_NAME} if value is not equal to "no"
dnl 
dnl Use this macro to help defining whether or not optional 
dnl features* should compiled.  For example:
dnl
dnl TAC_ARG_ENABLE_FEATURE(epetra, [Configure and build epetra], EPETRA, yes)
dnl 
dnl will test for --enable-epetra when configure is run.  If it is defined 
dnl and not set to "no" or not defined (default is "yes") then HAVE_EPETRA will
dnl be defined, if --enable-epetra is defined to be "no", HAVE_EPETRA will not
dnl be defined.
dnl
dnl *NOTE: epetra, aztecoo, komplex, ifpack, and other software found in
dnl subdirectories of Trilinos/packages are "packages" in their own right.
dnl However, these packages are also "features" of the larger package
dnl "Trilinos".  Therefore, when configuring from the Trilinos directory,
dnl it is appropriate to refer to these software packages as "features".
dnl
dnl This file was based on tac_arg_with_package.m4 by Mike Heroux
dnl @author James Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_FEATURE],
[
AC_ARG_ENABLE([$1],
AC_HELP_STRING([--enable-$1],[$2 (default is [$4])]),
ac_cv_use_$1=$enableval, ac_cv_use_$1=$4)

AC_MSG_CHECKING(whether to use [$1])

if test "X$ac_cv_use_$1" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([HAVE_$3],,[Define if want to build $1])
else
  AC_MSG_RESULT(no)
fi
])


AC_DEFUN([TAC_ENABLE_PYTHON],
[
AZ_PYTHON_DEFAULT( )  # Default to no python support
AZ_PYTHON_ENABLE( )   # Check for --enable-python, --disable-python
AZ_PYTHON_WITH( )     # Check for --with-python, --without-python

AC_MSG_CHECKING(whether we should build the python wrappers)
if test -n "$PYTHON"; then
  BUILD_PYTHON=yes
  AC_MSG_RESULT(yes)
  AM_CONDITIONAL(BUILD_PYTHON, true)

  # Ensure that we have python version 2.2 or greater (for distutils)
  AZ_PYTHON_VERSION_ENSURE( [2.2] )

  # Python compiler and linker flags
  AZ_PYTHON_CSPEC
  AZ_PYTHON_LSPEC

  # Check that Python.h is available
  save_CPPFLAGS=$CPPFLAGS
  CPPFLAGS="$save_CPPFLAGS $PYTHON_CSPEC"
  AC_LANG([C++])
  AC_CHECK_HEADER(
  [Python.h],
  break,
  AC_MSG_ERROR([You must have Python.h in order to build the Python support!!]))
  CPPFLAGS="$save_CPPFLAGS"

  # Check for Numeric
  #AC_PYTHON_MODULE(Numeric,yes)

  # If user specifies prefix, use it for the PYTHON_PREFIX
  if test "$prefix" != "$ac_default_prefix"; then
    PYTHON_PREFIX=$prefix
  fi
  if test "$exec_prefix" != "$ac_default_prefix"; then
    PYTHON_EXECPREFIX=$exec_prefix
  fi

else
  AC_MSG_RESULT(no)
  BUILD_PYTHON=no
  AM_CONDITIONAL(BUILD_PYTHON, false)
fi

# ------------------------------------------------------------------------
# If the python wrappers are to be built, then SWIG (Simple Wrapper
# Interface Generator) is required
# ------------------------------------------------------------------------

if test -n "$PYTHON"; then

  # Check for --with-swig[=path]
  AC_MSG_CHECKING(for --with-swig)
  AC_ARG_WITH(swig,
              [AC_HELP_STRING([--with-swig@<:@=SWIG@:>@],
                              [enable swig and set swig binary])],
              [AC_MSG_RESULT(yes)
               WITH_SWIG=yes
               if test X${withval} != Xyes ; then
                 SWIG=$withval
               fi],
              [AC_MSG_RESULT(no)
               AC_CHECK_PROG(WITH_SWIG,swig,yes,no)])

  # Report error if no swig found
  if test ${WITH_SWIG} = no; then
     AC_MSG_ERROR(
     [Python wrappers require swig (Simple Wrapper Interface Generator).
      See http://www.swig.org])
  fi

  # SWIG configuration
  AC_PROG_SWIG(1.3.23)
  SWIG_ENABLE_CXX
  SWIG_MULTI_MODULE_SUPPORT
  SWIG_PYTHON
fi
AM_CONDITIONAL(HAVE_SWIG,test X${WITH_SWIG} = Xyes)
AC_SUBST(GNU_HAVE_SWIG, ${WITH_SWIG})

])

dnl @synopsis AZ_PYTHON_DEFAULT
dnl @synopsis AZ_PYTHON_ENABLE
dnl @synopsis AZ_PYTHON_WITH
dnl @synopsis AZ_PYTHON_PATH
dnl @synopsis AZ_PYTHON_VERSION_ENSURE( [2.2] )
dnl @synopsis AZ_PYTHON_CSPEC
dnl @synopsis AZ_PYTHON_LSPEC
dnl
dnl @summary New and revised Python support.
dnl
dnl This file provides autoconf support for those applications that
dnl want to embed python. It supports all pythons >= 2.2 which is the
dnl first official release containing distutils. Version 2.2 of python
dnl was released December 21, 2001. Since it actually executes the
dnl python, cross platform configuration will probably not work. Also,
dnl most of the platforms supported are consistent until you look into
dnl MacOSX. The python included with it is installed as a framework
dnl which is a very different environment to set up the normal tools
dnl such as gcc and libtool to deal with. Therefore, once we establish
dnl which python that we are going to use, we use its distutils to
dnl actually compile and link our modules or applications.
dnl
dnl At this time, it does NOT support linking with Python statically.
dnl It does support dynamic linking.
dnl
dnl This set of macros help define $PYTHON, $PYTHON_USE, $PYTHON_CSPEC
dnl and $PYTHON_LSPEC. $PYTHON defines the full executable path for the
dnl Python being linked to and is used within these macros to determine
dnl if that has been specified or found. These macros do execute this
dnl python version so it must be present on the system at configure
dnl time.
dnl
dnl $PYTHON_USE is an automake variable that defines whether Python
dnl support should be included or not in your application.
dnl $PYTHON_CSPEC is a variable that supplies additional CFLAGS for the
dnl compilation of the application/shared library. $PYTHON_LSPEC is a
dnl variable that supplies additional LDFLAGS for linking the
dnl application/shared library.
dnl
dnl The following is an example of how to set up for python usage
dnl within your application in your configure.in:
dnl
dnl   AZ_PYTHON_DEFAULT( )
dnl   AZ_PYTHON_ENABLE( )             # Optional
dnl   AZ_PYTHON_WITH( )               # Optional
dnl   AZ_PYTHON_PATH( )               # or AZ_PYTHON_INSIST( )
dnl   # if $PYTHON is not defined, then the following do nothing.
dnl   AZ_PYTHON_VERSION_ENSURE( [2.2] )
dnl   AZ_PYTHON_CSPEC
dnl   AZ_PYTHON_LSPEC
dnl
dnl The AZ_PYTHON_DEFAULT sets the $PYTHON_USE to false. Thereby,
dnl excluding it if it was optional.
dnl
dnl The AZ_PYTHON_ENABLE looks for the optional configure parameters of
dnl --enable-python/--disable-python and establishes the $PYTHON and
dnl $PYTHON_USE variables accordingly.
dnl
dnl The AZ_PYTHON_WITH looks for the optional configure parameters of
dnl --with-python/--without-python and establishes the $PYTHON and
dnl $PYTHON_USE variables accordingly.
dnl
dnl The AZ_PYTHON_PATH looks for python assuming that none has been
dnl previously found or defined and issues an error if it does not find
dnl it. If it does find it, it establishes the $PYTHON and $PYTHON_USE
dnl variables accordingly. AZ_PYTHON_INSIST could be used here instead
dnl if you want to insist that Python support be included using the
dnl --enable-python or --with-python checks previously done.
dnl
dnl The AZ_PYTHON_VERSION_ENSURE issues an error if the Python
dnl previously found is not of version 2.2 or greater.
dnl
dnl Once that these macros have be run, we can use PYTHON_USE within
dnl the makefile.am file to conditionally add the Python support such
dnl as:
dnl
dnl Makefile.am example showing optional inclusion of directories:
dnl
dnl  if PYTHON_USE
dnl  plugins = plugins
dnl  src = src
dnl  else
dnl  plugins =
dnl  src =
dnl  endif
dnl
dnl  SUBDIRS = . $(plugins) $(src)
dnl
dnl Makefile.am example showing optional shared library build:
dnl
dnl  if PYTHON_USE
dnl  lib_LTLIBRARIES        = libElemList.la
dnl  libElemList_la_SOURCES = libElemList.c
dnl  libElemList_la_CFLAGS  = @PYTHON_CSPEC@
dnl  libElemList_la_LDFLAGS = @PYTHON_LSPEC@
dnl  endif
dnl
dnl Makefile.am example showing optional program build:
dnl
dnl  if PYTHON_USE
dnl  bin_PROGRAMS    = runFunc
dnl  runFunc_SOURCES = runFunc.c
dnl  runFunc_CFLAGS  = @PYTHON_CSPEC@
dnl  runFunc_LDFLAGS = @PYTHON_LSPEC@
dnl  endif
dnl
dnl The above compiles the modules only if PYTHON_USE was specified as
dnl true. Also, the else portion of the if was optional.
dnl
dnl @category InstalledPackages
dnl @author Robert White <kranki@mac.com>
dnl @author Dustin Mitchell <dustin@cs.uchicago.edu>
dnl @version 2005-01-14
dnl @license GPLWithACException

# AZ_PYTHON_DEFAULT( )
# -----------------
# Sets the default to not include Python support.

AC_DEFUN([AZ_PYTHON_DEFAULT],
[
    az_python_use=false
    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
])



# AZ_PYTHON_ENABLE( [path] )
# -----------------------------------------------------------------
# Handles the various --enable-python commands.
# Input:
#   $1 is the optional search path for the python executable if needed
# Ouput:
#   PYTHON_USE (AM_CONDITIONAL) is true if python executable found
#   and --enable-python was requested; otherwise false.
#   $PYTHON contains the full executable path to python if PYTHON_ENABLE_USE
#   is true.
#
# Example:
#   AZ_PYTHON_ENABLE( )
#   or
#   AZ_PYTHON_ENABLE( "/usr/bin" )

AC_DEFUN([AZ_PYTHON_ENABLE],
[
    AC_ARG_VAR([PYTHON],[Python Executable Path])

    # unless PYTHON was supplied to us (as a precious variable),
    # see if --enable-python[=PythonExecutablePath], --enable-python,
    # --disable-python or --enable-python=no was given.
    if test -z "$PYTHON"
    then
        AC_MSG_CHECKING(for --enable-python)
        AC_ARG_ENABLE(
            python,
            AC_HELP_STRING([--enable-python@<:@=PYTHON@:>@],
                [absolute path name of Python executable]
            ),
            [
                if test "$enableval" = "yes"
                then
                    # "yes" was specified, but we don't have a path
                    # for the executable.
                    # So, let's searth the PATH Environment Variable.
                    AC_MSG_RESULT(yes)
                    AC_PATH_PROG(
                        [PYTHON],
                        python,
                        [],
                        $1
                    )
                    if test -z "$PYTHON"
                    then
                        AC_MSG_ERROR(no path to python found)
                    fi
                    az_python_use=true
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                    AZ_PYTHON_PREFIX( )
                elif test "$enableval" = "no"
                then
                    AC_MSG_RESULT(no)
                    az_python_use=false
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                else
                    # $enableval must be the executable path then.
                    AC_SUBST([PYTHON], ["${enableval}"])
                    AC_MSG_RESULT($withval)
                    az_python_use=true
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                    AZ_PYTHON_PREFIX( )
                fi
            ],
            [
                # --with-python was not specified.
                AC_MSG_RESULT(no)
                az_python_use=false
                AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
            ]
        )
    fi

])



# AZ_PYTHON_CSPEC( )
# -----------------
# Set up the c compiler options to compile Python
# embedded programs/libraries in $PYTHON_CSPEC if
# $PYTHON has been defined.

AC_DEFUN([AZ_PYTHON_CSPEC],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -n "$PYTHON"
    then
        az_python_prefix=`${PYTHON} -c "import sys; print sys.prefix"`
        if test -z "$az_python_prefix"
        then
            AC_MSG_ERROR([Python Prefix is not known])
        fi
        az_python_execprefix=`${PYTHON} -c "import sys; print sys.exec_prefix"`
        az_python_version=`$PYTHON -c "import sys; print sys.version[[:3]]"`
        az_python_includespec="-I${az_python_prefix}/include/python${az_python_version}"
        if test x"$python_prefix" != x"$python_execprefix"; then
            az_python_execspec="-I${az_python_execprefix}/include/python${az_python_version}"
            az_python_includespec="${az_python_includespec} $az_python_execspec"
        fi
        az_python_ccshared=`${PYTHON} -c "import distutils.sysconfig; print distutils.sysconfig.get_config_var('CFLAGSFORSHARED')"`
        az_python_cspec="${az_python_ccshared} ${az_python_includespec}"
        AC_SUBST([PYTHON_CSPEC], [${az_python_cspec}])
        AC_MSG_NOTICE([PYTHON_CSPEC=${az_python_cspec}])
    fi
])



# AZ_PYTHON_INSIST( )
# -----------------
# Look for Python and set the output variable 'PYTHON'
# to 'python' if found, empty otherwise.

AC_DEFUN([AZ_PYTHON_PATH],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable not found])
    fi
])



# AZ_PYTHON_LSPEC( )
# -----------------
# Set up the linker options to link Python embedded
# programs/libraries in $PYTHON_LSPEC if $PYTHON
# has been defined.

AC_DEFUN([AZ_PYTHON_LSPEC],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -n "$PYTHON"
    then
        AZ_PYTHON_RUN([
import sys
import distutils.sysconfig
strUseFrameWork = "--enable-framework"
dictConfig = distutils.sysconfig.get_config_vars( )
strConfigArgs = dictConfig.get("CONFIG_ARGS")
strLinkSpec =  dictConfig.get('LDFLAGS')
if -1 ==  strConfigArgs.find(strUseFrameWork):
    strLibPL = dictConfig.get("LIBPL")
    if strLibPL and (strLibPL != ""):
        strLinkSpec += " -L%s" % (strLibPL)
    strSys = dictConfig.get("SYSLIBS")
    if strSys and (strSys != ""):
        strLinkSpec += " %s" % (strSys)
    strSHL = dictConfig.get("SHLIBS")
    if strSHL and (strSHL != ""):
        strLinkSpec += " %s" % (strSHL)
    # Construct the Python Library Name.
    strTmplte = " -lpython%d.%d"
    if (sys.platform == "win32") or (sys.platform == "os2emx"):
        strTmplte = " -lpython%d%d"
    strWrk = strTmplte % ( (sys.hexversion >> 24),
                            ((sys.hexversion >> 16) & 0xff))
    strLinkSpec += strWrk
else:
    # This is not ideal since it changes the search path
    # for Frameworks which could have side-effects on
    # other included Frameworks.  However, it is necessary
    # where someone has installed more than one frameworked
    # Python.  Frameworks are really only used in MacOSX.
    strLibFW = dictConfig.get("PYTHONFRAMEWORKPREFIX")
    if strLibFW and (strLibFW != ""):
        strLinkSpec += " -F%s" % (strLibFW)
strLinkSpec += " %s" % (dictConfig.get('LINKFORSHARED'))
print strLinkSpec
        ])
        AC_SUBST([PYTHON_LSPEC], [${az_python_output}])
        AC_MSG_NOTICE([PYTHON_LSPEC=${az_python_output}])
    fi
])



# AZ_PYTHON_PATH( )
# -----------------
# Look for Python and set the output variable 'PYTHON'
# to 'python' if found, empty otherwise.

AC_DEFUN([AZ_PYTHON_PATH],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    AC_PATH_PROG( PYTHON, python, [], $1 )
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable not found])
    else
        az_python_use=true
    fi
    AM_CONDITIONAL(PYTHON_USE, test "$az_python_use" = "true")
])



# AZ_PYTHON_PREFIX( )
# -------------------
# Use the values of $prefix and $exec_prefix for the corresponding
# values of PYTHON_PREFIX and PYTHON_EXEC_PREFIX.

AC_DEFUN([AZ_PYTHON_PREFIX],
[
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable Path is not known])
    fi
    ax_python_prefix=`${PYTHON} -c "import sys; print sys.prefix"`
    ax_python_execprefix=`${PYTHON} -c "import sys; print sys.exec_prefix"`
    AC_SUBST([PYTHON_PREFIX], ["${ax_python_prefix}"])
    AC_SUBST([PYTHON_EXECPREFIX], ["${ax_python_execprefix}"])
])



# AZ_PYTHON_RUN( PYTHON_PROGRAM )
# -----------------
# Run a Python Test Program saving its output
# in az_python_output and its condition code
# in az_python_cc.

AC_DEFUN([AZ_PYTHON_RUN],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable not found])
    else
        cat >conftest.py <<_ACEOF
$1
_ACEOF
        az_python_output=`$PYTHON conftest.py`
        az_python_cc=$?
        rm conftest.py
        if test -f "conftest.pyc"
        then
            rm conftest.pyc
        fi
    fi
])



# AZ_PYTHON_VERSION_CHECK( VERSION, [ACTION-IF-TRUE], [ACTION-IF-FALSE] )
# -----------------------------------------------------------------------------
# Run ACTION-IF-TRUE if the Python interpreter has version >= VERSION.
# Run ACTION-IF-FALSE otherwise.
# This test uses sys.hexversion instead of the string equivalant (first
# word of sys.version), in order to cope with versions such as 2.2c1.
# hexversion has been introduced in Python 1.5.2; it's probably not
# worth to support older versions (1.5.1 was released on October 31, 1998).

AC_DEFUN([AZ_PYTHON_VERSION_CHECK],
 [
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -n "$PYTHON"
    then
        AC_MSG_CHECKING([whether $PYTHON version >= $1])
        AZ_PYTHON_RUN([
import sys, string
# split strings by '.' and convert to numeric.  Append some zeros
# because we need at least 4 digits for the hex conversion.
minver = map(int, string.split('$1', '.')) + [[0, 0, 0]]
minverhex = 0
for i in xrange(0, 4): minverhex = (minverhex << 8) + minver[[i]]
if sys.hexversion >= minverhex:
    sys.exit( 0 )
else:
    sys.exit( 1 )
        ])
        if test $az_python_cc -eq 0
        then
            $2
        m4_ifvaln(
            [$3],
            [else $3]
        )
        fi
    fi
])



# AZ_PYTHON_VERSION_ENSURE( VERSION )
# -----------------
# Insure that the Python Interpreter Version
# is greater than or equal to the VERSION
# parameter.

AC_DEFUN([AZ_PYTHON_VERSION_ENSURE],
[
    AZ_PYTHON_VERSION_CHECK(
        [$1],
        [AC_MSG_RESULT(yes)],
        [AC_MSG_ERROR(too old)]
    )
])



# AZ_PYTHON_WITH( [path] )
# -----------------------------------------------------------------
# Handles the various --with-python commands.
# Input:
#   $1 is the optional search path for the python executable if needed
# Ouput:
#   PYTHON_USE (AM_CONDITIONAL) is true if python executable found
#   and --with-python was requested; otherwise false.
#   $PYTHON contains the full executable path to python if PYTHON_USE
#   is true.
#
# Example:
#   AZ_PYTHON_WITH( )
#   or
#   AZ_PYTHON_WITH("/usr/bin")

AC_DEFUN([AZ_PYTHON_WITH],
[
    AC_ARG_VAR([PYTHON],[Python Executable Path])

    # unless PYTHON was supplied to us (as a precious variable),
    # see if --with-python[=PythonExecutablePath], --with-python,
    # --without-python or --with-python=no was given.
    if test -z "$PYTHON"
    then
        AC_MSG_CHECKING(for --with-python)
        AC_ARG_WITH(
            python,
            AC_HELP_STRING([--with-python@<:@=PYTHON@:>@],
                [absolute path name of Python executable]
            ),
            [
                if test "$withval" = "yes"
                then
                    # "yes" was specified, but we don't have a path
                    # for the executable.
                    # So, let's searth the PATH Environment Variable.
                    AC_MSG_RESULT(yes)
                    AC_PATH_PROG(
                        [PYTHON],
                        python,
                        [],
                        $1
                    )
                    if test -z "$PYTHON"
                    then
                        AC_MSG_ERROR(no path to python found)
                    fi
                    az_python_use=true
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                    AZ_PYTHON_PREFIX( )
                elif test "$withval" = "no"
                then
                    AC_MSG_RESULT(no)
                    az_python_use=false
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                else
                    # $withval must be the executable path then.
                    AC_SUBST([PYTHON], ["${withval}"])
                    AC_MSG_RESULT($withval)
                    az_python_use=true
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                    AZ_PYTHON_PREFIX( )
                fi
            ],
            [
                # --with-python was not specified.
                AC_MSG_RESULT(no)
                az_python_use=false
                AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
            ]
        )
    fi

])

dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_pkg_swig.html
dnl
AC_DEFUN([AC_PROG_SWIG],[
	if test -z "$SWIG"; then
        	AC_PATH_PROG([SWIG],[swig])
	fi
        if test -z "$SWIG" ; then
                AC_MSG_WARN([cannot find 'swig' program. You should look at http://www.swig.org])
                SWIG='echo "Error: SWIG is not installed. You should look at http://www.swig.org" ; false'
        else
	    AC_MSG_NOTICE([SWIG executable is '$SWIG'])
	    if test -n "$1" ; then
                AC_MSG_CHECKING([for SWIG version])
                [swig_version=`$SWIG -version 2>&1 | grep 'SWIG Version' | sed 's/.*\([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/g'`]
                AC_MSG_RESULT([$swig_version])
                if test -n "$swig_version" ; then
                        # Calculate the required version number components
                        [required=$1]
                        [required_major=`echo $required | sed 's/[^0-9].*//'`]
                        if test -z "$required_major" ; then
                                [required_major=0]
                        fi
                        [required=`echo $required | sed 's/[0-9]*[^0-9]//'`]
                        [required_minor=`echo $required | sed 's/[^0-9].*//'`]
                        if test -z "$required_minor" ; then
                                [required_minor=0]
                        fi
                        [required=`echo $required | sed 's/[0-9]*[^0-9]//'`]
                        [required_patch=`echo $required | sed 's/[^0-9].*//'`]
                        if test -z "$required_patch" ; then
                                [required_patch=0]
                        fi
                        # Calculate the available version number components
                        [available=$swig_version]
                        [available_major=`echo $available | sed 's/[^0-9].*//'`]
                        if test -z "$available_major" ; then
                                [available_major=0]
                        fi
                        [available=`echo $available | sed 's/[0-9]*[^0-9]//'`]
                        [available_minor=`echo $available | sed 's/[^0-9].*//'`]
                        if test -z "$available_minor" ; then
                                [available_minor=0]
                        fi
                        [available=`echo $available | sed 's/[0-9]*[^0-9]//'`]
                        [available_patch=`echo $available | sed 's/[^0-9].*//'`]
                        if test -z "$available_patch" ; then
                                [available_patch=0]
                        fi
                        if test $available_major -ne $required_major \
                                -o $available_minor -ne $required_minor \
                                -o $available_patch -lt $required_patch ; then
                                AC_MSG_WARN([SWIG version >= $1 is required.  You have $swig_version.  You should look at http://www.swig.org])
                                SWIG='echo "Error: SWIG version >= $1 is required.  You have '"$swig_version"'.  You should look at http://www.swig.org" ; false'
                        else
                                SWIG_LIB=`$SWIG -swiglib`
                                AC_MSG_NOTICE([SWIG library directory is '$SWIG_LIB'])
                        fi
                else
                        AC_MSG_WARN([cannot determine SWIG version])
                        SWIG='echo "Error: Cannot determine SWIG version.  You should look at http://www.swig.org" ; false'
                fi
            fi
        fi
        AC_SUBST([SWIG_LIB])
])

# SWIG_ENABLE_CXX()
#
# Enable SWIG C++ support.  This affects all invocations of $(SWIG).
AC_DEFUN([SWIG_ENABLE_CXX],[
        AC_REQUIRE([AC_PROG_SWIG])
        AC_REQUIRE([AC_PROG_CXX])
        SWIG="$SWIG -c++"
])

# SWIG_MULTI_MODULE_SUPPORT()
#
# Enable support for multiple modules.  This effects all invocations
# of $(SWIG).  You have to link all generated modules against the
# appropriate SWIG runtime library.  If you want to build Python
# modules for example, use the SWIG_PYTHON() macro and link the
# modules against $(SWIG_PYTHON_LIBS).
#
AC_DEFUN([SWIG_MULTI_MODULE_SUPPORT],[
        AC_REQUIRE([AC_PROG_SWIG])
#        SWIG="$SWIG -noruntime"
#        SWIG="$SWIG -c"
	SWIG="$SWIG"
])

# SWIG_PYTHON([use-shadow-classes = {no, yes}])
#
# Checks for Python and provides the $(SWIG_PYTHON_CPPFLAGS),
# and $(SWIG_PYTHON_OPT) output variables.
#
# $(SWIG_PYTHON_OPT) contains all necessary SWIG options to generate
# code for Python.  Shadow classes are enabled unless the value of the
# optional first argument is exactly 'no'.  If you need multi module
# support (provided by the SWIG_MULTI_MODULE_SUPPORT() macro) use
# $(SWIG_PYTHON_LIBS) to link against the appropriate library.  It
# contains the SWIG Python runtime library that is needed by the type
# check system for example.
AC_DEFUN([SWIG_PYTHON],[
        AC_REQUIRE([AC_PROG_SWIG])
#        AC_REQUIRE([AC_PYTHON_DEVEL])
#        test "x$1" != "xno" || swig_shadow=" -noproxy"
#        AC_SUBST([SWIG_PYTHON_OPT],[-python$swig_shadow])
        AC_SUBST([SWIG_PYTHON_OPT],[-python])
        AC_SUBST([SWIG_PYTHON_CPPFLAGS],[$PYTHON_CPPFLAGS])
])


dnl @synopsis AC_LIB_WAD
dnl
dnl This macro searches for installed WAD library.
dnl
AC_DEFUN([AC_LIB_WAD],
[
        AC_REQUIRE([AC_PYTHON_DEVEL])
        AC_ARG_ENABLE(wad,
        AC_HELP_STRING([--enable-wad], [enable wad module]),
        [
                case "${enableval}" in
                        no)     ;;
                        *)      if test "x${enableval}" = xyes;
                                then
                                        check_wad="yes"
                                fi ;;
                esac
        ], [])

        if test -n "$check_wad";
        then
                AC_CHECK_LIB(wadpy, _init, [WADPY=-lwadpy], [], $PYTHON_LDFLAGS $PYTHON_EXTRA_LIBS)
                AC_SUBST(WADPY)
        fi
])

dnl @synopsis TAC_ARG_WITH_LIBDIRS
dnl
dnl Test for --with-libdirs="-Llibdir1 -Llibdir2".  if defined, 
dnl prepend "-Llibdir1 -Llibdir2" to LDFLAGS
dnl
dnl Use this macro to facilitate addition of directories to library search path.
dnl 
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_LIBDIRS],
[
AC_MSG_CHECKING([whether additional library search paths defined])
AC_ARG_WITH(libdirs,
AC_HELP_STRING([--with-libdirs], 
[additional directories containing libraries: will prepend to search here for libraries, use -Ldir format]),
[
LDFLAGS="${withval} ${LDFLAGS}"
AC_MSG_RESULT([${withval}])
],
AC_MSG_RESULT(no)
)
])


dnl @synopsis TAC_ARG_WITH_INCDIRS
dnl
dnl Test for --with-incdirs="-Iincdir1 -Iincdir2".  if defined, prepend 
dnl "-Iincdir1 -Iincdir2" to CPPFLAGS
dnl
dnl Use this macro to facilitate addition of directories to include file search path.
dnl 
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_INCDIRS],
[
AC_MSG_CHECKING([whether additional include search paths defined])
AC_ARG_WITH(incdirs,
AC_HELP_STRING([--with-incdirs], 
[additional directories containing include files: will prepend to search here for includes, use -Idir format]),
[
CPPFLAGS="${withval} ${CPPFLAGS}"
AC_MSG_RESULT([${withval}])
],
AC_MSG_RESULT(no)
)
])


