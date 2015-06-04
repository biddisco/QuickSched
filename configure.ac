
# This file is part of QuickSched.
# Coypright (c) 2013 pedro.gonnet@durham.ac.uk.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# init the project
AC_INIT([QUICKSCHED],[0.1.0])
AC_CONFIG_SRCDIR([src/qsched.c])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE

# Add local macro collection
AC_CONFIG_MACRO_DIR([m4])

# Stop default CFLAGS from anyone except the environment.
: ${CFLAGS=""}

# generate header file
AM_CONFIG_HEADER(config.h)

# Find and test the compiler.
AX_CHECK_ENABLE_DEBUG
AC_PROG_CC
AM_PROG_CC_C_O

# Enable POSIX and platform extension preprocessor macros.
AC_USE_SYSTEM_EXTENSIONS

# Add libtool support.
LT_INIT

# Check for compiler version and vendor.
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

# Need C99 and inline support.
AC_PROG_CC_C99
AC_C_INLINE

# Only optimize if allowed, otherwise assume user will set CFLAGS as
# appropriate.
AC_ARG_ENABLE([optimization],
   [AS_HELP_STRING([--enable-optimization],
     [Enable compile time optimization flags for host @<:@default=yes@:>@]
   )],
   [enable_opt="$enableval"],
   [enable_opt="yes"]
)

if test "$enable_opt" = "yes" ; then

   # Add code optimisation flags and tuning to host. This is a funny macro
   # that does not like CFLAGS being already set. Work around that as we have
   # at least set it to "", so it is set.
   ac_test_CFLAGS="no"
   old_CFLAGS="$CFLAGS"
   AX_CC_MAXOPT
   ac_test_CFLAGS="yes"
   CFLAGS="$old_CFLAGS $CFLAGS"

   # Check SSE & AVX support (some overlap with AX_CC_MAXOPT).
   AX_EXT
   if test "$SIMD_FLAGS" != ""; then
       CFLAGS="$CFLAGS $SIMD_FLAGS"
   fi
fi

# Add address sanitizer options to flags, if requested. Only useful for GCC
# version 4.8 and later.
AC_ARG_ENABLE([sanitizer],
   [AS_HELP_STRING([--enable-sanitizer],
     [Enable memory error detection using address sanitizer @<:@default=no@:>@]
   )],
   [enable_san="$enableval"],
   [enable_san="no"]
)

if test "$enable_san" = "yes"; then
   if test "$ax_cv_c_compiler_vendor" = "gnu"; then
      AX_COMPARE_VERSION( $ax_cv_c_compiler_version, [ge], [4.8.0],
                          [enable_san="yes"], [enable_san="no"] )
      if test "$enable_san" = "yes"; then
         CFLAGS="$CFLAGS -fsanitize=address -fno-omit-frame-pointer"
      fi
   else
      AC_MSG_WARN([Compiler does not support address sanitizer option])
   fi
fi

# Autoconf stuff
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_HEADER_STDC

# Check for pthreads.
AX_PTHREAD([LIBS="$PTHREAD_LIBS $LIBS" CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
    CC="$PTHREAD_CC" LDFLAGS="$LDFLAGS $PTHREAD_LIBS $LIBS"]
    AC_DEFINE(HAVE_PTHREAD,1,[Define if you have POSIX threads libraries and header files.]),
    AC_MSG_ERROR([Could not find a working version of
    the pthread library. Make sure you have the library and header files installed
    or use CPPFLAGS and LDFLAGS if the library is installed in a
    non-standard location.]))
    
# Check for OpenMP, is required.
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)
enable_openmp="no"
if test -z "$OPENMP_CFLAGS"; then
   AC_MSG_ERROR(Compiler does not support OpenMP, 1)
else
   CFLAGS="$CFLAGS $OPENMP_CFLAGS"
   enable_openmp="yes"
   AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])
fi

# Check for timing functions needed by cycle.h
AC_C_INLINE
AC_HEADER_TIME
AC_CHECK_HEADERS([sys/time.h c_asm.h intrinsics.h mach/mach_time.h])
AC_CHECK_TYPE([hrtime_t],[AC_DEFINE(HAVE_HRTIME_T, 1, [Define to 1 if hrtime_t is defined in <sys/time.h>])],,[#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif])
AC_CHECK_FUNCS([gethrtime read_real_time time_base_to_time clock_gettime mach_absolute_time])
AC_MSG_CHECKING([for _rtc intrinsic])
rtc_ok=yes
AC_TRY_LINK([#ifdef HAVE_INTRINSICS_H
#include <intrinsics.h>
#endif], [_rtc()], [AC_DEFINE(HAVE__RTC,1,[Define if you have the UNICOS _rtc() intrinsic.])], [rtc_ok=no])
AC_MSG_RESULT($rtc_ok)

# Check for feenableexcept, should be present in C99, also checks for -lm.
AC_CHECK_LIB([m],[feenableexcept])

# Check for git, needed for revision stamps.
AC_PATH_PROG([GIT_CMD], [git])
AC_SUBST([GIT_CMD])

# make the documentation
DX_INIT_DOXYGEN(libquicksched,doc/Doxyfile,doc/)

# .in files.
AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile doc/Makefile doc/Doxyfile])

# Report general configuration.
AC_MSG_RESULT([
   Compiler: $CC
     vendor: $ax_cv_c_compiler_vendor
    version: $ax_cv_c_compiler_version
      flags: $CFLAGS
   OpenMP enabled: $enable_openmp
])

# generate output.
AC_OUTPUT
