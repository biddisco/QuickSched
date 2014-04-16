
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
AC_INIT(src/qsched.c)
VERSION="0.1.0"
AC_SUBST(VERSION)
AM_INIT_AUTOMAKE(QUICKSCHED,VERSION)
AC_GNU_SOURCE
AC_PROG_LIBTOOL

# my own macro collection
AC_CONFIG_MACRO_DIR([m4])

# generate header file
AM_CONFIG_HEADER(config.h)

# compiler settings
#CFLAGS="-Wall $(CFLAGS)"

# find and test the compiler
AM_PROG_CC_C_O
AC_PROG_CC_C99
AC_LANG_C
# AX_CC_MAXOPT
AX_FUNC_POSIX_MEMALIGN
AX_GCC_ARCHFLAG([no])
AX_EXT

# autoconf stuff
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_HEADER_STDC

# Check for pthreads
ACX_PTHREAD([LIBS="$PTHREAD_LIBS $LIBS" CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
    CC="$PTHREAD_CC" LDFLAGS="$PTHREAD_LIBS $LIBS"
    AC_DEFINE(HAVE_PTHREAD,1,[Define if you have POSIX threads libraries and header files.])],
    AC_MSG_ERROR([Could not find a working version of
    the pthread library. Make sure you have the library and header files installed
    or use CPPFLAGS and LDFLAGS if the library is installed in a
    non-standard location.]))
    
# Check for OpenMP
AX_OPENMP
AC_SUBST(OPENMP_CFLAGS)
if test -z "${OPENMP_CFLAGS}"; then
   echo $OPENMP_CFLAGS
   AC_MSG_ERROR(Compiler does not support OpenMP, 1)
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


# make the documentation
DX_INIT_DOXYGEN(libchebfun,doc/Doxyfile,doc/)

# generate the Makefiles
AC_OUTPUT(Makefile src/Makefile examples/Makefile doc/Makefile)
