
# This file is part of SWIFT.
# Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
#                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

# Add the source directory and debug to CFLAGS
AM_CFLAGS = -I../src -DCPU_TPS=2.67e9 -DTIMERS 

AM_LDFLAGS = -lcblas

# Set-up the library
bin_PROGRAMS = test test_bh test_bh_sorted test_fmm_sorted
if HAVECBLAS
bin_PROGRAMS += test_qr 
endif

# Sources for test
test_SOURCES = test.c
test_CFLAGS = $(AM_CFLAGS)
test_LDADD =  ../src/.libs/libquicksched.a

# Sources for test_qr
test_qr_SOURCES = test_qr.c
test_qr_CFLAGS = $(AM_CFLAGS)
test_qr_LDADD =  ../src/.libs/libquicksched.a -llapacke -llapacke -lblas

# Sources for test_bh
test_bh_SOURCES = test_bh.c
test_bh_CFLAGS = $(AM_CFLAGS)
test_bh_LDADD =  ../src/.libs/libquicksched.a

# Sources for test_bh_sorted
test_bh_sorted_SOURCES = test_bh_sorted.c
test_bh_sorted_CFLAGS = $(AM_CFLAGS)
test_bh_sorted_LDADD =  ../src/.libs/libquicksched.a

# Sources for test_fmm_sorted
test_fmm_sorted_SOURCES = test_fmm_sorted.c
test_fmm_sorted_CFLAGS = $(AM_CFLAGS)
test_fmm_sorted_LDADD =  ../src/.libs/libquicksched.a

