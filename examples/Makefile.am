
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

# Automake stuff
AUTOMAKE_OPTIONS=gnu

# Add the source directory and debug to CFLAGS
AM_CFLAGS = -g -Wall -Werror -I../src $(OPENMP_CFLAGS) -DCPU_TPS=2.67e9 -DTIMERS \
    # -fsanitize=address -fno-omit-frame-pointer

AM_LDFLAGS = -lm # -fsanitize=address

# Set-up the library
bin_PROGRAMS = test test_qr test_bh

# Sources for test
test_SOURCES = test.c
test_CFLAGS = $(AM_CFLAGS)
test_LDADD =  ../src/.libs/libquicksched.a

# Sources for test_qr
test_qr_SOURCES = test_qr.c
test_qr_CFLAGS = $(AM_CFLAGS)
test_qr_LDADD =  ../src/.libs/libquicksched.a -llapacke -llapack -lblas

# Sources for test_bh
test_bh_SOURCES = test_bh.c
test_bh_CFLAGS = $(AM_CFLAGS)
test_bh_LDADD =  ../src/.libs/libquicksched.a

