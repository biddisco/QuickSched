
# This file is part of Quicksched.
# Coypright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
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

# Add the debug flag to the whole thing
AM_CFLAGS = -g -O3 -Wall -Werror -ffast-math -fstrict-aliasing -ftree-vectorize \
    -funroll-loops $(SIMD_FLAGS) $(OPENMP_CFLAGS) -DTIMERS \
    # -fsanitize=address -fno-omit-frame-pointer

# Assign a "safe" version number
AM_LDFLAGS = -version-info 0:0:0

# Build the libquicksched library
lib_LTLIBRARIES = libquicksched.la
libquicksched_la_SOURCES = qsched.c queue.c

# List required headers
include_HEADERS = atomic.h lock.h queue.h qsched.h task.h res.h error.h

# Private headers.
noinst_HEADERS = cycle.h quicksched.h
