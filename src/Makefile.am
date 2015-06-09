
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
AM_CFLAGS = -DTIMERS

# Assign a "safe" version number
AM_LDFLAGS = -version-info 0:0:0

# The git command, if available.
GIT_CMD = @GIT_CMD@

# Build the libquicksched library
lib_LTLIBRARIES = libquicksched.la
libquicksched_la_SOURCES = qsched.c queue.c

# List required headers
include_HEADERS = atomic.h lock.h queue.h qsched.h task.h res.h error.h

# Private headers.
noinst_HEADERS = cycle.h quicksched.h

# Versioning. If any sources change then update the version.h file with
# the current git revision and package version.
version.h: version.h.in $(AM_SOURCES) $(include_HEADERS) $(noinst_HEADERS)
	if test "X$(GIT_CMD)" != "X"; then \
	    GIT_REVISION=`git describe --abbrev=8  --always --tags --dirty`; \
	    sed -e "s,@PACKAGE_VERSION\@,$(PACKAGE_VERSION)," \
	        -e "s,@GIT_REVISION\@,$${GIT_REVISION}," version.h.in > version.h; \
	fi

#  Make sure version.h is built first.
BUILT_SOURCES = version.h

#  And distribute the built files.
EXTRA_DIST = version.h version.h.in

