# This file is part of iotgnss.
#
# iotgnss is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# iotgnss is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with iotgnss.  If not, see <http:#www.gnu.org/licenses/>.

CC=g++
CFLAGS=-Wall -g -std=c++11

all: iotgnss

iotgnss: iotgnss.o io.o
	$(CC) $(CFLAGS) -o iotgnss iotgnss.o io.o

iotgnss.o: iotgnss.cpp
	$(CC) $(CFLAGS) -c iotgnss.cpp

io.o: io.cpp
	$(CC) $(CFLAGS) -c io.cpp

clean:
	rm -f iotgnss iotgnss.o io.o
