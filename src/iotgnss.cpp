// This file is part of iotgnss.
//
// iotgnss is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// iotgnss is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with iotgnss.  If not, see <http://www.gnu.org/licenses/>.

#include "io.h"

#include <iostream>
#include <vector>

int main()
{
	IoSxnsr testdata("../matlab/snapshotReceiver/testdata/sxnsr_test_data_flight.stream");
	std::vector<float> s;
	auto n = testdata.readSamplesRealFloat(204800, s);
	std::cout << n << " samples read" << std::endl;
	return 0;
}
