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

#include "codes_gpsl1ca.h"

#include "io.h"

#include <iostream>

GpsL1caCodes::GpsL1caCodes()
{
    unsigned int N = IoCodes::readCodes("src/codes/gps_l1_ca_codes.txt", cacodes);
    if (N < 32)
    {
        std::cerr << "GPS L1 C/A codes: Unexpected number of spreading codes read: " << N << std::endl;
    }
}

const code_t* GpsL1caCodes::getPrimaryCode(unsigned int svId)
{
    if (svId <= cacodes.size())
    {
        return &cacodes[svId-1];
    }
    else
    {
        return nullptr;
    }
    
}
