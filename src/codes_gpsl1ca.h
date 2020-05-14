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

#pragma once

#include "codes.h"

//! Class for obtaining GPS L1 C/A spreading codes
class GpsL1caCodes : public Codes
{
public:
    //! Constructor
    GpsL1caCodes();
    //! Destructor
    ~GpsL1caCodes(){};

    //! Returns the spreading code for the given satellite
    const code_t* getPrimaryCode(unsigned int svId);

private:
    //! spreading codes
    std::vector<code_t> cacodes;
};
