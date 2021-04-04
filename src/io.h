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

#include <fstream>
#include <vector>

//! Class for reading sample data from an SX-NSR device
class IoSxnsr
{
public:
    //! Constructor with path to the file to be read
    IoSxnsr(const std::string &filename);
    //! Destructor
    ~IoSxnsr();

    //! Read N samples and return them in a real-valued float vector. Returns the number of samples read
    std::size_t readSamplesRealFloat(const std::size_t N, std::vector<float> &smplBuf);

private:
    //! handle of the opened file
    std::ifstream fileHandle;

    //! allocate in vector smplBuf at least N entries
    template <typename T>
    bool allocateSmplVector(std::vector<T> &smplBuf, const std::size_t N);
};

//! Class to parse spreading codes from files
class IoCodes
{
public:
    //! Reads the spreading codes from a file. Returns number of codes parsed
    static unsigned int readCodes(const std::string &filename, std::vector<code_t> &codes);
};
