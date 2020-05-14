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

IoSxnsr::IoSxnsr(const std::string &filename)
{
    fileHandle = std::ifstream(filename, std::ifstream::in | std::ifstream::binary);
}

IoSxnsr::~IoSxnsr()
{
    if (fileHandle.is_open())
    {
        // close the file handle
        fileHandle.close();
    }
}

std::size_t IoSxnsr::readSamplesRealFloat(const std::size_t N, std::vector<float> &smplBuf)
{
    std::size_t numSmplRead = 0;

    if (!fileHandle.fail() && fileHandle.is_open())
    {
        // pre-allocate memory for the samples (if not yet done already)
        smplBuf.clear();
        if (allocateSmplVector(smplBuf, N))
        {
            std::istreambuf_iterator<char> buf(fileHandle);

            unsigned char b = 0;
            while (!fileHandle.eof() && numSmplRead < N)
            {
                if ((numSmplRead % 4) == 0)
                {
                    b = *buf++;
                }

                // interpret the two leas significant bits as a 2's complement integer
                unsigned char bm = b & 0x3;
                float s = ( bm == 0x3 ? -1.0f :
                            bm == 0x2 ? -2.0f :
                            bm == 0x1 ?  1.0f :
                                         0.0f );
                smplBuf.emplace_back(s);
                b >>= 2;
                ++numSmplRead;
            }
        }
    }

    return numSmplRead;
}

template <typename T>
bool IoSxnsr::allocateSmplVector(std::vector<T> &smplBuf, const std::size_t N)
{
    if (smplBuf.capacity() < N)
    {
        smplBuf.reserve(N);
    }

    return smplBuf.capacity() >= N;
}

unsigned int IoCodes::readCodes(const std::string &filename, std::vector<std::vector<bool>> &codes)
{
    std::ifstream f(filename, std::ios::in);

    unsigned int numCodesParsed = 0;

    if (f && !f.fail() && f.is_open())
    {
        codes.clear();

        // parse the input codes file line by line
        std::string line;
        while (f >> line)
        {
            // reserve space for the expected number of chips
            std::vector<bool> code(line.length());
            auto itCode = code.begin();
            for (const auto c : line)
            {
                if ('0' <= c && c <= '1')
                {
                    *itCode++ = (c - '0' == 1);
                }
            }

            codes.emplace_back(code);
            ++numCodesParsed;
        }
    }

    return numCodesParsed;
}
