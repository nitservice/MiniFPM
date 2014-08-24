/*************************************************************************
 *
 *  This file is part of MiniFPM
 *  Copyright (C) 2014 Jan Marburger
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Description:
 *   Matrix class using SparseLib++
 *
 *************************************************************************/


#ifndef MATRIXIML_H
#define MATRIXIML_H

#include "mfpmMatrix.h"
#include "Array.h"

int solve3x3(
    const mfpmMatrix& a11,
    const mfpmMatrix& a12,
    const mfpmMatrix& a13,
    const mfpmMatrix& a21,
    const mfpmMatrix& a22,
    const mfpmMatrix& a23,
    const mfpmMatrix& a31,
    const mfpmMatrix& a32,
    const mfpmMatrix& a33,
    Array& res1,
    Array& res2,
    Array& res3
);


int solve2x2(
    const mfpmMatrix& a11,
    const mfpmMatrix& a12,
    const mfpmMatrix& a21,
    const mfpmMatrix& a22,
    Array& res1,
    Array& res2
);

int solve1x1(
    const mfpmMatrix& A,
    Array& res1
);


#endif
