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
 *   Doamin integration methods for scalars
 *
 *************************************************************************/

#ifndef DOMAININTEGRATION_H
#define DOMAININTEGRATION_H


#include "pointSet.h"
#include "scalarField.h"

// Integration of boundary values
double bndIntegration(const pointSet& mesh, const scalarField& sf, int bnd);

// Integration of domain values
double bgGridIntegration(const pointSet& mesh, const scalarField& val);





#endif
