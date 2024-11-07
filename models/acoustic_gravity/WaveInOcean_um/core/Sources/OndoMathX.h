// This file is part of OndomathX, a lightweight C++ template library
// for PDE simulation using spectral finite lelements.
//
// Copyright (C) 2023 Sebastien Imperiale <sebastien.imperiale@inria.fr >
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. A copy of the MPL can be obtain at
// http://mozilla.org/MPL/2.0/.
#pragma once

// Definition of types for indices and real values
#include <array>
#include <complex>
#include <cassert>

#include "omp.h"

namespace OndoMathX
{
	typedef size_t Index;

	typedef double Real;

	typedef std::complex<Real> Complex;

	typedef std::array<Real, 3> RealVector;

	typedef std::array<std::array<Real, 3>, 3> RealMatrix3x3;

	typedef std::array<std::array<Real, 2>, 2> RealMatrix2x2;

	constexpr Real Zero_Machine = 1e-100;
}

#include "Utility/Utility.h"
#include "Field/Field.h"
#include "Mesh/Mesh.h"
#include "FEM/FEM.h"
