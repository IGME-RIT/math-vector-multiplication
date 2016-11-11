/*
Title: Vector Mathematics
File Name: Vector2D.h
Copyright © 2016
Author: Andrew Litfin
Written under the supervision of David I. Schwartz, Ph.D., and
supported by a professional development seed grant from the B. Thomas
Golisano College of Computing & Information Sciences
(https://www.rit.edu/gccis) at the Rochester Institute of Technology.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include <iostream>
#include <math.h>

#include "helpers.h"

// A two-dimensional vector.
// Has an x and y component.
struct Vector2D
{
	float x, y;

	// Gives the zero vector.
	// Equivalent to Vector2D(0, 0)
	Vector2D();

	Vector2D(float x, float y);
};

// Takes a vector v to -v such that v + (-v) = 0.
Vector2D operator-(Vector2D v);

// Componentwise vector addition
Vector2D operator+(Vector2D l, Vector2D r);

// Equivalent to l + (-r)
Vector2D operator-(Vector2D l, Vector2D r);

// Scalar multiplication of a vector and a scalar
Vector2D operator*(float s, Vector2D v);

// Scalar multiplication is commutative
Vector2D operator*(Vector2D v, float s);

// Equivalent to (1/s)*v
Vector2D operator/(Vector2D v, float s);

// Check for exact equality (will return false on vectors that are "very close," but not exact)
// Also, be wary of floating-point error.
bool operator==(Vector2D l, Vector2D r);
// Negation of ==
bool operator!=(Vector2D l, Vector2D r);

// Standard dot product on R^n
float Dot(Vector2D l, Vector2D r);

// Projection of a onto b
Vector2D Project(Vector2D a, Vector2D b);
// Rejection of a from b
Vector2D Reject(Vector2D a, Vector2D b);

// Calculates |v|
float Magnitude(Vector2D v);
// Calculates 1/|v|
float MagInverse(Vector2D v);
// Calculates 1/|v| using FastInvSqrt
float MagFastInv(Vector2D v);
// Calculates |v|^2
// Useful for when lengths need to be compared, as it avoids a square root operation
float MagSquared(Vector2D v);

// Overloads the << operator for use primarily with std::cout
std::ostream& operator<<(std::ostream& os, Vector2D v);
