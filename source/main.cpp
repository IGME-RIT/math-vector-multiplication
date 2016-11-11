/*
Title: Intermediate Vector Math
File Name: main.cpp
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

// The bread and butter of games programming is vector mathematics and linear algebra.
// The vast majority of the mathematics used in games falls under this category,
//  though it is not the only math used (e.g. discrete math, calculus).
// This tutorial series will take you through the basics of vector math.
// Future tutorials in this series will detail other aspects of linear algebra, particularly matrices.
// The exposition follows that of "Foundations of Game Engine Development" (Volume 1),
//  by Eric Lengyel.
// This file provides exposition, definitions, and explanations, and all other files implement vector classes
//  as you would see them in most game engines. The implementations are modeled after those of Eric Lengyel
//  in FGED, Volume 1 and the Tombstone Engine, though they are my own.
// Only Vector2D.h is heavily annotated, the others are mostly identical.

// This tutorial will explain some intermediate vector math concepts.
// We will cover
//  1) How to find the Magnitude (also called the length, or norm) of a vector
//  2) How to calculate the dot product of two vectors and its geometric interpretation
//  3) How to calculate the cross product of two 3D vectors and its geometric interpretation
//  4) How to find and use angles between vectors

#include "../header/helpers.h"
#include "../header/Vector2D.h"
#include "../header/Vector3D.h"
#include "../header/Vector4D.h"

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <ctime>

int main()
{
	// Required for the random functions in helpers to work.
	srand(static_cast<unsigned>(time(0)));

	// Magnitude (Norm)
	// ----------------
	{
		// The magnitude of a vector is its length (remember how vectors have both length and direction).
		// For our purposes, it can be calculated just by using the Pythagorean theorem.
		// For example, if a = (1, 1), then its length is |a| = sqrt(1^2 + 1^2) = sqrt(1 + 1) = sqrt(2).
		if (Magnitude(Vector2D(1, 1)) == sqrtf(2))
		{
			std::cout << "|a| = sqrt(2)\n";
		}

		// In 3D and higher dimensional space, the magnitude of a vector can be found by extending the Pythagorean Theorem.
		// Say for example we want to find the magnitude of an arbitrary vector v=(x, y, z).
		// The vector's shadow on the xy plane has distance sqrt(x^2 + y^2).
		// Then the right triangle with base length sqrt(x^2 + y^2) and height z has hypotenuse of length
		//  sqrt(sqrt(x^2 + y^2)^2 + z^2) = sqrt(x^2 + y^2 + z^2).
		// This can be extended any finite number of times to get the general formula for magnitude of a Euclidean vector:
		//  Let v=(v_1, v_2, ..., v_n) \in R^n. Then |v| = sqrt(Sum_{i=1}^{n}(v_i^2))

		// Often times we are less concerned with the magnitude of a vector and more with its direction.
		// We can divide a vector by its magnitude to get a "unit vector," or a "normalized vector."

		Vector2D nonUnitLen = { 1, 1 };
		Vector2D unitLen = nonUnitLen / Magnitude(nonUnitLen);
		std::cout << "Vector unitLen = " << unitLen << " has magnitude " << Magnitude(unitLen) << '\n';

		// A common way to speed up normalization that does not use sqrt is a bit of amazing code called FastInvSqrt.
		// It takes a number x and returns 1/sqrt(x) using black magic and floating point bit-level hackery.
		// However, it is not quite as accurate.
		Vector2D almostUnitLen = MagFastInv(nonUnitLen) * nonUnitLen;
		std::cout << "Vector almostUnitLen = " << almostUnitLen << " has magnitude " << Magnitude(almostUnitLen) << '\n';
		std::cout << "This is an absolute error of " << 1 - Magnitude(almostUnitLen) << ", or a relative error of "
			<< (1 - Magnitude(almostUnitLen)) / 100 << "%\n";
	}

	// Dot (Inner) Product
	// -------------------
	{
		// The dot product is a binary operator that takes two vectors from R^n and returns a real number.
		// Intuitively, it represents how much the vectors "point in the same direction."
		// In a more exact sense, the dot product of two vectors in Euclidean space is
		//   Dot(x, y) = |x||y|cos(theta), where theta is the angle between the two vectors.
		// In writing, the dot product is represented by a single dot between the vectors (eg x . y),
		//  and when using the general inner product, as angled brackets (eg <x, y>).
		// For example, say x = (1, 0) and y = (0, 1)
		Vector2D x = { 1, 0 };
		Vector2D y = { 0, 1 };
		// The angle between these vectors is clearly 90 degrees. Then
		std::cout << "Dot(x, y) = " << Dot(x, y) << '\n';
		// because cos(90*) = 0.

		// There is a relation between the dot product and magnitude, namely for all vectors x, |x|^2 = Dot(x, x),
		//  since the angle between a vector and itself is zero, and hence cos(0) = 1.
		// Indeed, if you look at the code for the Magnitude function, that is how the magnitude is calculated, as |x| = sqrt(Dot(x, x))

		// Now if you examine the code for Dot(x, y), you will see nowhere do I compute magnitudes or cosines.
		// This is because of an alternate, equivalent, way to compute the dot product that arises from the law of cosines.
		// The law of cosines states that for any triangle, c^2 = a^2 + b^2 - 2*a*b*cos(C), where C is the angle opposite side c.
		// This can be reinterpreted in the context of vectors.  Say u and v are vectors in R^n.
		// Then |u - v|^2 = |u|^2 + |v|^2 - 2|u||v|cos(theta)
		//                = |u|^2 + |v|^2 - 2*Dot(u, v),
		// hence Dot(u, v) = 1/2 (|u|^2 + |v|^2 - |u - v|^2)
		//                 = 1/2 Sum_{i=1}^{n}(u_i^2 + v_i^2 - (u_i - v_i)^2)
		//                 = 1/2 Sum_{i=1}^{n}(u_i^2 + v_i^2 - u_i^2 + 2*u_i*v_i - v_i^2)
		//                 = 1/2 Sum_{i=1}^{n}(2*u_i*v_i)
		//                 = Sum_{i=1}^{n}(u_i*v_i)
		// Which, in 2 dimensions, means Dot(u, v) = u_1*v_1 + u_2*v_2, and similar for higher dimensions.
		// This means that rather than calculating square roots and cosines for dot products,
		//  instead we can just add and multiply the components.
		Vector2D u = { 4, -4 };
		Vector2D v = { 0, -1 };
		// The angle between u and v is 45 degrees, as can be seen by graphing the vectors.
		// cos(45)=1/sqrt(2), so we can manually place that value in.
		if (Dot(u, v) - Magnitude(u) * Magnitude(v) / sqrtf(2.0f) < FLT_EPSILON)
		{
			std::cout << "The formulations are equivalent.\n";
		}
		// Now you may ask "why manually place the sqrt in instead of just letting the computer calculate cos(45)?" Two reasons:
		//  1) Trigonometric functions are, in general, slow. Even slower than square roots. In games, speed is your primary concern.
		//  2) This is one of those cases where floating point calculations really makes a difference. See the following:
		float veryCloseToZero = Dot(u, v) - Magnitude(u) * Magnitude(v) / sqrtf(2.0f);
		// and compare it to
		float lessCloseToZero = Dot(u, v) - Magnitude(u) * Magnitude(v) * cosf(3.14159265f / 4.0f);
		// The calculations are equivalent, and yet
		std::cout << "veryCloseToZero: " << veryCloseToZero << " lessCloseToZero: " << lessCloseToZero << "\n";
		// This goes to show that for example, when you know an angle of rotation beforehand (that is, it is hard-coded rather than procedural),
		//  it is worth it to pre-compute the values and use the computed values instead.
	}

	// Cross Product
	// -------------
	{
		// The cross product of two vectors returns a vector that is perpendicular to both input vectors (that is, the dot product is zero)
		//  and has magnitude equal to the area of the parallelogram enclosed by the two input vectors.
		// The cross product is special in that it is only well-defined in three dimensions.
		// (Technically it is also well-defined in seven dimensions, but that is irrelevant to our purposes.)

		// Since we know from geometry that the area of a parallelogram is just base * height, we have that
		// |a*b| = |a|*|b|sin(theta), where theta is the angle between a and b.

		/*

		Cross(a, b)
		^
		|
		|
		|    b
		|    ^- - - - - - - - /
		|   /
		|  /     Area =     /
		| /  |Cross(a, b)|
		|/                /
		|----------------> a

		*/

		// There are two forumlations that we will be concerned with.
		// First is the explicit form, which can be found in the source code for any of the Vector structs.
		// The second is the skew-symmetric form, which will not be covered until after matrices have been introduced.

		// To illustrate that the cross product is perpendicular to any two input vectors, a and b are random vectors in the positive unit cube.
		Vector3D a = { randFloat(0, 1), randFloat(0, 1), randFloat(0, 1) };
		Vector3D b = { randFloat(0, 1), randFloat(0, 1), randFloat(0, 1) };
		Vector3D c = Cross(a, b);
		std::cout << "a = " << a << "\nb = " << b << "\nc = " << c << '\n';
		if (Dot(a, c) < FLT_EPSILON && Dot(b, c) < FLT_EPSILON)
		{
			std::cout << "a and b are orthogonal to c.\n";
		}

		// The cross product distributes over addition as follows:
		// a*(b + c) = a*b + a*c
		if (Magnitude(Cross(a, b + c) - (Cross(a, b) + Cross(a, c))) < FLT_EPSILON)
		{
			std::cout << "The cross product distributes over addition.\n";
		}

		// The following three properties are not critical to understanding the cross product, but nevertheless interesting.

		// Vector triple identity
		// a*(b*c) = (a.c)b - (a.b)c
		// The proof can be easily brute forced by expanding both sides.
		// There is a slightly easier way to prove it, however, it uses the skew-symmetric form of the cross product, which requires matrices.
		if (Magnitude(Cross(a, Cross(b, c)) - (Dot(a, c)*b - Dot(a, b)*c)) < FLT_EPSILON)
		{
			std::cout << "The cross product satisfies the vector triple identity.\n";
		}

		// Lagrange identity
		// |a*b|^2 = |a|^2*|b|^2 - (a.b)^2
		// This is a step along the way to proving the above relation between the magnitude of the cross product and the area of a parallelogram.
		// First, |a*b|^2 = (a_y*b_z - a_z*b_y)^2 + (a_z*b_x - a_x*b_z)^2 + (a_x*b_y - a_y*b_x)^2
		//                = (a_y*b_z)^2 - 2*a_y*a_z*b_y*b_z + (a_z*b_y)^2
		//                + (a_z*b_x)^2 - 2*a_x*a_z*b_x*b_z + (a_x*b_z)^2
		//                + (a_x*b_y)^2 - 2*a_x*a_y*b_x*b_y + (a_y*b_x)^2
		// Adding (a_x*b_x)^2 + (a_y*b_y)^2 + (a_z*b_z)^2 to the front and subtracting it from the back yields
		//                = (a_x^2 + a_y^2 + a_z^2)*(b_x^2 + b_y^2 + b_z^2) - (a_x*b_x + a_y*b_y + a_z*b_z)^2
		//                = |a|^2*|b|^2 - (a.b)^2
		// Since (a.b)^2 = |a|^2*|b|^2*cos^2(theta), then
		//                = |a|^2*|b|^2 - |a|^2*|b|^2*cos^2(theta)
		//                = |a|^2*|b|^2*(1 - cos^2(theta))
		//                = |a|^2*|b|^2*sin^2(theta) [By the Pythagorean Theorem, cos^2 + sin^2 = 1]
		// Taking square roots yields |a*b| = |a|*|b|*sin(theta).
		// If |a| is considered the length of the base of the parallelogram, then |b|*sin(theta) is the height, hence |a*b| is the area.
		if (MagSquared(Cross(a, b)) - (MagSquared(a)*MagSquared(b) - powf(Dot(a, b), 2.0f)) < FLT_EPSILON)
		{
			std::cout << "The Lagrange identity holds.\n";
		}

		// The cross product also satisfies what is known as the Jacobi Identity, which is only really useful if you're interested in higher algebras.
		// It is included for completeness.
		// The Jacobi identity is satisfied if, for a binary operation * together with addition, we have that for all a, b, and c
		//  a*(b*c) + c*(a*b) + b*(c*a) = 0.
		// Note that * need not be associative.
		// The cross product satisfies the Jacobi Identity, as can be seen by expanding each cross product.
		// Again, because of floating point inaccuracies, we must slightly alter how we check to make sure that it indeed is satisfied.
		if (Magnitude(Cross(a, Cross(b, c)) + Cross(c, Cross(a, b)) + Cross(b, Cross(c, a))) < FLT_EPSILON)
		{
			std::cout << "The cross product satisfies the Jacobi Identity.\n";
		}
	}

	// Angles
	// ------
	{
		// Recall that a.b was initially defined as |a||b|cos(theta).
		// This means that, since we have an alternate, but equivlant, way of finding a.b as the sum of products,
		//  then we can actually solve the equation for theta.
		// First, though, suppose that either a or b is the zero vector, such that |a|=0 or |b|=0. Then the equation is trivially true for any value of theta.
		// This is what we call a "degenerate case," and we are not interested in such cases. A consequence is that we call the zero vector "orthogonal"
		//  to all other vectors, even though it has lost all sense of geometric meaning.

		// So, suppose then that neither |a|=0 nor |b|=0.  Then theta = acos((a.b)/(|a||b|)) = acos(a^.b^),
		//  which has a unique solution in [0, pi].

		// A common use case for these angles is determining if a vector is in the positive or negative half-space separated by a plane normal to another vector.
		// Or, in plain English, whether one vector points in the same direction as another.

		// For example, consider the three vectors
		Vector2D a = { 1, 0 };
		Vector2D b = { 1, 1 };
		Vector2D c = { -1, 1 };
		// I would highly recommend you grab a piece of paper and pencil and follow along.

		// Vectors a and b both lie on the right side of the y axis.  We can confirm this by finding theta:
		float thetaab = acosf(Dot(a, b) / (Magnitude(a)*Magnitude(b)));
		std::cout << "Theta is " << thetaab << ", or " << thetaab*(180.0f / 3.14159265f) << " degrees.\n";
		if (thetaab < (3.14159265f / 2.0f))
		{
			std::cout << "Since theta is less than 90 degrees, then a and b point roughly \"in the same direction.\"\n";
		}

		// However, on the other hand,
		float thetaac = acosf(Dot(a, c) / (Magnitude(a)*Magnitude(c)));
		if (thetaac > (3.14159265f / 2.0f))
		{
			std::cout << "The angle between a and c is " << thetaac*(180.0f / 3.14159265f) << " degrees.\n"
				<< "Since it's greater than 90 degrees, a and c point roughly \"in the opposite direction.\"\n";
		}


		// Now, if all we want to know is that the vectors are in the same half-space, and it's not important by how much,
		//  we can utilize a quirk of the cosine function to speed up the calculation.
		// Say for example you have a character in 3D space with some orientation
		float yawAngle = randFloat(0, 2 * 3.14159265f);
		Vector3D forward = { cosf(yawAngle), sinf(yawAngle), 0 };
		// and
		Vector3D right = { sinf(yawAngle), -cosf(yawAngle), 0};
		// This of course implies
		Vector3D up = { 0, 0, 1 };
		// This is a right-handed basis for the character.

		std::cout << "The character has yaw angle " << yawAngle*(180.0f / 3.14159265f) << " degrees, giving an orientation basis of\n"
			<< "forward: " << forward << "\nright: " << right << "\nand up " << up << "\n";

		// Now say a projectile is heading towards our character with global velocity
		Vector3D projectileVelocity = 10 * Vector3D(randFloat(-1, 1), randFloat(-1, 1), randFloat(-1, 1));
		while (projectileVelocity == Vector3D(0, 0, 0))
		{
			// Do this to avoid any weirdness that could come about if projectileVelocity is exactly zero.
			projectileVelocity = 10 * Vector3D(randFloat(-1, 1), randFloat(-1, 1), randFloat(-1, 1));
		}

		std::cout << "The projectile has velocity " << projectileVelocity << "\n";
		
		// We want to know, for whatever reason, if the projectile is coming at our character from their right or left side.
		// This amounts to answering the question is theta in [-90, 90] or in [90, 180] U [-90, -180]?
		// However, look at the behavior of the cosine function on those intervals.
		// On [-90, 90], cosine is positive, and on [90, 180] U [-90, -180], it is negative.
		// So, we can determine which interval theta lies in by simply checking if cos(theta) is positive or negative.
		// But since cos(theta) = a.b / |a||b| implies |a||b|cos(theta) = a.b and magnitudes are always positive or zero,
		//  that means that cos(theta) is positive exactly when a.b is positive and cos(theta) is negative exactly when a.b is negative.
		// So we can determine if the projectile is coming at our character from the left or right simply by checking
		//  if the dot product of the velocity with the character's right vector is negative or positive, respectively!

		float dot = Dot(right, projectileVelocity);
		if (dot < 0)
		{
			std::cout << "The projectile is coming at the character from their left side.\n";
		}
		else if (dot > 0)
		{
			std::cout << "The projectile is coming at the character from their right side.\n";
		}
		else
		{
			std::cout << "The projectile is coming straight at the character from neither their left nor right side.\n";
		}

	}

	std::cout << "Press Enter to continue . . . ";
	std::cin.get();
	return 0;
}
