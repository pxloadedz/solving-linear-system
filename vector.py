#Most common vector operations for applications of linear algebra.
#Be familiar with the arithmetic of vectors and the geometry underlying vector algebra.
#Next, start exploring how we can use these basic operations to solve different types of problems.

from decimal import Decimal
from math import *

class Vector(object):
	def __init__(self, coordinates):
		try:
			if not coordinates:
				raise ValueError
			for i, s in enumerate(coordinates):
				coordinates[i] = round(Decimal(s), 3)
			self.coordinates = tuple(coordinates)
			self.dimension = len(coordinates)

		except ValueError:
			raise ValueError('The coordinates must not be empty')

		except TypeError:
			raise TypeError('The coordinates must be an iterable')


	def __str__(self):
		return 'Vector: {}'.format(self.coordinates)

	def __eq__(self, v):
		return self.coordinates == v.coordinates

	def minus(self, v):
		new_coordinates = [x-y for x, y in zip(self.coordinates, v.coordinates)]
		return Vector(new_coordinates)

	def plus(self, v):
		new_coordinates = [x+y for x, y in zip(self.coordinates, v.coordinates)]
		return Vector(new_coordinates)

	def times_scale(self, c):
		c = Decimal(c)
		new_coordinates = [c*x for x in self.coordinates]
		return Vector(new_coordinates)

	def magnitude(self):
		return Decimal(sqrt(sum(x**Decimal('2') for x in self.coordinates)))

	def normalized(self):
		try:
			scale = 1/self.magnitude()
			return self.times_scale(scale)

		except ZeroDivisionError:
			raise Exception('Cannot normalize zero vector')

	def dot(self, v):
		return sum(x*y for x, y in zip(self.coordinates, v.coordinates))

	#def angle_rad(self, v):
	#	return acos(self.dot_product(v) / (self.magnitude() * v.magnitude()))

	#def angle_deg(self, v):
	#	return degrees(acos(self.dot_product(v) / (self.magnitude() * v.magnitude())))

	def angle_with(self, v, in_degrees=False):
	# angle = arccos (dot product of the two normalized vectors)
		try:
			u1 = self.normalized()
			u2 = v.normalized()
			angle_in_rad = acos(u1.dot(u2))

			if in_degrees:
				deg_per_rad = Decimal('180') / pi
				return angle_in_rad * deg_per_rad
			else:
				return angle_in_rad

		except Exception as e:
			if str(e) == 'Cannot normalize zero vector':
				raise Exception('Cannot compute an angle with the zero vector')
			else:
				raise e

	#def is_zero(self):
	#	check = []
	#	for x in self.coordinates:
	#		if x == 0:
	#			check.append(x)
	#	if len(check) == self.dimension:
	#		return "zero"
	#	else:
	#		return "not zero"

	def is_zero(self, tolerance=1e-10):
		return self.magnitude() < tolerance

	def parallel_orthogonal(self, v, tolerance=1e-10):
		if self.is_zero() or v.is_zero():
			return "both parallel and orthogonal"
		elif abs(abs(self.dot(v)/(self.magnitude()*v.magnitude())) - 1) < 1e-6:
		#self.dot(v)/(self.magnitude()*v.magnitude()) is 1 or -1
			return "parallel"
		elif abs(self.dot(v)) < tolerance:
			return "orthogonal"
		else:
			return "neither parallel nor orthogonal"

	def component_parallel_to(self, b):
	#returns vector that is the projection of v onto b
		try:
			unit_b = b.normalized()
			return unit_b.times_scale(self.dot(unit_b))

		except Exception as e:
			if str(e) == "Cannot normalize zero vector":
				raise Exception("No unique parallel component")
			else:
				raise e

	def component_orthogonal_to(self, b):
	#returns vector that is a component of v orthogonal to b
		try:
			return self.minus(self.component_parallel_to(b))

		except Exception as e:
			if str(e) == "No unique parallel component":
				raise Exception("No unique orthogonal component")
			else:
				raise e

	def cross(self, w):
	#returns vector that is a cross product of self and w
		a = list(self.coordinates)
		b = list(w.coordinates)
		if self.dimension == 2:
			a.append(0)
			b.append(0)
		x = a[1]*b[2] - b[1]*a[2]
		y = -(a[0]*b[2] - b[0]*a[2])
		z = a[0]*b[1] - b[0]*a[1]
		new_coordinates = [x, y, z]
		return Vector(new_coordinates)

	def area_of_parallelogram(self, w):
		x_prod = self.cross(w)
		return x_prod.magnitude()

	def area_of_triangle(self, w):
		return 0.5*self.area_of_parallelogram(w)

