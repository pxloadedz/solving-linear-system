from decimal import Decimal, getcontext

from vector import Vector 

getcontext().prec = 30

class Line(object):

	NO_NONZERO_ELTS_FOUND_MSG = "No nonzero elements found"

	def __init__(self, normal_vector=None, constant_term=None):
	#can easily find the direction vector of the line from the normal vector if need to
		self.dimension = 2

		if not normal_vector:
			all_zeros = [0]*self.dimension
			normal_vector = Vector(all_zeros)
		self.normal_vector = normal_vector

		if not constant_term:
			constant_term = Decimal('0')
		self.constant_term = Decimal(constant_term)

		self.set_basepoint()


	def set_basepoint(self):
		try:
			n = self.normal_vector.coordinates
			c = self.constant_term
			basepoint_coords = [0]*self.dimension

			# if line is y=2 then basepoint coord is [0,2] normal vector is [0, 1], constant = 2

			initial_index = Line.first_nonzero_index(n)
			initial_coefficient = n[initial_index]

			basepoint_coords[initial_index] = Decimal(c)/Decimal(initial_coefficient)
			self.basepoint = Vector(basepoint_coords)

		except Exception as e:
			if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
				self.basepoint = None
			else:
				raise e


	def __str__(self):

		num_decimal_places = 3

		def write_coefficient(coefficient, is_initial_term=False):
			coefficient = round(coefficient, num_decimal_places)
			if coefficient % 1 == 0:
				coefficient = int(coefficient)

			output = ''

			if coefficient < 0:
				output += '-'
			if coefficient > 0 and not is_initial_term:
				output += '+'

			if not is_initial_term:
				output += ' '

			if abs(coefficient) != 1:
				output += '{}'.format(abs(coefficient))

			return output

		n = self.normal_vector.coordinates

		try:
			initial_index = Line.first_nonzero_index(n)
			terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
					 for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
			output = ' '.join(terms)

		except Exception as e:
			if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
				output = '0'
			else:
				raise e

		constant = round(self.constant_term, num_decimal_places)
		if constant % 1 == 0:
			constant = int(constant)
		output += ' = {}'.format(constant)

		return output

	@staticmethod
	def first_nonzero_index(iterable):
		for k, item in enumerate(iterable):
			if not MyDecimal(item).is_near_zero():
				return k
		raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)


	def is_parallel(self, l):
		n1 = self.normal_vector
		n2 = l.normal_vector
		return n1.parallel_orthogonal(n2) in ["parallel", "both parallel and orthogonal"]

	def __eq__(self, l):
	# two parallel lines are equal if the vector connecting one point on each line is
	# orthogonal to the lines' normal vectors

		#if both normal vectors are zero and their constant terms are the same
		#then the lines are equal
		if self.normal_vector.is_zero():
			if not ell.normal_vector_is_zero():
				return False
			else:
				diff = self.constant_term - l.constant_term
				return MyDecimal(diff).is_near_zero()
		elif l.normal_vector.is_zero():
			return False

		if not self.is_parallel(l):
			return False

		v = self.basepoint.minus(l.basepoint)
		return v.parallel_orthogonal(self.normal_vector) and v.parallel_orthogonal(l.normal_vector) \
		in ["orthogonal", "both parallel and orthogonal"]
		#checking if basepoint difference is orthogonal to both normal vectors is redundant
		#because we already know that the lines are parallel

	def compute_intersection(self, l):
	#unless self and l are parallel, they will have a single intersection point
		if self.is_parallel(l):
			if self == l:
				return "infinite intersection" #or return self
			else:
				return "no intersection"
		n1 = list(self.normal_vector.coordinates)
		n2 = list(l.normal_vector.coordinates)
		k1 = self.constant_term
		k2 = l.constant_term
		x = (n2[1]*k1 - n1[1]*k2)/(n1[0]*n2[1]-n1[1]*n2[0])
		y = (-(n2[0])*k1 + n1[0]*k2)/(n1[0]*n2[1]-n1[1]*n2[0])
		return [x, y]


class MyDecimal(Decimal):
	def is_near_zero(self, eps=1e-10):
		return abs(self) < eps

l1 = Line(Vector([4.046, 2.836]), 1.21)
l2 = Line(Vector([10.115, 7.09]), 3.025)

l3 = Line(Vector([7.204, 3.182]), 8.68)
l4 = Line(Vector([8.172, 4.114]), 9.883)

l5 = Line(Vector([1.182, 5.562]), 6.744)
l6 = Line(Vector([1.773, 8.343]), 9.525)

print(l1.compute_intersection(l2))
print(l3.compute_intersection(l4))
print(l5.compute_intersection(l6))