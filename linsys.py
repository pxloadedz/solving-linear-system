from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):
        p_row1 = self.planes[row1]
        p_row2 = self.planes[row2]
        self.planes[row1] = p_row2
        self.planes[row2] = p_row1

    def multiply_coefficient_and_row(self, coefficient, row):
        p_row = self.planes[row]
        p_row.normal_vector = p_row.normal_vector.times_scale(coefficient)
        p_row.constant_term = coefficient * p_row.constant_term
        p_row.set_basepoint()
        #need to set basepoint again because normal vector and constant term has changed
        #it did not automatically update basepoint because did not go through init fn

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):

        normal_vector_to_add = self.planes[row_to_add].normal_vector
        normal_vector_to_add = normal_vector_to_add.times_scale(coefficient)
        self.planes[row_to_be_added_to].normal_vector = self.planes[row_to_be_added_to].normal_vector.plus(normal_vector_to_add)

        constant_to_add = self.planes[row_to_add].constant_term
        constant_to_add *= coefficient
        self.planes[row_to_be_added_to].constant_term += constant_to_add

        self.planes[row_to_be_added_to].set_basepoint()
        #can also set self[row_to_be_added] = new Plane()

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector.coordinates)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices


    def compute_triangular_form(self):
        #compute the triangular form in gaussian elimination that has the same
        #solution set as the original system

        system = deepcopy(self) #so we don't mess up the original system
        
        num_equations = len(self)
        num_variables = self.dimension
        i = 0
        while i < num_equations:
            j = 0
            while j < num_variables:
                eq_n = system[i].normal_vector.coordinates
                c = eq_n[j] #coefficient of var j in the row i
                if MyDecimal(c).is_near_zero():
                    row_to_swap = system.found_var_underneath(j, i)
                    if row_to_swap != 999:
                        system.swap_rows(i, row_to_swap)
                    else:
                        j += 1
                else:
                    break
            if j < num_variables:
                system.clear_var_underneath(j, i)
            i += 1
        return system

    def found_var_underneath(self, var, row):
        num_equations = len(self)
        next_row = row + 1
        while next_row < num_equations:
            eq_n = self[next_row].normal_vector.coordinates
            c = eq_n[var]
            if not MyDecimal(c).is_near_zero():
                return next_row
            else:
                next_row += 1
        return 999

    def clear_var_underneath(self, var, row):
        num_equations = len(self)
        while row < num_equations:
            row_to_clear = self.found_var_underneath(var, row)
            if row_to_clear != 999:
                eq_n = self[row].normal_vector.coordinates
                eq_n_c = eq_n[var]
                eq_clear_n = self[row_to_clear].normal_vector.coordinates
                eq_clear_n_c = eq_clear_n[var]
                coefficient = -(eq_clear_n_c/eq_n_c)
                self.add_multiple_times_row_to_row(coefficient, row, row_to_clear)
            else:
                row = row_to_clear

    def clear_var_above(self, var, row):
        above_row = row - 1
        while above_row > -1:
            row_to_clear = self[above_row]
            eq_clear_n = row_to_clear.normal_vector.coordinates
            if not MyDecimal(eq_clear_n[var]).is_near_zero():
                coefficient = -(eq_clear_n[var])
                self.add_multiple_times_row_to_row(coefficient, row, above_row)
            above_row -= 1

    def compute_rref(self):
        #compute reduced row echelon form
        tf = self.compute_triangular_form()
        pivot_var = tf.indices_of_first_nonzero_terms_in_each_row()
        for i in range(len(pivot_var)):
            eq_n = tf[i].normal_vector.coordinates
            pivot_index = pivot_var[i]
            if pivot_index != -1:
                c = eq_n[pivot_index]
                if c != 1:
                    tf.multiply_coefficient_and_row(1/c, i) 
                tf.clear_var_above(pivot_index, i)
        return tf


    def compute_ge_solution_with_parametrization(self):
        #compute unique sol, no sol, or infinitely many sol
        rref = self.compute_rref()
        zero_nv = Vector(['0']*self.dimension)
        last_eq = rref[-1]
        if last_eq.normal_vector == zero_nv and not MyDecimal(last_eq.constant_term).is_near_zero():
            return "no solutions"

        pivot_vars = rref.indices_of_first_nonzero_terms_in_each_row()
        #calc number of non-zero pivots
        nz_pivot = 0
        for i in pivot_vars:
            if i != -1:
                nz_pivot += 1
        if nz_pivot == self.dimension:
            solutions = []
            for i in range(len(pivot_vars)):
                solutions.append(rref[i].constant_term)
            return solutions
        else:
            return rref.compute_parametrization()

    def compute_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        free_variable_indices = set(range(num_variables)) - set(pivot_indices)

        direction_vectors = []

        for free_var in free_variable_indices:
            vector_coords = [0]*num_variables
            vector_coords[free_var] = 1
            for i, p in enumerate(self.planes):
                pivot_var = pivot_indices[i]
                if pivot_var == -1:
                    break
                vector_coords[pivot_var] = -p.normal_vector.coordinates[free_var]
            direction_vectors.append(vector_coords)

        basepoint_coords = [0]*num_variables
        for i, p in enumerate(self.planes):
            pivot_var = pivot_indices[i]
            if pivot_var == -1:
                break
            basepoint_coords[pivot_var] = p.constant_term

        return 'Basepoint = {}; Direction Vectors = {}'.format(basepoint_coords, direction_vectors)

    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
    #test if new item x has the same dimension as self
    #if True, then set x as the plane at index i in the self linear system
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


p1 = Plane(normal_vector=Vector(['0.786','0.786', '0.588']), constant_term='-0.714')
p2 = Plane(normal_vector=Vector(['-0.138','-0.138','0.244']), constant_term='0.319')
s = LinearSystem([p1, p2])
print(s.compute_ge_solution_with_parametrization())

p1 = Plane(normal_vector=Vector(['8.631','5.112','-1.816']), constant_term='-5.113')
p2 = Plane(normal_vector=Vector(['4.315','11.132','-5.27']), constant_term='-6.775')
p3 = Plane(normal_vector=Vector(['-2.158','3.01','-1.727']), constant_term='-0.831')
s = LinearSystem([p1, p2, p3])
print(s.compute_ge_solution_with_parametrization())

p1 = Plane(normal_vector=Vector(['0.935','1.76','-9.365']), constant_term='-9.955')
p2 = Plane(normal_vector=Vector(['0.187','0.352','-1.873']), constant_term='-1.991')
p3 = Plane(normal_vector=Vector(['0.374','0.704','-3.746']), constant_term='-3.982')
p4 = Plane(normal_vector=Vector(['-0.561','-1.056','5.619']), constant_term='5.973')
s = LinearSystem([p1, p2, p3, p4])
print(s.compute_ge_solution_with_parametrization())

