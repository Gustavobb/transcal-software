import numpy as np

class structure:

    def __init__(self, elements, knot_numbers, restriction, elements_number):

        self.elements = elements
        self.knot_numbers = knot_numbers
        self.elements_number = elements_number
        self.restriction = restriction
        
        answ = self.calculate_rigity_matrix()
        self.rigity_matrix_no_outline, self.rigity_matrix = answ[0], answ[1]
        self.force_matrix = self.calculate_force_matrix()
        
        # jacobi ou gauss seidel
        answ = self.methods(10000, 1e-6, "gauss seidel")
        self.u_vector, self.err = answ[0], answ[1]

        answ = self.calculate_deformation()
        self.deformation, self.tension = answ[0], answ[1]

        self.internal_force = self.calculate_internal_force()
        self.reaction = self.calculate_reaction()

    def calculate_rigity_matrix(self):

        rigity_matrix = np.zeros((self.knot_numbers*2, self.knot_numbers*2))

        for element in self.elements:
            rigity_matrix[element.dof[0][0]:element.dof[0][1] + 1, element.dof[0][0]:element.dof[0][1] + 1] += element.rigity_matrix[0:2, 0:2]
            rigity_matrix[element.dof[0][0]:element.dof[0][1] + 1, element.dof[1][0]:element.dof[1][1] + 1] += element.rigity_matrix[0:2, 2:4]
            rigity_matrix[element.dof[1][0]:element.dof[1][1] + 1, element.dof[0][0]:element.dof[0][1] + 1] += element.rigity_matrix[2:4, 0:2]
            rigity_matrix[element.dof[1][0]:element.dof[1][1] + 1, element.dof[1][0]:element.dof[1][1] + 1] += element.rigity_matrix[2:4, 2:4]

        return rigity_matrix, self.apply_outline_condition(rigity_matrix)

    def calculate_force_matrix(self):

        force_matrix = np.zeros((self.knot_numbers*2, 1))
        
        for i in range(self.elements_number - 1):
            force_matrix[self.elements[i].dof[0][0]] = float(self.elements[i].knots[0].force[0])
            force_matrix[self.elements[i].dof[0][1]] = float(self.elements[i].knots[0].force[1])
            force_matrix[self.elements[i].dof[1][0]] = float(self.elements[i].knots[1].force[0])
            force_matrix[self.elements[i].dof[1][1]] = float(self.elements[i].knots[1].force[1])
        
        return self.apply_outline_condition(force_matrix)
    
    def apply_outline_condition(self, matrix):

        for i in range(len(self.restriction) - 1, -1, -1):
            if matrix.shape[1] != 1: matrix = np.delete(matrix, int(self.restriction[i]), 1)
            matrix = np.delete(matrix, int(self.restriction[i]), 0)

        return matrix
        
    def methods(self, itr, tol, method):
        
        u_vector = np.zeros(len(self.force_matrix))
        x_vector = np.zeros(len(self.force_matrix))
        loop = True
        err = float('inf')

        for iteration in range(itr):
            for x in range(len(x_vector)):

                for j in range(len(x_vector)):
                    if x != j: x_vector[x] += self.rigity_matrix[x][j]*u_vector[j]
                    
                x_vector[x] = (self.force_matrix[x] - x_vector[x])/self.rigity_matrix[x][x]

                if x_vector[x] != 0: err = abs((x_vector[x] - u_vector[x])/x_vector[x])*100
                if err < tol:
                    loop = False
                    break

                if method == "gauss seidel": u_vector[x] = x_vector[x]

            if not loop: break

            if method == "jacobi":
                for j in range(len(x_vector)): u_vector[j] = x_vector[j]

        return_matrix = np.zeros((self.knot_numbers*2, 1))
        counter = 0

        for i in range(self.knot_numbers*2):
            if i in self.restriction: return_matrix[i] = 0
            else: 
                return_matrix[i] = u_vector[counter]
                counter += 1

        return return_matrix, err

    def calculate_deformation(self):

        e = np.zeros(self.elements_number)
        sigma = np.zeros(self.elements_number)

        for i in range(self.elements_number):
            dof_1 = self.elements[i].dof[0]
            dof_2 = self.elements[i].dof[1]
            u1, u2 = self.u_vector[dof_1[0]], self.u_vector[dof_2[0]]
            v1, v2 = self.u_vector[dof_1[1]], self.u_vector[dof_2[1]]
            dot_product = np.array([-self.elements[i].cos, -self.elements[i].sin, self.elements[i].cos, self.elements[i].sin]) @ np.array([u1, v1, u2, v2])
            
            e[i] = 1/self.elements[i].lenght * dot_product[0]
            sigma[i] = self.elements[i].elasticity * e[i]

        return e, sigma

    def calculate_internal_force(self): 
        
        internal_force = np.zeros(self.elements_number)
        
        for i in range(self.elements_number): internal_force[i] = self.tension[i] * self.elements[i].area

        return internal_force

    def calculate_reaction(self):

        reaction = self.rigity_matrix_no_outline.dot(self.u_vector)
        return_reaction = []

        for i in self.restriction: return_reaction.append(reaction[int(i)][0])

        return return_reaction   