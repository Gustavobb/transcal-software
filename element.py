from math import sqrt
import numpy

class element:

    def __init__(self, elasticity, area, knots, dof):

        self.knots = knots
        self.dof = dof
        self.area = area
        self.elasticity = elasticity
        self.lenght = sqrt((self.knots[0].x - self.knots[1].x)**2 + (self.knots[0].y - self.knots[1].y)**2)
        self.cos = (knots[1].x - knots[0].x)/self.lenght
        self.sin = (knots[1].y - knots[0].y)/self.lenght
        self.rigity_matrix = self.calculate_rigity_matrix()

    def calculate_rigity_matrix(self):
        
        matrix = numpy.array(
            [[self.cos**2, self.cos * self.sin, -self.cos**2, -self.cos * self.sin],
            [self.cos * self.sin, self.sin**2, -self.cos * self.sin, -self.sin**2],
            [-self.cos**2, -self.cos * self.sin, self.cos**2, self.cos * self.sin],
            [-self.cos * self.sin, -self.sin**2, self.cos * self.sin, self.sin**2]])

        return (self.elasticity * self.area)/self.lenght * matrix