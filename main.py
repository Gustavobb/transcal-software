from funcoesTermosol import importa, plota, geraSaida
from knot import knot
from element import element
from structure import structure
import numpy as np

[knot_numbers, knot_coordinates, element_numbers, incidence, load_numbers, forces, restriction_numbers, restriction] = importa('entry.xlsx')

knot_list = []
element_list = []

for i in range(knot_numbers):

    # melhor codigo do mundo
    knot_list.append(knot(knot_coordinates[0][i], knot_coordinates[1][i], [forces[2 * i], forces[2*i + 1]]))

for i in range(element_numbers):

    incidence0 = (int(incidence[i][0]) - 1)*2
    incidence1 = (int(incidence[i][1]) - 1)*2
    element_list.append(element(incidence[i][2], incidence[i][3], [knot_list[int(incidence[i][0]) - 1], knot_list[int(incidence[i][1]) - 1]], [[incidence0, incidence0 + 1], [incidence1, incidence1 + 1]]))

object_structure = structure(element_list, knot_numbers, restriction, element_numbers)

geraSaida(object_structure.reaction, object_structure.u_vector, object_structure.deformation, object_structure.internal_force, object_structure.tension)

knot_coordinates_after = np.zeros((2, len(knot_coordinates[0])))

for i in range(int(len(object_structure.u_vector)/2)):
    knot_coordinates_after[0][i] = knot_coordinates[0][i] - object_structure.u_vector[2*i]
    knot_coordinates_after[1][i] = knot_coordinates[1][i] - object_structure.u_vector[2*i + 1]

plota(knot_coordinates_after, incidence)