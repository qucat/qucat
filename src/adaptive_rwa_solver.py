"""
author: Brian Baker
email:  bbaker.1796@gmail.com

The 'adaptive-RWA solver' module provides the routines needed to solve
for the density matrix for an open quantum system driven by an external 
oscillatory filed. The approach is to determine the
rotating frame deemed best suited to solve the steady state master equation in, and
then update that frame iteratively until one is converged upon. 
"""
from __future__ import division
from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
import qutip as qt  ## Import qutip packages

import multiprocessing as mp
import time

from numpy import linalg as nla
from scipy.sparse import linalg as la
from scipy.sparse.linalg import spsolve
from functools import partial


class DrivenOpenSystem:
    """Class for an open quantum system subject to a sinusoidal drive tone. Bare system Hamiltonian and drive Hamiltonian are both
       of type QObj.
       Expected parameters:
            hamiltonian: (QObj) system Hamiltonian
            decoherence_rates: (array or list) of system decoherence rates 
            drive_hamiltonian: (QObj) drive Hamiltonian that couples to external sinusoidal field
            angular_frequency: (float) frequency of modulation

    """

    def __init__(
        self, hamiltonian, decoherence_rates, drive_hamiltonian, angular_frequency
    ):

        if isinstance(hamiltonian, qt.Qobj) and isinstance(drive_hamiltonian, qt.Qobj):
            self.hamiltonian = hamiltonian
            self.v = drive_hamiltonian
        else:
            raise TypeError(
                "The system Hamiltonian and drive Hamiltonian must both be a quantum object"
            )

        self.dim = self.hamiltonian.shape[0]
        self.gamma_table = decoherence_rates
        self.omega = angular_frequency
        self.evals, self.evecs = self.hamiltonian.eigenstates()

        self.jump_ops = self.jump_operators(self.gamma_table)
        # This initial steady-state is the thermal equilibrium state (bootstrapping)
        self.rho_s = qt.steadystate(self.hamiltonian, self.jump_ops)
        self.density_matrix = [
            [
                (self.evecs[a].dag() * self.rho_s * self.evecs[b])[0]
                for b in range(self.dim)
            ]
            for a in range(self.dim)
        ]
        self.density_matrix = np.asarray(self.density_matrix)

        self.integer_list = [0 for k in range(self.dim)]
        self.drive_term_list = []
        self.rotating_frame_hamiltonian = self.hamiltonian

    def adaptive_rwa_solve(
        self, cutoff_matrix_element=0.0, it_max=10, neglected_term_info=False
    ):
        """run steady-state master equation iteratively until the same terms in drive Hamiltonian are kept
           returns: steady-state density matrix; information on neglected terms and a measure of the deviation from
                    the exact solution based on the relevance of these terms (optional return if user wants this information) 
           @param cutoff_matrix_element: (float) smallest matrix element size desired to be considered in determing drive term relevance
           @param it_max: (int) maximum iteration count to avoid convergence issues
           @param neglected_term_info: (bool) set to True if want to return both a list of neglected terms and the measure of error 

        """
        evecs = self.evecs
        evals = self.evals
        iteration_count = 0
        while iteration_count <= it_max:

            (
                self.rotating_frame_hamiltonian,
                self.integer_list,
                new_drive_term_list,
                neglected_term_list,
                error_parameter,
            ) = self.determine_frame(cutoff_matrix_element)

            if new_drive_term_list == self.drive_term_list:
                break
            self.drive_term_list = new_drive_term_list[:]

            # Solve the steady-state master equation in chosen rotating frame
            self.rho_s = qt.steadystate(self.rotating_frame_hamiltonian, self.jump_ops)
            # Update density matrix to determine new rotating frame
            self.density_matrix = [
                [(evecs[a].dag() * self.rho_s * evecs[b])[0] for b in range(self.dim)]
                for a in range(self.dim)
            ]
            self.density_matrix = np.asarray(self.density_matrix)
            iteration_count += 1

        if iteration_count > it_max:
            warnings.warn("maximum iteration count exceeded")

        if neglected_term_info:
            return self.rho_s, neglected_term_list, error_parameter
        else:
            return self.rho_s

    def jump_operators(self, gamma_table):
        """Sets up list of collapse operators for input into master equation
           returns: list of Lindblad jump operators
           @param gamma_table: (array or list of floats) table of decoherence rates between all eigenstates of system Hamiltonian
        """
        evecs = self.evecs
        evals = self.evals

        L = []

        for k in range(self.dim):
            for j in range(self.dim):
                L.append(np.sqrt(abs(gamma_table[j][k])) * evecs[k] * evecs[j].dag())

        return L

    def determine_frame(self, cutoff_matrix_element):
        """Construct a rotating frame transformation with corresponding integers, Hamiltonian, and drive terms
           returns: rotating frame Hamiltonian, h; list of chosen integers; list of chosen drive terms; list of neglected terms; error parameter
           @param cutoff_matrix_element: (float) smallest matrix element size desired to be considered in determing drive term relevance
        """

        evecs = self.evecs
        evals = self.evals
        v = self.v

        ranking, integers = self.rank_and_assign(cutoff_matrix_element)

        h_bare = 0
        for j in range(self.dim):
            h_bare += (evals[j]) * (evecs[j] * evecs[j].dag()) - self.omega * integers[
                j
            ] * (evecs[j] * evecs[j].dag())

        h_drive = 0
        error_parameter = 0
        number_of_transitions = int(self.dim * (self.dim - 1) / 2)
        drive_terms = []
        neglected_terms = []

        for j in range(number_of_transitions):
            index = ranking[j][1]
            if index[0] != index[1]:
                if integers[max(index)] == integers[min(index)] + 1:
                    h_drive += (
                        evecs[min(index)]
                        * evecs[max(index)].dag()
                        * (evecs[min(index)].dag() * (v * evecs[max(index)]))
                    )
                    drive_terms.append(
                        evecs[min(index)]
                        * evecs[max(index)].dag()
                        * (evecs[min(index)].dag() * (v * evecs[max(index)]))
                    )
                else:
                    neglected_terms.append(
                        evecs[min(index)]
                        * evecs[max(index)].dag()
                        * (evecs[min(index)].dag() * (v * evecs[max(index)]))
                    )
                    error_parameter += ranking[j][0] ** 2
            else:
                break
        error_parameter = np.sqrt(error_parameter)

        return (
            h_bare + h_drive + h_drive.dag(),
            integers,
            drive_terms,
            neglected_terms,
            error_parameter,
        )

    def rank_and_assign(self, cutoff_matrix_element):
        """rank individual drive terms based on their calculated relevence parameter and assign integers
           returns: transition_rank, a list containing (in order of decreasing relevance) all relevance parameters for each term and the indices of each term.
                    integers, a list containing the integer assigned for each state
           @param cutoff_matrix_element: (float) smallest matrix element size desired to be considered in determing drive term relevance
        """

        L0 = qt.liouvillian(self.rotating_frame_hamiltonian, self.jump_ops)

        relevance_table = [
            [
                self.calculate_first_order_correction(cutoff_matrix_element, n, m, L0)
                for m in range(self.dim)
            ]
            for n in range(self.dim)
        ]
        relevance_table = np.asarray(relevance_table)

        number_of_transitions = int(self.dim * (self.dim - 1) / 2)
        transition_rank = [None for i in range(number_of_transitions)]
        # This loop ranks drive terms according to relevance
        for rank in range(number_of_transitions):
            max_ranked_indices = np.where(relevance_table == relevance_table.max())
            indices = [max_ranked_indices[0][0], max_ranked_indices[1][0]]
            transition_rank[rank] = [relevance_table.max(), indices]
            relevance_table[indices[0]][indices[1]] = relevance_table[indices[1]][
                indices[0]
            ] = 0

        # This graphical algorithm assigns an integer to each eigenstate of the Hamiltonian based on the ranking from above
        integer_list = [None for i in range(self.dim)]
        # START ALGORITHM
        # initialize first term into a graph
        first_index = transition_rank[0][1]
        graph_list = [[first_index[0], first_index[1]]]
        integer_list[max(first_index)] = 1
        integer_list[min(first_index)] = 0
        # assign subsequent terms
        for i in range(1, number_of_transitions):
            # if no more non-zero relevance parameters, then break
            if transition_rank[i][0] == 0.0:
                break
            else:
                index = transition_rank[i][1]
                # scenario (i) neither states have been incorporated into the graph
                if integer_list[index[0]] == integer_list[index[1]] == None:
                    integer_list[max(index)] = 1
                    integer_list[min(index)] = 0
                    # place them in a new graph
                    graph_list.append([index[0], index[1]])
                # scenario (ii) one of the states has been incorporated, but not the other
                elif integer_list[index[0]] == None:
                    if index[0] > index[1]:
                        integer_list[index[0]] = integer_list[index[1]] + 1
                    else:
                        integer_list[index[0]] = integer_list[index[1]] - 1
                    # find which graph component to put the state in (the component the other state is in)
                    for k, graph in enumerate(graph_list):
                        if index[1] in graph:
                            # place that state in that graph component
                            graph_list[k].append(index[0])
                            break
                elif integer_list[index[1]] == None:
                    if index[0] > index[1]:
                        integer_list[index[1]] = integer_list[index[0]] - 1
                    else:
                        integer_list[index[1]] = integer_list[index[0]] + 1
                    for k, graph in enumerate(graph_list):
                        if index[0] in graph:
                            graph_list[k].append(index[1])
                            break
                # scenario (iii) both states have already been incorporated in the graph
                else:
                    # find the graph components where these states have been placed
                    for k, graph in enumerate(graph_list):
                        overlap = list(set(index) & set(graph))
                        # subscenario: the states are in the same graph component, hence a cycle so nothing can do
                        if len(overlap) == 2:
                            break
                        # subscenario: the states are in two disjoint graph components
                        elif len(overlap) == 1:
                            fixed_index = overlap[0]
                            shift_index = list(set(index) - set(graph))[0]
                            old_integer = integer_list[shift_index]
                            if shift_index > fixed_index:
                                new_integer = integer_list[fixed_index] + 1
                            else:
                                new_integer = integer_list[fixed_index] - 1
                            shift_amount = new_integer - old_integer
                            # merge one graph component into the other
                            for j, graph2 in enumerate(graph_list):
                                if shift_index in graph2:
                                    for m, index2 in enumerate(graph2):
                                        integer_list[index2] = (
                                            integer_list[index2] + shift_amount
                                        )
                                    graph_list[k] = graph_list[k] + graph2
                                    graph_list.pop(j)
                                    break
                            break
                        else:
                            continue
                    continue
        # Just in case, if a state was not assigned an integer due to not participating in dynamics, set its integer to 0
        for i, integer in enumerate(integer_list):
            if integer == None:
                integer_list[i] = 0
        ## END algorithm
        return transition_rank, integer_list

    def calculate_first_order_correction(self, cutoff_matrix_element, n, m, L0):
        """Calculates the first order correction to the steady-state density matrix due to drive term,
           and thereby calculating the corresponding relevance parameter
           returns: relevance parameter, Delta_{nm}
           @param cutoff_matrix_element: (float)
           @param n,m: (int) the indices of the drive term
           @param L0: (QObj) the Liouvillian superoperator L_0 needed to solve for the correction, independent of indices
                             so not needed to be re-evaluated every time this function is called
        """
        if n >= m:
            return 0.0

        evecs = self.evecs
        evals = self.evals

        # ignore drive terms whose matrix elements are beneath a specificied cutoff for speed-up.
        v_nm = (evecs[n].dag() * (self.v * evecs[m]))[0][0][0]
        if abs(v_nm) <= cutoff_matrix_element:
            return 0.0

        k = self.integer_list
        rho_s_vectorform = np.reshape(
            self.density_matrix, (self.dim ** 2, 1), order="F"
        )

        V_nm = evecs[n] * evecs[m].dag() * (evecs[n].dag() * (self.v * evecs[m]))
        L_nm = qt.liouvillian(V_nm)
        # b = np.dot(L_nm.full(),rho_0)
        b = (L_nm * rho_s_vectorform).data
        omega_of_k = (k[n] - k[m] + 1) * self.omega

        A = 1j * omega_of_k * qt.identity(self.dim ** 2).data - L0.data

        # A = A.full()
        # del_rho = la.lstsq(A,b,rcond = 1e-6)[0]

        if omega_of_k == 0:
            del_rho = la.lsmr(A, b)[0]
        else:
            del_rho = spsolve(A, b)

        return nla.norm(del_rho)
