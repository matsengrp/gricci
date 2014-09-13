"""
SAGE functions for computing the Ollivier-Ricci curvature by solving the
integer linear problem.

Original code by Pascal Romon 2014 (pascal.romon@u-pem.fr)
Tidied and wrapped by Erick Matsen 2014 (http://matsen.fhcrc.org/)
"""

from sage.all import QQ, MixedIntegerLinearProgram, matrix
from collections import namedtuple


def ric_unif_rw(g, source=0, target=1, t1=1, t2=4):
    """
    Calculate the Ollivier-Ricci curvature for vertices source and target of
    graph g under the standard random walk (selecting moves uniformly).
    t1 and t2 are the numerator and denominator of t; due to linear behaviour
    for small t, t=1/4 is sufficient here.
    """
    N = g.order()
    D = g.distance_matrix()
    ds = g.degree(source)
    dt = g.degree(target)
    # t=t1/t2 is the laziness.
    t = float(t1)/t2
    # Matrices are multiplied by mass_denominator in order to be
    # integer-valued.
    mass_denominator = t2*ds*dt
    # Set up linear program.
    p = MixedIntegerLinearProgram()
    # x[i,j] is the amount of mass that goes from i to j.
    # It is constrained to be nonnegative, which are the only inequalities for
    # our LP.
    x = p.new_variable(nonnegative=True)
    # Maximize the negative of the mass transport.
    p.set_objective(
        -p.sum(D[i, j]*x[i, j] for i in range(N) for j in range(N)))

    # Mass distribution at j (multiplied by t2*ds*dt) of one time-step of a
    # discrete lazy random walk starting at i.
    # This is (a multiple of a) random walk $m$ as defined in Ollivier's paper,
    # where it would be denoted $m_i(j)$.
    def m(i, j):
        if i == j:
            return (t2-t1)*ds*dt
        elif D[i, j] == 1:
            return t1*ds*dt/g.degree(i)
        else:
            return 0

    # Just handy for taking a look.
    def mass_vector(i):
        return [QQ(m(i, j))/mass_denominator for j in range(N)]

    # The equality constraints simply state that the mass starts in
    # $m_source...
    for i in range(N):
        p.add_constraint(p.sum(x[i, j] for j in range(N)) == m(source, i))
    # and finishes in $m_target$.
    for j in range(N):
        p.add_constraint(p.sum(x[i, j] for i in range(N)) == m(target, j))

    p.solve()
    W1 = -QQ(p.solve())/mass_denominator
    kappa = 1 - W1/D[source, target]
    Result = namedtuple('Result', ['coupling', 'W1', 'kappa', 'ric'])
    return Result(
        coupling={
            k: QQ(v/mass_denominator)
            for (k, v) in p.get_values(x).items() if v > 0},
        W1=W1,
        kappa=kappa,
        ric=QQ(kappa/t))


def ricci_list(g, pair_list):
    return[(s, t, ric_unif_rw(g, source=s, target=t).ric)
           for (s, t) in pair_list]


def ricci_matrix(g):
    N = g.num_verts()
    m = matrix(QQ, N)
    for i in range(N):
        for j in range(i+1, N):
            r = ric_unif_rw(g, source=i, target=j).ric
            m[i, j] = r
            m[j, i] = r
    return m
