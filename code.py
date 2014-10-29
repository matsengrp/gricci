"""
SAGE functions for computing the Ollivier-Ricci curvature by solving the
integer linear problem.

Original code by Pascal Romon 2014 (pascal.romon@u-pem.fr)
Tidied, wrapped and made faster by Erick Matsen 2014 (http://matsen.fhcrc.org/)
"""

from sage.all import QQ, ZZ, MixedIntegerLinearProgram, matrix
from collections import OrderedDict, namedtuple


def ric_unif_rw(g, source=0, target=1, t1=1, t2=4):
    """
    Calculate the Ollivier-Ricci curvature for vertices source and target of
    graph g under the standard random walk (selecting moves uniformly).
    t1 and t2 are the numerator and denominator of t; due to linear behaviour
    for small t, t=1/4 is sufficient here.
    """

    d = OrderedDict.fromkeys(
        [source, target] +
        list(g.neighbor_iterator(source)) +
        list(g.neighbor_iterator(target)))
    # relevant_verts will have unique items starting with source, then target,
    # then the neighbors of source and target.
    relevant_verts = d.keys()
    N = len(relevant_verts)
    D = matrix(ZZ, N, N)
    for i in range(N):
        for j in range(i, N):
            dist = g.distance(relevant_verts[i], relevant_verts[j])
            D[i, j] = dist
            D[j, i] = dist
    ds = g.degree(source)
    dt = g.degree(target)
    # t=t1/t2 is the laziness.
    t = float(t1)/t2
    # Matrices are multiplied by mass_denominator in order to be
    # integer-valued.
    mass_denominator = t2*ds*dt
    # Set up linear program.
    p = MixedIntegerLinearProgram()
    # Note that here and in what follows, i and j are used as the vertices that
    # correspond to relevant_verts[i] and relevant_verts[j].
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
            return t1*ds*dt/g.degree(relevant_verts[i])
        else:
            return 0

    # Just handy for taking a look.
    def mass_vector(i):
        return [QQ(m(i, j))/mass_denominator for j in range(N)]

    # The equality constraints simply state that the mass starts in
    # m_source, which is relevant_verts[0]...
    for i in range(N):
        p.add_constraint(p.sum(x[i, j] for j in range(N)) == m(0, i))
    # and finishes in m_target, which is relevant_verts[1].
    for j in range(N):
        p.add_constraint(p.sum(x[i, j] for i in range(N)) == m(1, j))

    p.solve()
    W1 = -QQ(p.solve())/mass_denominator
    # Below D[0, 1] is Dist(source, target) by def of relevant_verts.
    kappa = 1 - W1/D[0, 1]
    Result = namedtuple('Result',
                        ['verts', 'coupling', 'dist', 'W1', 'kappa', 'ric'])
    return Result(
        verts=relevant_verts,
        coupling={
            k: QQ(v/mass_denominator)
            for (k, v) in p.get_values(x).items() if v > 0},
        dist=D[0, 1],
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
