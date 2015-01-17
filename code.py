"""
SAGE functions for computing the Ollivier-Ricci curvature on graphs by solving
the corresponding integer linear problem.

Reference:
Ollivier Y (2009) Ricci curvature of Markov chains on metric spaces.
J Funct Anal 256: 810–864.

Original code by Pascal Romon 2014 (pascal.romon@u-pem.fr)
Tidied, wrapped, and made faster by Erick Matsen 2014-2015
(http://matsen.fredhutch.org/)
"""

from sage.all import QQ, ZZ, MixedIntegerLinearProgram, matrix
from collections import OrderedDict, namedtuple


def ric_unif_prior_rw(g, source=0, target=1):
    """
    Calculate the coarse Ollivier-Ricci curvature for vertices source and
    target of graph g under the the Metropolis-Hastings setup for the uniform
    prior on nodes.
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

    # Mass transport from a vertex x works as follows: From x, propose a
    # uniform neighbor y. Prior MH ratio for proposing a move from x to y is
    # min[1, (g(y -> x) / g(x -> y))] In our case, g(y -> x) = 1/d(y) and g(x
    # -> y) = 1/d(x) giving MH ratio min[1, d(y) / d(x)] The amount of mass
    # moving from x to a neighboring y is min[1, d(y) / d(x)] / d(x). We can
    # then keep things integral by multiplyng by d(x)^2, such that mass
    # movement is min[d(x), d(y)], and the amount of mass that stays put is
    # d(x)^2 - sum_i min[d(x), d(yi)]. For example, no mass stays put if every
    # neighbor of x has lower degree.

    # In the above, we don't know what x is a priori (it could be either source
    # or target), so we multiply the amount of mass by the product of both
    # squares:
    mass_denominator = ds*ds*dt*dt

    def m(x, y):
        """
        This function assumes that x is either source or target.
        """
        assert(x == source or x == target)
        # rn is for "remaining numerator" when we multiply the ratio by
        # mass_denominator. The rest gets absorbed into the calculation as
        # described just above.
        rn = ds*ds if x == target else dt*dt
        dx = g.degree(x)
        dy = g.degree(y)
        if x == y:
            tot = sum(min(dx, g.degree(z)) for z in g.neighbor_iterator(x))
            return rn*(dx*dx - tot)
        elif g.has_edge(x, y):
            return rn*min(dx, dy)
        else:
            return 0

    # Set up linear program.
    p = MixedIntegerLinearProgram()
    # Note that here and in what follows, i and j are used as the vertices that
    # correspond to relevant_verts[i] and relevant_verts[j].
    # a[i,j] is the amount of mass that goes from i to j.
    # It is constrained to be nonnegative, which are the only inequalities for
    # our LP.
    a = p.new_variable(nonnegative=True)
    # Maximize the negative of the mass transport.
    p.set_objective(
        -p.sum(D[i, j]*a[i, j] for i in range(N) for j in range(N)))
    # The equality constraints simply state that the mass starts in the places
    # it has diffused to from source...
    for i in range(N):
        p.add_constraint(
            p.sum(a[i, j] for j in range(N)) ==
            m(source, relevant_verts[i]))
    # and finishes in the places it has diffused to from target.
    for j in range(N):
        p.add_constraint(
            p.sum(a[i, j] for i in range(N)) ==
            m(target, relevant_verts[j]))

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
            for (k, v) in p.get_values(a).items() if v > 0},
        dist=D[0, 1],
        W1=W1,
        kappa=kappa,
        ric=QQ(kappa/pm))


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
