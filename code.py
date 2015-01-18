"""
SAGE functions for computing the Ollivier-Ricci curvature on graphs by solving
the corresponding integer linear problem.

Reference:
Ollivier Y (2009) Ricci curvature of Markov chains on metric spaces.
J Funct Anal 256: 810–864.

Original code by Pascal Romon 2014 (pascal.romon@u-pem.fr)
Tidied, wrapped, generalized and made faster by Erick Matsen 2014-2015
(http://matsen.fredhutch.org/)
"""

from sage.all import QQ, ZZ, MixedIntegerLinearProgram, matrix
from collections import OrderedDict, namedtuple


def ricci(walk, g, source=0, target=1, verbose=False):
    """
    Calculate the coarse Ollivier-Ricci curvature for vertices source and
    target of graph g. There are two options: lazy_unif and unif_prior_mh.

    Under `lazy_unif`, we use the lazy random walk with probability of movement
    pm=pm1/pm2 (expressed as rational to keep results rational).
    The results are reported in an ordered tuple including entries `kappa`,
    which is the coarse Ricci curvature as defined in Ollivier (2009), and
    `ric`, which is kappa divided by pm. This normalized version only depends
    on the graph.

    Under `unif_prior_mh`, we use the the Metropolis-Hastings setup for the
    uniform prior on nodes. There are no free parameters, so in this case we
    just take `ric` to be `kappa`.
    """

    # Yes, this will become an Enum when they are available.
    assert(walk == 'lazy_unif' or walk == 'unif_prior_mh')

    # Use an OrderedDict to get uniques but preserve order.
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

    # -- lazy_unif --

    # Mass distribution at y (multiplied by pm2*ds*dt) of one time-step of a
    # discrete lazy random walk starting at x.
    # This is (a multiple of a) random walk $m$ as defined in Ollivier's paper,
    # where it would be denoted $m_x(y)$.
    def m_lazy_unif(x, y):
        if x == y:
            return (pm2-pm1)*ds*dt
        elif g.has_edge(x, y):
            # Note that g.degree(x) is either ds or dt.
            return pm1*ds*dt/g.degree(x)
        else:
            return 0

    if walk == 'lazy_unif':
        # pm=pm1/pm2 is the probability of movement.
        pm1 = 1
        pm2 = 4
        # Matrices are multiplied by mass_denominator in order to be
        # integer-valued.
        mass_denominator = pm2*ds*dt
        m = m_lazy_unif

    # -- unif_prior_mh --

    # Mass transport from a vertex x works as follows: From x, propose a
    # uniform neighbor y. Uniform-prior MH ratio for proposing a move from x to
    # y is min[1, (g(y -> x) / g(x -> y))] In our case, g(y -> x) = 1/d(y) and
    # g(x -> y) = 1/d(x) giving MH ratio min[1, d(x) / d(y)]. Thus we always
    # move to a higher degree node, and sometimes to a lower. The amount of
    # mass moving from x to a neighboring y is min[1, d(x) / d(y)] / d(x). We
    # can then keep things integral by multiplying by d(x) d(y), such that mass
    # movement is min[d(x), d(y)], and the amount of mass that stays put is
    # d(x) d(y) - sum_z min[d(x), d(z)] where z ranges over the neighbors of x.
    # For example, no mass stays put if every neighbor of x has higher degree.

# NOTE: below does not yet reflect the above.

    def m_unif_prior_mh(x, y):
        """
        This function assumes that x is either source or target.
        """
        assert(x == source or x == target)
        # rn is for "remaining numerator" when we multiply the ratio by
        # mass_denominator. The rest gets absorbed into the calculation as
        # described just above.
        rn = ds*ds if x == target else dt*dt
        dx = g.degree(x)
        if x == y:
            tot = sum(min(dx, g.degree(z)) for z in g.neighbor_iterator(x))
            return rn*(dx*dx - tot)
        elif g.has_edge(x, y):
            return rn*min(dx, g.degree(y))
        else:
            return 0

    if walk == 'unif_prior_mh':
        # In the above, we don't know what x is a priori (it could be either
        # source or target), so we multiply the amount of mass by the product
        # of both squares:
        mass_denominator = ds*ds*dt*dt
        m = m_unif_prior_mh

    # Just handy for taking a look.
    def mass_vector(x):
        v = [(y, QQ(m(x, y))/mass_denominator) for y in relevant_verts]
        return [(xm, mass) for (xm, mass) in v if mass > 0]
    if verbose:
        print "source mass: ", mass_vector(source)
        print "target mass: ", mass_vector(target)

    # Set up linear program.
    p = MixedIntegerLinearProgram()
    # Note that here and in what follows, i and j are used as the vertices that
    # correspond to relevant_verts[i] and relevant_verts[j].
    # a[i,j] is the (unknown) amount of mass that goes from i to j.
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

    def relevant_vert_pair_labels(src, dst):
        return (relevant_verts[src], relevant_verts[dst])

    p.solve()
    W1 = -QQ(p.solve())/mass_denominator
    # Below D[0, 1] is Dist(source, target) by def of relevant_verts.
    kappa = 1 - W1/D[0, 1]
    if walk == 'lazy_unif':  # Normalize out the laziness.
        ric = QQ(kappa * pm2/pm1)
    else:  # Nothing to normalize out.
        ric = kappa
    Result = namedtuple('Result',
                        ['verts', 'coupling', 'dist', 'W1', 'kappa', 'ric'])
    return Result(
        verts=relevant_verts,
        coupling={
            relevant_vert_pair_labels(*k): QQ(v/mass_denominator)
            for (k, v) in p.get_values(a).items() if v > 0},
        dist=D[0, 1],
        W1=W1,
        kappa=kappa,
        ric=ric)


def ricci_list(walk, g, pair_list):
    return[(s, t, ricci(walk, g, source=s, target=t).ric)
           for (s, t) in pair_list]


def ricci_matrix(walk, g):
    N = g.num_verts()
    m = matrix(QQ, N)
    for i in range(N):
        for j in range(i+1, N):
            r = ricci(walk, g, source=i, target=j).ric
            m[i, j] = r
            m[j, i] = r
    return m
