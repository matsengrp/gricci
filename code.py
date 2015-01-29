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

from sage.all import QQ, ZZ, MixedIntegerLinearProgram, matrix, LCM
from collections import OrderedDict, namedtuple

Problem = namedtuple(
    'Problem',
    ['src', 'm_src',  # Source node and mass distribution.
     'tgt', 'm_tgt',  # Target node and mass distribution.
     'denom'  # A denominator: multiplying by this makes everything integral.
     ])

Result = namedtuple(
    'Result',
    ['coupling',  # How mass gets moved about.
     'dist',  # Distance between src and tgt.
     'W1',  # Wasserstein distance between random walks.
     'kappa',  # Coarse Ricci curvature as in Ollivier 2009.
     'ric'  # Curvature potentially normalized by a factor.
     ])


def ricci_gen(g, problem):
    # relevant_verts will have unique items starting with source, then target,
    # then the neighbors of source and target.
    # Use an OrderedDict to get uniques but preserve order.
    relevant_verts = OrderedDict.fromkeys(
        [problem.src, problem.tgt] +
        list(g.neighbor_iterator(problem.src)) +
        list(g.neighbor_iterator(problem.tgt))).keys()
    N = len(relevant_verts)
    D = matrix(ZZ, N, N)
    for i in range(N):
        for j in range(i, N):
            dist = g.distance(relevant_verts[i], relevant_verts[j])
            D[i, j] = dist
            D[j, i] = dist

    # We multiply through by problem.denom, so this checks for unit mass:
    assert(problem.denom == sum(problem.m_src[z] for z in relevant_verts))
    assert(problem.denom == sum(problem.m_tgt[z] for z in relevant_verts))

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
            problem.m_src[relevant_verts[i]])
    # and finishes in the places it has diffused to from target.
    for j in range(N):
        p.add_constraint(
            p.sum(a[i, j] for i in range(N)) ==
            problem.m_tgt[relevant_verts[j]])

    W1 = -QQ(p.solve())/problem.denom

    dist_src_tgt = D[0, 1]  # By def'n of relevant_verts.
    kappa = 1 - W1/dist_src_tgt
    return Result(
        coupling={
            (relevant_verts[i], relevant_verts[j]): QQ(v/problem.denom)
            for ((i, j), v) in p.get_values(a).items() if v > 0},
        dist=dist_src_tgt,
        W1=W1,
        kappa=kappa,
        ric=kappa)


def careful_div(i, j):
    """
    Integer division, making sure that j divides i.
    """
    assert(i >= j and i % j == 0)
    return int(i/j)


def lurw_mass(g, x, denom, pm1, pm2):
    """
    Under `lurw`, we use the lazy random walk with probability of movement
    pm=pm1/pm2 (expressed as rational to keep results rational) starting at x.
    This mass distribution is reported as it would be when we multiply
    everything by pm2/pm1. This is (a multiple of a) random walk $m$ as defined
    in Ollivier's paper, where it would be denoted $m_x(y)$.
    """
    m = [0] * g.order()
    mass = careful_div(denom*pm1, g.degree(x)*pm2)
    for y in g.neighbor_iterator(x):
        m[y] = mass
    m[x] = denom - sum(m)
    return m


def upmh_mass(g, x, denom):
    """
    One step from the MH process drawing from uniform distribution on vertices.
    Mass transport from a vertex x works as follows: From x, propose a
    uniform neighbor y. Uniform-prior MH ratio for proposing a move from x to
    y is min[1, (g(y -> x) / g(x -> y))] In our case, g(y -> x) = 1/d(y) and
    g(x -> y) = 1/d(x) giving MH ratio min[1, d(x) / d(y)]. Thus we always
    move to a lower degree node, and sometimes to a higher. The amount of
    mass moving from x to a neighboring y is
    min[1, d(x) / d(y)] / d(x)  (*)
    and the amount of mass that stays put is
    1 - sum_z min[1, d(x) / d(z)] / d(x)
    where z ranges over the neighbors of x.
    """
    m = [0] * g.order()
    dx = g.degree(x)
    for y in g.neighbor_iterator(x):
        m[y] = careful_div(
            min(denom, careful_div(dx*denom, g.degree(y))),
            dx)
    m[x] = denom - sum(m)
    return m


def ricci(walk, g, source=0, target=1, verbose=False):
    """
    Calculate the coarse Ollivier-Ricci curvature for vertices source and
    target of graph g. There are two options: lurw and upmh.
    For lurw, the results are reported in an ordered tuple including
    entries `kappa`, which is the coarse Ricci curvature as defined in Ollivier
    (2009), and `ric`, which is kappa divided by the probability of movement.
    This normalized version only depends on the graph. Under `upmh`,
    we use the the Metropolis-Hastings setup for the uniform prior on nodes.
    There are no free parameters, so in this case we just take `ric` to be
    `kappa`.
    """
    ds = g.degree(source)
    dt = g.degree(target)

    if walk == 'lurw':
        # pm=pm1/pm2 is the probability of movement.
        # Here we take 1/4, but we just have to be < 1.
        pm1 = 1
        pm2 = 4
        denom = pm2*ds*dt

        problem = Problem(
            src=source,
            m_src=lurw_mass(g, source, denom, pm1, pm2),
            tgt=target,
            m_tgt=lurw_mass(g, target, denom, pm1, pm2),
            denom=denom)

    elif walk == 'upmh':
        # Let lcm = the least common multiple of degrees of the relevant nodes.
        # Multiplying (*) above by ds*dt*lcm, we get the amount of mass moving
        # from x to y as
        # (ds * dt / d(x)) min[lcm, d(x) * lcm / d(y)],
        # both terms of which are integral (recall x is either src or tgt).

        lcm = LCM(
            [g.degree(z) for z in g.neighbor_iterator(source)] +
            [g.degree(z) for z in g.neighbor_iterator(target)])
        denom = ds*dt*lcm

        problem = Problem(
            src=source,
            m_src=upmh_mass(g, source, denom),
            tgt=target,
            m_tgt=upmh_mass(g, target, denom),
            denom=denom)

    else:
        assert(False)

    if verbose:
        print problem

    result = ricci_gen(g, problem)

    if walk == 'lurw':  # Normalize out the laziness.
        return result._replace(ric=QQ(result.kappa * pm2/pm1))

    return result


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
