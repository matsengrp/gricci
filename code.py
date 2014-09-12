#
#    SAGE worksheet for computing the Ollivier-Ricci curvature by solving the integer linear problem
#
#    Applies to graphs (a list is provided), defines the random walk, and then the constraints
#    depending on a lazyness parameter t
#    returns the cost W1, the fixed t curvature and the asymptotic curvature
#
#    It can be adapted to deal with graphs with edges of non constant lengths,
#    or to different random walks
#
#    Pascal Romon 2014 — pascal.romon@u-pem.fr
#

# Here is a list of graphs
tetrahedron=Graph({2:[0,2,3], 2:[0,1,3], 3:[1,2,0],  0:[1,2,3]});
cube=graphs.HexahedralGraph();
cubex=graphs.HexahedralGraph(); # cube with diagonals from 0 and 1 added
cubex.add_edges([ (0,2),(0,5),(0,7),(1,3),(1,4),(1,6) ]);
octahedron=graphs.OctahedralGraph();
dodecahedron=graphs.DodecahedralGraph();
icosahedron=graphs.IcosahedralGraph();
# generic triangulations indexed by source degree and target degree
G34=Graph({0:[1,2,3],1:[0,2,3,4],2:[0,1,4],3:[0,1,4],4:[1,2,3]});
G44=Graph({0:[1,2,3,4],1:[0,2,3,5],2:[0,1,4,5],3:[0,1,4,5],4:[0,2,3],5:[1,2,3]});
G45=Graph({0:[1,2,3,4],1:[0,2,3,5,6],2:[0,1,4,6],3:[0,1,4,5],4:[0,2,3],5:[1,3,6],6:[1,2,5]});
G46=Graph({0:[1,2,3,4],1:[0,2,3,5,6,7],2:[0,1,4,7],3:[0,1,4,5],4:[0,2,3],5:[1,3,6],6:[1,5,7],7:[1,2,6]});
G55=Graph({0:[1,2,3,4,5],1:[0,2,3,6,7],2:[0,1,4,7],3:[0,1,5,6],4:[0,2,5],5:[0,3,4],6:[1,3,7],7:[1,2,6]});
G56=Graph({0:[1,2,3,4,5],1:[0,2,3,6,7,8],2:[0,1,4,8],3:[0,1,5,6],4:[0,2,5],5:[0,3,4],6:[1,3,7],7:[1,6,8],8:[1,2,7]});
# regular tilings
R4=Graph({ 0:[1,2,4,6],1:[7,9,11],2:[3,11],4:[3,5],6:[5,7],8:[7,9],10:[9,11]}); # square lattice
R6=Graph({ 0:[1,2,6], 1:[9,13], 2:[3,15], 4:[3,5],6:[5,7],8:[7,9],10:[9,11],12:[11,13],14:[13,15] }); # hexagonal lattice
# semiregular tilings
SnubSquare=Graph({ 0:[1,2,3,5,6],1:[2,3,9,10],2:[7,8,12],3:[4,11],5:[4,6],7:[6,12],8:[9,12],10:[9,11]}); # two types of edges: 0-1 between triangles and 1-2 between square and triangle
g=SnubSquare; N=g.order(); g.show()

# matrix of distances
D=g.distance_matrix(); view(D)

# t=t1/t2 is the lazyness; matrices are multiplied by t2, and by the degrees of source and target, in order to be integer-valued
# due to linear behaviour for small t, t=1/4 is sufficient here
t1=1; t2=4; t=t1/t2
# Set up linear program.
p=MixedIntegerLinearProgram()
# x[i,j] is the amount of mass that goes from i to j.
# It is constrained to be nonnegative, which gives the only inequalities for our LP.
x=p.new_variable(nonnegative=True)
# Maximize the negative of the mass transport.
p.set_objective(-p.sum(D[i,j]*x[i,j] for i in [0..N-1] for j in [0..N-1]))

# source and target vertices for Ollivier-Ricci, ds,dt their degrees we often
# choose 0 and 1, except for the snub square where we also pick 1 and 2
source=0; target=1; ds=g.degree(source); dt=g.degree(target);

# Mass distribution at j (multiplied by t2*ds*dt) of one time-step of a
# discrete lazy random walk starting at i.
# This is (a multiple of a) random walk $m$ as defined in Ollivier's paper,
# where it would be denoted $m_i(j)$.
def m(i,j):
    if i==j: return (t2-t1)*ds*dt
    elif D[i,j]==1: return t1*ds*dt/g.degree(i)
    else: return 0

# The equality constraints simply state that the mass starts in $m_source$ and
# finishes in $m_target$.
for i in [0..N-1]:
     p.add_constraint( p.sum( x[i,j] for j in [0..N-1] ) == m(source,i) )
for j in [0..N-1]:
     p.add_constraint( p.sum( x[i,j] for i in [0..N-1] ) == m(target,j) )

# RESULTS
# cost, kappa(t), ric
W1=-QQ(p.solve())/(t2*ds*dt); kappa=1-W1; ric=kappa/t; (W1,kappa,ric)

# optimal coupling
X=matrix(QQ,N,N) # preparing for the optimal coupling matrix X
for i in [0..N-1]: # filling X in
    for j in [0..N-1]:
        X[i,j] = p.get_values(x[i,j])/(t2*ds*dt)

view(X)

