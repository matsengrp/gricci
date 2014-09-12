# By Pascal Romon 2014 (pascal.romon@u-pem.fr)

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


