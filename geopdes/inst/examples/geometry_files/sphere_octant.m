function srf = sphere_octant(radius)
circ = nrbcirc(radius,[0 0 0],0, pi/2);
srf = nrbrevolve(circ,[0 0 0],[1.0 0 0], pi/2);
end