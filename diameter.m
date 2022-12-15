function D = diameter(Simplex)
Z = squareform(pdist(Simplex));
D = max(max(Z));