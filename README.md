# Conormal Complexes and Rings

Ardila, Denham and Huh develop the conormal fan and ring of a matroid. They use this to prove the log concavity of the h-vector of the broken circuit complex of a matroid. Work of Partida and Nathanson study the simplicial complex associated to the conormal fan, the conormal complex. This code computes the conormal complex and ring of a matroid. It also provides tools to work the building blocks of these objects, the biflats and biflags of a matroid.

For example, we can compute the conormal complex and check if it is flag:
```
  sage: G = Graph([['a','b'], ('b','a'), ['a','c'], ['b','c'],
  ....: ['b','d'], ['c','d']], multiedges=True)
  sage: M = Matroid(G)
  sage: C = conormal_complex(M)
  sage: C
  Simplicial complex with 30 vertices and 257 facets
  sage: C.is_flag_complex()
  False
```
We can also compute the conormal ring and find a Gröbner basis for it:
```
sage: G = Graph([['a','b'], ('b','a'), ['a','c'], ['b','c'], [
....: 'b','d'], ['c','d']], multiedges=True)
sage: M = Matroid(G)
sage: R = conormal_ring(M)
sage: R.defining_ideal().groebner_basis()
Polynomial Sequence with 256 Polynomials in 32 Variables
```
