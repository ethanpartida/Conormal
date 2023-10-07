def biflats(M):
    """
    Return a list containing all biflats of a loopless
    and coloopless matroid `M`.
    
    A biflat of a matroid consists of a pair `F|G` of a
    proper flat `F` of `M` and a proper coflat `G` of Mdual
    such that `F` union `G` is the groundset of `M` and 
    `F` intersect `G` is not the groundset of `M`. 
    
    INPUT:
    
    - ``M`` -- A loopless and coloopless matroid.
    
    OUTPUT:
    
    - ``biflats`` -- A list containing all biflats of M.
    
    EXAMPLES:
    sage: M = matroids.Wheel(3)
    sage: biflats(M)
    [(frozenset({0}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({1}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({2}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({3}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({4}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({5}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({0, 1, 3}), frozenset({2, 4, 5})),
    (frozenset({0, 2, 5}), frozenset({1, 3, 4})),
    (frozenset({1, 2, 4}), frozenset({0, 3, 5})),
    (frozenset({3, 4, 5}), frozenset({0, 1, 2})),
    (frozenset({0, 1, 3}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({0, 2, 5}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({0, 4}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({1, 2, 4}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({1, 5}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({2, 3}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({3, 4, 5}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({3})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({4})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({5})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({0})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({1})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({2})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({1, 3, 4})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({0, 3, 5})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({2, 3})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({2, 4, 5})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({0, 4})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({1, 5})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({0, 1, 2}))]

    """
    if M.loops() or M.coloops():
        raise Exception('M must be loopless and coloopless') 
    Mdual = M.dual()

    #Create a list of all pairs of flats and coflats
    possible_biflats = [(X,Y) for i in range(1,M.rank()+1) for j in
                        range(1,Mdual.rank()+1)
                        for X in M.flats(i) for Y in Mdual.flats(j)]

    #Filter the list to only include valid biflats
    biflats = [X for X in possible_biflats 
              if X[0].union(X[1]) == M.groundset()
              and not (X[0] == M.groundset() and X[1] == M.groundset())]
    return biflats

def biflat_comparison(X, Y):
    """
    Return if biflat `X` is less than or equal to biflat `Y`.
    
    Biflats form a partial order under the relations X_0|X_1 <= Y_0|Y_1
    if X_0 is a subset of Y_0 and Y_1 is a subset of X_1.

    INPUT:

    - ``X`` -- A biflat
    - ``Y`` -- A biflat
    
    OUTPUT:

    - ``B`` -- A boolean determining if X<= Y

    EXAMPLES:

    sage: M = matroids.Wheel(3)
    sage: BF = biflats(M)
    sage: BF[0]
    (frozenset({0}), frozenset({0, 1, 2, 3, 4, 5}))
    sage: BF[6]
    (frozenset({0, 1, 3}), frozenset({2, 4, 5}))
    sage: biflat_comparison(BF[0],BF[6])
    True

    
    """
    return X[0].issubset(Y[0]) and Y[1].issubset(X[1])

def biflat_poset(M):
    r"""
    Return the poset of biflats of `M`.

    The biflats of `M` form a partial order under the relations 
    X_0|X_1 <= Y_0|Y_1 if X_0 is a subset of Y_0 and Y_1 is a 
    subset of X_1. This function returns the biflats of `M` ordered
    under this partial order.
   
    INPUT:
    - ``M`` -- A loopless and coloopless matroid

    OUTPUT:
    - ``P`` -- A poset of biflats

    EXAMPLE:

    sage: M = matroids.Wheel(3)
    sage: biflat_poset(M)
    Finite poset containing 30 elements

    """
    
    biflat_list = biflats(M)
    P = Poset((biflat_list, lambda x,y : biflat_comparison(x,y)))
    return P

def check_biflag(c, M):
    r"""
    Return if a tuple `c` of biflats is a biflag of `M`.

    A tuple of biflats `c` of biflats is a biflag if `c` forms a chain
    in the poset of biflats of `M` and `c` satisfies a "gap" condition.
    See Definition 2.14 of Ardila, Denham and Huh's 
    "Lagrangian Geometry of Matroids".
    
    .. WARNING::

    This code checks if `c` is a bichain using the passed in order.
    If `c` forms a bichain undering a non-trivial permutation of its
    elements, this code will still return `False`.

    INPUT:
    
    - ``c`` -- A tuple of biflats of `M`
    - ``M`` -- A loopless and coloopless matroid
    
    OUTPUT:
    
    - ``B`` -- A boolean determining if `c` is a biflag of `M`

    EXAMPLES:

    sage: M = matroids.Wheel(3)
    sage: P = biflat_poset(M)
    sage: P.chains()[6]
    [(frozenset({0}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({0, 1, 3}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({0, 1, 3}), frozenset({2, 4, 5})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({2, 4, 5})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({2}))]
    sage: check_biflag(P.chains()[6],M)
    False
    sage: P.chains()[8]
    [(frozenset({0}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({0, 1, 3}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({0, 1, 3}), frozenset({2, 4, 5})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({5}))]
    sage: check_biflag(P.chains()[8],M)
    True
    """
    
    #Check if c forms a chain in the poset of biflats.
    for i in range(0,len(c)-1):
        if not biflat_comparison(c[i],c[i+1]):
            return False
        
    #Compute the union of the intersections of the biflats of c
    UofI = set().union(*[x[0].intersection(x[1]) for x in c])

    #Check the gap condition
    return  UofI != M.groundset()

#Returns the conormal complex of M
def conormal_complex(M):
    r"""
    Return the conormal complex of a matroid `M`.

    The conormal complex is the simplicial complex whose faces are
    the biflags of `M`. This is the simplicial complex associated to
    the simplical fan, the conormal fan.

    INPUT:

    - ``M`` -- A loopless and coloopless matroid

    OUTPUT:

    - ``C`` -- A simplicial complex whose faces are biflags of `M`.

    EXAMPLES:
    
    sage: M = matroids.Wheel(3)
    sage: conormal_complex(M)
    Simplicial complex with 30 vertices and 360 facets

    """
    
    P = biflat_poset(M) 

    #Find the inclusion-maximal biflags of `M`. By Lemma 2.1 of
    #"Lagrangian Combinatorics of Matroids", maximal biflags
    #all have length equal to the groundset of `M`-2.
    max_biflags = [tuple(c) for c in
            P.chains().elements_of_depth_iterator(len(M.groundset())-2)
                   if check_biflag(c,M)]

    #Return the simplicial complex with facets equal
    #to the maximal biflags.
    return SimplicialComplex(max_biflags)


def pretty_biflat(X, M):
    r"""
    Returns a readable version of a biflat `X` of a matroid `M`.

    Biflats are stored as tuples of frozensets. While nice for 
    manipulation, this format is not nicely readable. This
    function returns a readable string representing 'X'.

    INPUT:
    
    - `X` -- A biflat of `M`
    - `M` -- A loopless and coloopless matroid

    OUTPUT:

    - `s` -- A string representing the biflat `X`

    EXAMPLES:

    sage: M = matroids.Wheel(3)
    sage: x = biflats(M)[16]
    sage: x
    (frozenset({3, 4, 5}), frozenset({0, 1, 2, 3, 4, 5}))
    sage: pretty_biflat(x,M)
    '345uE'


    """
    if X[0] == M.groundset():
        top = 'E'
    else:
        top = ''.join([str(x) for x in sorted(X[0])])
    if X[1] == M.groundset():
        bot = 'E'
    else:
        bot = ''.join([str(x) for x in sorted(X[1])])
    return top + 'u' + bot


def pretty_biflag(c, M):
    r"""
    Returns a readable version of a biflag `c` of a matroid `M`.

    Biflags are stored as tuples of biflats, which are in turn
    strored as tuples of frozensets. While nice for 
    manipulation, this format is not nicely readable. This
    function returns a readable string representing 'c'.

    INPUT:
    
    - `c` -- A tuple of biflats of `M`
    - `M` -- A loopless and coloopless matroid

    OUTPUT:

    - `s` -- A string representing the biflag `c`

    EXAMPLES:
    
    sage: M = matroids.Wheel(3)
    sage: c = biflat_poset(M).chains()[6]
    sage: c
    [(frozenset({0}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({0, 1, 3}), frozenset({0, 1, 2, 3, 4, 5})),
    (frozenset({0, 1, 3}), frozenset({2, 4, 5})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({2, 4, 5})),
    (frozenset({0, 1, 2, 3, 4, 5}), frozenset({2}))]
    sage: pretty_biflag(c,M)
    ['0uE', '013uE', '013u245', 'Eu245', 'Eu2']

    
    """
    return [pretty_biflat(x, M) for x in c]

def conormal_ring(M, R=None):
    r"""
    Returns the conormal ring of a matroid `M` with `R` coefficients. 

    The conormal ring of `M` is the Chow ring of the tropical fan,
    the conormal fan of `M`. It is a quotient of the polynomial ring
    whose variables are the biflats of `M` along with 'E|empty' and
    'empty|E'. The relations are squarefree monomial relations coming
    from the Stanley--Reisner relations of the conormal complex along
    with linear relations.

    INPUT:
    
    - `M` -- A loopless and coloopless matroid
    - `R` -- A ring of coefficients (default is ZZ).

    OUTPUT:

    - `S` -- The conormal ring of M.

    EXAMPLES:
    
    sage: M = matroids.Uniform(2,3)
    sage: conormal_ring(M)
    Quotient of Multivariate Polynomial Ring in AEuempty, A0uE, A1uE, A2uE, AemptyuE over Integer Ring by the ideal (A0uE*A1uE, A0uE*A2uE, A1uE*A2uE, AEuempty + A0uE, AEuempty + A1uE, AEuempty + A2uE, AemptyuE, AemptyuE, AemptyuE)

    
    """
    
    #Setup
    Mdual= M.dual()
    if M.loops() and  M.coloops():
        raise Exception("M must be a loopless and coloopless matroid.")
    import itertools
    if R is None:
        R = ZZ
    biflat_list = biflats(M)
    E = list(M.groundset())
    
    #Find biflats containing groundset elements
    flats_containing = {x: [] for x in E}
    dual_flats_containing = {x: [] for x in E}
    for i,(F,G) in enumerate(biflat_list):
        for x in F:
            if F != M.groundset():
                flats_containing[x].append(i)
        for y in G:
            if G != Mdual.groundset():
                dual_flats_containing[y].append(i)
                
    #Find all the collecions of flats which aren't biflags.
    #We don't use `biflat_poset` because we want to keep
    #track of the index of each biflat in `biflat_list`.
    poset = Poset((list(enumerate(biflat_list)), lambda x,y : biflat_comparison(x[1],y[1])))
    biflags = set(frozenset([x[0] for x in c]) for c in poset.chains()
             if check_biflag([x[1] for x in c],M))
    nonbiflags = SimplicialComplex(biflags).minimal_nonfaces()

    # Create the ambient polynomial ring
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    try:
        names = ['AEuempty']
        names += ['A{}'.format(''.join(pretty_biflat(F,M))) for F in biflat_list]
        names += ['AemptyuE']
        P = PolynomialRing(R, names)
    except ValueError: 
        P = PolynomialRing(R, 'A', len(biflat_list)+2)
        names = P.variable_names()
        
    gens = P.gens()
    
    # Create the Stanley--Reisner relations
    I = []
    for F in nonbiflags:
        x = 1
        for i in F:
            #This plus one is there because we appended AE_empty to the front of the gens
            x = x*gens[i+1]
        I.append(x)
        
    #Create the linear relations
    #The list appendation is getting the E|empty empty|E terms
    J = [sum(gens[i] for i in [x+1 for x in flats_containing[x]] +[0])
         for x in E]
    J_bar = [sum(gens[i] for i in [x+1 for x in dual_flats_containing[x]]+[len(biflat_list)+1])
             for x in E]

    #Return the conormal ring
    return P.quotient(Ideal(I+J+J_bar))
