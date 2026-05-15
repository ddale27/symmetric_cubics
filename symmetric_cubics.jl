# To deal with characters other than the trivial one
using Oscar
using LinearAlgebra

# Groups we are using to act
S5 = symmetric_group(5)
groups = map(representative,conjugacy_classes_subgroups(S5))

# First get the elements of the unit circle that the group can even map into
orders = (Int∘order∘representative).(conjugacy_classes(S5))
#roots_of_unity = [GAP.Globals.E(i) for i in orders]
#roots_of_unity = GAP.Obj(roots_of_unity)
#roots_of_unity = GAP.Globals.List(roots_of_unity, l->GAP.Obj([GAP.Obj([l])]))
#C = GAP.Globals.Group(g)
C = cyclic_group(lcm(orders))

# Make rings
#R = QQ
R, ζ = cyclotomic_field(lcm(orders))
S, (x,y,z,w) = polynomial_ring(R,[:x,:y,:z,:w])
# Sylvester ring -- to calculate the space of invariant cubics
Sylv, v = polynomial_ring(R, :v => 0:4)

M3 = monomials_of_degree(Sylv,3)
φ = hom(Sylv, S, [x,y,z,w,-sum([x,y,z,w])])

for (j,G) in enumerate(groups)
    group_id = GAP.Globals.IdGroup(GAP.Obj(G))

    # Make directory for group
    groupdir = string("S5n", j)
    try
        mkdir(string(groupdir))
    catch
    end

    # Put a description of the group in this directory
    idfile = open(string(groupdir, "/", "group_id"), "w")
    write(idfile, string(group_id, "\n"))
    write(idfile, "Generators:\n")
    [write(idfile, string(gen, "\n")) for gen in gens(G)]
    close(idfile)

    ## Calculate orbits
    #orbs = orbits(gset(G, on_indeterminates, M3))

    # Get all characters of G
    characters = GAP.Globals.AllHomomorphisms(GAP.Obj(G),GAP.Obj(C))

    # Calculate all eigenspaces
    eigenspaces = Vector{Vector{typeof(x)}}([])

    for chi in characters
        # Just lazily doing discrete log to figure out exponents of the
        # cyclotomics
        weights = g->ζ^GAP.Globals.DLog(GAP.Obj(C[1]), chi(GAP.Obj(g)))

        # Take the associated weighted sums
        weighted_sums = Vector{typeof(x)}([])
        for f in collect(M3)
            weighted_sum = sum([weights(g)^-1*on_indeterminates(f,g) for g in collect(G)])
            push!(weighted_sums, weighted_sum)
        end

        push!(eigenspaces, weighted_sums)
    end
    # Do linear algebra to figure out bases
    V = [vector_space(R, space) for space in eigenspaces]
    B = [map(linmap, basis(space)) for (space, linmap) in V]

    # Clear denominators to avoid problems later
    denoms = Vector{ZZRingElem}([])
    for space in B
        for f in space
            append!(denoms, denominator.(coefficients(f)))
        end
    end
    B = [lcm(denoms).*space for space in B]

    # Write Sylvester system to file
    sylvfile = open(string(groupdir, "/", "sylvester"), "w")
    for space in B
        write(sylvfile, string("Eigenspace of dim ", length(space), "\n"))

        for f in space
            write(sylvfile, string(f, "\n"))
        end
    end
    close(sylvfile)


    # Now we can finally get the ideals we need
    I = [φ(ideal(if isempty(space) Sylv(0) else space end)) for space in B]
    # Do linear algebra to figure out a basis
    V = [vector_space(R,gens(i)) for i in I]
    B = [map(linmap, basis(space)) for (space, linmap) in V]

    # Clear denominators to avoid problems later
    denoms = Vector{ZZRingElem}([])
    for space in B
        for f in space
            append!(denoms, denominator.(coefficients(f)))
        end
    end
    B = [lcm(denoms).*space for space in B]

    # Write standard system to file
    stanfile = open(string(groupdir, "/", "standard"), "w")
    for space in B
        write(stanfile, string("Eigenspace of dim ", length(space), "\n"))

        for f in space
            write(stanfile, string(f, "\n"))
        end
    end
    close(stanfile)
end
