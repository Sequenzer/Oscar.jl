################################################################################
#
#  Subgroup function
#
################################################################################

function _as_subgroup_bare(G::T, H::GapObj) where T <: GAPGroup
  return _oscar_group(H, G)
end

function _as_subgroup(G::GAPGroup, H::GapObj)
  H1 = _as_subgroup_bare(G, H)
  return H1, hom(H1, G, x -> group_element(G, x.X), x -> group_element(H1, x.X); is_known_to_be_bijective = false)
end

"""
    sub(G::GAPGroup, gens::AbstractVector{<:GAPGroupElem}; check::Bool = true)
    sub(gens::GAPGroupElem...)

Return two objects: a group `H`, that is the subgroup of `G`
generated by the elements `x,y,...`, and the embedding homomorphism of `H`
into `G`. The object `H` has the same type of `G`, and it has no memory of the
"parent" group `G`: it is an independent group.

If `check` is set to `false` then it is not checked whether each element of
`gens` is an element of `G`.

# Examples
```jldoctest
julia> G = symmetric_group(4); H, _ = sub(G,[cperm([1,2,3]),cperm([2,3,4])]);

julia> H == alternating_group(4)
true
```
"""
function sub(G::GAPGroup, gens::AbstractVector{S}; check::Bool = true) where S <: GAPGroupElem
  @assert elem_type(G) == S
  if check
    @req all(x -> parent(x) === G || x in G, gens) "not all elements of gens lie in G"
  end
  elems_in_GAP = GapObj([x.X for x in gens])
  H = GAP.Globals.SubgroupNC(G.X, elems_in_GAP)::GapObj
  return _as_subgroup(G, H)
end

function sub(gens::GAPGroupElem...)
   @req length(gens) > 0 "Empty list"
   l = collect(gens)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]), l, check = false)
end

"""
    is_subset(H::T, G::T) where T <: GAPGroup

Return `true` if `H` is a subset of `G`, otherwise return `false`.

# Examples
```jldoctest
julia> g = symmetric_group(300); h = derived_subgroup(g)[1];

julia> is_subset(h, g)
true

julia> is_subset(g, h)
false
```
"""
function is_subset(H::T, G::T) where T <: GAPGroup
   return all(h -> h in G, gens(H))
end

"""
    is_subgroup(H::T, G::T) where T <: GAPGroup

Return (`true`,`f`) if `H` is a subgroup of `G`, where `f` is the embedding
homomorphism of `H` into `G`, otherwise return (`false`,`nothing`).

If you do not need the embedding then better call
[`is_subset(H::T, G::T) where T <: GAPGroup`](@ref).
"""
function is_subgroup(H::T, G::T) where T <: GAPGroup
   if !is_subset(H, G)
      return (false, nothing)
   else
      # We do not call `_as_subgroup` because we want to store `H`.
      return (true, hom(H, G,
                        x -> group_element(G, x.X),
                        x -> group_element(H, x.X);
                        is_known_to_be_bijective = false))
   end
end

"""
    embedding(H::T, G::T) where T <: GAPGroup

Return the embedding morphism of `H` into `G`.
An exception is thrown if `H` is not a subgroup of `G`.
"""
function embedding(H::T, G::T) where T <: GAPGroup
   a, f = is_subgroup(H, G)
   @req a "H is not a subgroup of G"
   return f
end

@doc """
    trivial_subgroup(G::GAPGroup)

Return the trivial subgroup of `G`,
together with its embedding morphism into `G`.
"""
@gapattribute trivial_subgroup(G::GAPGroup) = _as_subgroup(G, GAP.Globals.TrivialSubgroup(G.X)::GapObj)


###############################################################################
#
#  Index
#
###############################################################################

"""
    index(::Type{I} = ZZRingElem, G::T, H::T) where I <: IntegerUnion where T <: Union{GAPGroup, GrpAbFinGen}

Return the index of `H` in `G`, as an instance of `I`.
"""
index(G::T, H::T) where T <: Union{GAPGroup, GrpAbFinGen} = index(ZZRingElem, G, H)

function index(::Type{I}, G::T, H::T) where I <: IntegerUnion where T <: GAPGroup
   i = GAP.Globals.Index(G.X, H.X)::GapInt
   @req (i !== GAP.Globals.infinity) "index() not supported for subgroup of infinite index, use is_finite()"
   return I(i)
end

###############################################################################
#
#  subgroups computation
#
###############################################################################

# convert a GAP list of subgroups into a vector of Julia groups objects
function _as_subgroups(G::T, subs::GapObj) where T <: GAPGroup
  res = Vector{T}(undef, length(subs))
  for i = 1:length(res)
    res[i] = _as_subgroup_bare(G, subs[i]::GapObj)
  end
  return res
end


"""
    normal_subgroups(G::Group)

Return the vector of normal subgroups of `G` (see [`is_normal`](@ref)).
"""
@gapattribute normal_subgroups(G::GAPGroup) =
  _as_subgroups(G, GAP.Globals.NormalSubgroups(G.X))

"""
    subgroups(G::Group)

Return the vector of all subgroups of `G`.
"""
function subgroups(G::GAPGroup)
  # TODO: this is super inefficient. Slightly better would be to return an iterator
  # which iterates over the (elements of) the conjugacy classes of subgroups
  return _as_subgroups(G, GAP.Globals.AllSubgroups(G.X))
end

"""
    maximal_subgroups(G::Group)

Return the vector of maximal subgroups of `G`.
"""
@gapattribute maximal_subgroups(G::GAPGroup) =
  _as_subgroups(G, GAP.Globals.MaximalSubgroups(G.X))

"""
    maximal_normal_subgroups(G::Group)

Return the vector of maximal normal subgroups of `G`,
i.e., of those proper normal subgroups of `G` that are maximal
among the proper normal subgroups.
"""
@gapattribute maximal_normal_subgroups(G::GAPGroup) =
  _as_subgroups(G, GAP.Globals.MaximalNormalSubgroups(G.X))

"""
    minimal_normal_subgroups(G::Group)

Return the vector of minimal normal subgroups of `G`,
i.e., of those nontrivial normal subgroups of `G` that are minimal
among the nontrivial normal subgroups.
"""
@gapattribute minimal_normal_subgroups(G::GAPGroup) =
  _as_subgroups(G, GAP.Globals.MinimalNormalSubgroups(G.X))

"""
    characteristic_subgroups(G::Group)

Return the list of characteristic subgroups of `G`,
i.e., those subgroups that are invariant under all automorphisms of `G`.
"""
@gapattribute characteristic_subgroups(G::GAPGroup) =
  _as_subgroups(G, GAP.Globals.CharacteristicSubgroups(G.X))

@doc raw"""
    center(G::Group)

Return the center of `G`, i.e.,
the subgroup of all $x$ in `G` such that $x y$ equals $y x$ for every $y$
in `G`, together with its embedding morphism into `G`.
"""
@gapattribute center(G::GAPGroup) = _as_subgroup(G, GAP.Globals.Centre(G.X))

@doc raw"""
    centralizer(G::Group, H::Group)

Return the centralizer of `H` in `G`, i.e.,
the subgroup of all $g$ in `G` such that $g h$ equals $h g$ for every $h$
in `H`, together with its embedding morphism into `G`.
"""
function centralizer(G::T, H::T) where T <: GAPGroup
  return _as_subgroup(G, GAP.Globals.Centralizer(G.X, H.X))
end

@doc raw"""
    centralizer(G::Group, x::GroupElem)

Return the centralizer of `x` in `G`, i.e.,
the subgroup of all $g$ in `G` such that $g$ `x` equals `x` $g$,
together with its embedding morphism into `G`.
"""
function centralizer(G::GAPGroup, x::GAPGroupElem)
  return _as_subgroup(G, GAP.Globals.Centralizer(G.X, x.X))
end

const centraliser = centralizer

################################################################################
#
#  is_normal_subgroup, is_characteristic_subgroup, is_solvable, is_nilpotent
#
################################################################################

"""
    is_maximal_subgroup(H::T, G::T; check::Bool = true) where T <: GAPGroup

Return whether `H` is a maximal subgroup of `G`, i. e.,
whether `H` is a proper subgroup of `G` and there is no proper subgroup of `G`
that properly contains `H`.

If `check` is set to `false` then it is not checked
whether `H` is a subgroup of `G`.
If `check` is not set to `false` then an exception is thrown
if `H` is not a subgroup of `G`.

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> is_maximal_subgroup(sylow_subgroup(G, 2)[1], G)
true

julia> is_maximal_subgroup(sylow_subgroup(G, 3)[1], G)
false

julia> is_maximal_subgroup(sylow_subgroup(G, 3)[1], sylow_subgroup(G, 2)[1])
ERROR: ArgumentError: H is not a subgroup of G
```
"""
function is_maximal_subgroup(H::T, G::T; check::Bool = true) where T <: GAPGroup
  # In earlier times, `is_maximal` returned `false` if `H` was not a subgroup
  # if `G`, but at that time `G` was the first argument.
  # In order to avoid wrong results due to the reordering of arguments,
  # we throw an exception if `H` is not a subgroup of `G`.
  # (Just in case that you think about removing this exception.)
  if check
    @req is_subset(H, G) "H is not a subgroup of G"
  end
  if order(G) // order(H) < 100
    t = right_transversal(G, H)[2:end] #drop the identity
    return all(x -> order(sub(G, vcat(gens(H), [x]))[1]) == order(G), t)
  end
  return any(M -> is_conjugate(G, M, H), maximal_subgroup_reps(G))
end

"""
    is_normalized_by(H::T, G::T) where T <: GAPGroup

Return whether the group `H` is normalized by `G`, i.e.,
whether `H` is invariant under conjugation with elements of `G`.

Note that `H` need not be a subgroup of `G`.
To test whether `H` is a normal subgroup of `G`,
use [`is_normal_subgroup`](@ref).

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> is_normalized_by(sylow_subgroup(G, 2)[1], G)
false

julia> is_normalized_by(derived_subgroup(G)[1], G)
true

julia> is_normalized_by(derived_subgroup(G)[1], sylow_subgroup(G, 2)[1])
true
```
"""
is_normalized_by(H::T, G::T) where T <: GAPGroup = GAPWrap.IsNormal(G.X, H.X)

"""
    is_normal_subgroup(H::T, G::T) where T <: GAPGroup

Return whether the group `H` is a normal subgroup of `G`, i.e., whether `H`
is a subgroup of `G` that is invariant under conjugation with elements of `G`.

(See [`is_normalized_by`](@ref) for an invariance check only.)

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> is_normal_subgroup(sylow_subgroup(G, 2)[1], G)
false

julia> is_normal_subgroup(derived_subgroup(G)[1], G)
true

julia> is_normal_subgroup(derived_subgroup(G)[1], sylow_subgroup(G, 2)[1])
false
```
"""
function is_normal_subgroup(H::T, G::T) where T <: GAPGroup
  return is_subset(H, G) && is_normalized_by(H, G)
end

"""
    is_characteristic_subgroup(H::T, G::T; check::Bool = true) where T <: GAPGroup

Return whether the subgroup `H` of `G` is characteristic in `G`,
i.e., `H` is invariant under all automorphisms of `G`.

If `check` is set to `false` then it is not checked
whether `H` is a subgroup of `G`.
If `check` is not set to `false` then an exception is thrown
if `H` is not a subgroup of `G`.

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> is_characteristic_subgroup(derived_subgroup(G)[1], G)
true

julia> is_characteristic_subgroup(sylow_subgroup(G, 3)[1], G)
false

julia> is_characteristic_subgroup(sylow_subgroup(G, 3)[1], sylow_subgroup(G, 2)[1])
ERROR: ArgumentError: H is not a subgroup of G
```
"""
function is_characteristic_subgroup(H::T, G::T; check::Bool = true) where T <: GAPGroup
  if check
    @req is_subset(H, G) "H is not a subgroup of G"
  end
  return GAPWrap.IsCharacteristicSubgroup(G.X, H.X)
end

"""
    is_solvable(G::GAPGroup)

Return whether `G` is solvable,
i.e., whether [`derived_series`](@ref)(`G`)
reaches the trivial subgroup in a finite number of steps.
"""
@gapattribute is_solvable(G::GAPGroup) = GAP.Globals.IsSolvableGroup(G.X)::Bool

"""
    is_nilpotent(G::GAPGroup)

Return whether `G` is nilpotent,
i.e., whether the lower central series of `G` reaches the trivial subgroup
in a finite number of steps.
"""
@gapattribute is_nilpotent(G::GAPGroup) = GAP.Globals.IsNilpotentGroup(G.X)::Bool

"""
    is_supersolvable(G::GAPGroup)

Return whether `G` is supersolvable,
i.e., `G` is finite and has a normal series with cyclic factors.
"""
@gapattribute is_supersolvable(G::GAPGroup) = GAP.Globals.IsSupersolvableGroup(G.X)::Bool

################################################################################
#
#  Quotient functions
#
################################################################################

function quo(G::FPGroup, elements::Vector{S}) where S <: GAPGroupElem
  @assert elem_type(G) == S
  if GAP.Globals.HasIsWholeFamily(G.X) && GAPWrap.IsWholeFamily(G.X)
    # For a *full* free or f.p. group, GAP can handle this via its `\/'.
    elems_in_gap = GapObj([x.X for x in elements])
    Q = FPGroup((G.X)/elems_in_gap)
    function proj(x::FPGroupElem)
      return group_element(Q,GAP.Globals.MappedWord(x.X,
               GAPWrap.GeneratorsOfGroup(G.X), GAPWrap.GeneratorsOfGroup(Q.X)))
    end
    return Q, hom(G,Q,proj)
  else
    # Currently GAP's `\/' does not support a list of group elements
    # as the second argument,
    # but forming the quotient modulo a normal subgroup may work.
    return quo(G, normal_closure(G, sub(G, elements)[1])[1])
  end
end

"""
    quo([::Type{Q}, ]G::T, elements::Vector{elem_type(G)})) where {Q <: GAPGroup, T <: GAPGroup}

Return the quotient group `G/N`, together with the projection `G` -> `G/N`,
where `N` is the normal closure of `elements` in `G`.

See [`quo(G::T, N::T) where T <: GAPGroup`](@ref)
for information about the type of `G/N`.
"""
function quo(G::T, elements::Vector{S}) where T <: GAPGroup where S <: GAPGroupElem
  @assert elem_type(G) == S
  if length(elements) == 0
    H1 = trivial_subgroup(G)[1]
  else
    elems_in_gap = GapObj([x.X for x in elements])
    H = GAP.Globals.NormalClosure(G.X,GAP.Globals.Group(elems_in_gap))::GapObj
    @assert GAPWrap.IsNormal(G.X, H)
    H1 = _as_subgroup_bare(G, H)
  end
  return quo(G, H1)
end

function quo(::Type{Q}, G::T, elements::Vector{S}) where {Q <: GAPGroup, T <: GAPGroup, S <: GAPGroupElem}
  F, epi = quo(G, elements)
  if !(F isa Q)
    map = isomorphism(Q, F)
    F = codomain(map)
    epi = compose(epi, map)
  end
  return F, epi
end

"""
    quo([::Type{Q}, ]G::T, N::T) where {Q <: GAPGroup, T <: GAPGroup}

Return the quotient group `G/N`, together with the projection `G` -> `G/N`.

If `Q` is given then `G/N` has type `Q` if possible,
and an exception is thrown if not.

If `Q` is not given then the type of `G/N` is not determined by the type of `G`.
- `G/N` may have the same type as `G` (which is reasonable if `N` is trivial),
- `G/N` may have type `PcGroup` (which is reasonable if `G/N` is finite and solvable), or
- `G/N` may have type `PermGroup` (which is reasonable if `G/N` is finite and non-solvable).
- `G/N` may have type `FPGroup` (which is reasonable if `G/N` is infinite).

An exception is thrown if `N` is not a normal subgroup of `G`.

# Examples
```jldoctest
julia> G = symmetric_group(4)
Permutation group of degree 4 and order 24

julia> N = pcore(G, 2)[1];

julia> typeof(quo(G, N)[1])
PcGroup

julia> typeof(quo(PermGroup, G, N)[1])
PermGroup
```
"""
function quo(G::T, N::T) where T <: GAPGroup
  mp = GAP.Globals.NaturalHomomorphismByNormalSubgroup(G.X, N.X)::GapObj
  # The call may have found out new information about `G.X`,
  # for example that `G.X` is finite.
#FIXME: The GAP function should deal with this situation.
  GAP.Globals.UseSubsetRelation(G.X, N.X)
  cod = GAP.Globals.ImagesSource(mp)::GapObj
  S = elem_type(G)
  S1 = _get_type(cod)
  codom = S1(cod)
  mp_julia = __create_fun(mp, codom, S)
  return codom, hom(G, codom, mp_julia)
end

function quo(::Type{Q}, G::T, N::T) where {Q <: GAPGroup, T <: GAPGroup}
  F, epi = quo(G, N)
  if !(F isa Q)
    map = isomorphism(Q, F)
    F = codomain(map)
    epi = compose(epi, map)
  end
  return F, epi
end

"""
    maximal_abelian_quotient([::Type{Q}, ]G::GAPGroup) where Q <: Union{GAPGroup, GrpAbFinGen}

Return `F, epi` such that `F` is the largest abelian factor group of `G`
and `epi` is an epimorphism from `G` to `F`.

If `Q` is given then `F` has type `Q` if possible,
and an exception is thrown if not.

If `Q` is not given then the type of `F` is not determined by the type of `G`.
- `F` may have the same type as `G` (which is reasonable if `G` is abelian),
- `F` may have type `PcGroup` (which is reasonable if `F` is finite), or
- `F` may have type `FPGroup` (which is reasonable if `F` is infinite).

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> F, epi = maximal_abelian_quotient(G);

julia> order(F)
2

julia> domain(epi) === G && codomain(epi) === F
true

julia> typeof(F)
PcGroup

julia> typeof(maximal_abelian_quotient(free_group(1))[1])
FPGroup

julia> typeof(maximal_abelian_quotient(PermGroup, G)[1])
PermGroup
```
"""
function maximal_abelian_quotient(G::GAPGroup)
  map = GAP.Globals.MaximalAbelianQuotient(G.X)::GapObj
  F = GAPWrap.Range(map)::GapObj
  S1 = _get_type(F)
  F = S1(F)
  return F, GAPGroupHomomorphism(G, F, map)
end

function maximal_abelian_quotient(::Type{Q}, G::GAPGroup) where Q <: Union{GAPGroup, GrpAbFinGen}
  F, epi = maximal_abelian_quotient(G)
  if !(F isa Q)
    map = isomorphism(Q, F)
    F = codomain(map)
    epi = compose(epi, map)
  end
  return F, epi
end

has_maximal_abelian_quotient(G::GAPGroup) = GAPWrap.HasMaximalAbelianQuotient(G.X)

function set_maximal_abelian_quotient(G::T, val::Tuple{GAPGroup, GAPGroupHomomorphism{T}}) where T <: GAPGroup
  return GAPWrap.SetMaximalAbelianQuotient(G.X, val[2].map)
end


function __create_fun(mp, codom, ::Type{S}) where S
  function mp_julia(x::S)
    el = GAPWrap.Image(mp, x.X)
    return group_element(codom, el)
  end
  return mp_julia
end

"""
    epimorphism_from_free_group(G::GAPGroup)

Return an epimorphism `epi` from a free group `F == domain(epi)` onto `G`,
where `F` has the same number of generators as `G` and such that for each `i`
it maps `gen(F,i)` to `gen(G,i)`.

A useful application of this function is expressing an element of `G` as
a word in its generators.

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> epi = epimorphism_from_free_group(G)
Group homomorphism
  from free group
  to permutation group of degree 4 and order 24

julia> pi = G([2,4,3,1])
(1,2,4)

julia> w = preimage(epi, pi);

julia> map_word(w, gens(G))
(1,2,4)
```
"""
function epimorphism_from_free_group(G::GAPGroup)
  mfG = GAP.Globals.EpimorphismFromFreeGroup(G.X)
  fG = FPGroup(GAPWrap.Source(mfG))
  return Oscar.GAPGroupHomomorphism(fG, G, mfG)
end

################################################################################
#
#  Derived subgroup and derived series
#
################################################################################

"""
    derived_subgroup(G::GAPGroup)

Return the derived subgroup of `G`, i.e.,
the subgroup generated by all commutators of `G`.
"""
@gapattribute derived_subgroup(G::GAPGroup) =
  _as_subgroup(G, GAP.Globals.DerivedSubgroup(G.X))

@doc raw"""
    derived_series(G::GAPGroup)

Return the vector $[ G_1, G_2, \ldots ]$,
where $G_1 =$ `G` and $G_{i+1} =$ `derived_subgroup`$(G_i)$.
"""
@gapattribute derived_series(G::GAPGroup) = _as_subgroups(G, GAP.Globals.DerivedSeries(G.X))


################################################################################
#
#  Intersection
#
################################################################################

@doc raw"""
    intersect(V::T...) where T <: Group
    intersect(V::AbstractVector{T}) where T <: Group

If `V` is $[ G_1, G_2, \ldots, G_n ]$,
return the intersection $K$ of the groups $G_1, G_2, \ldots, G_n$,
together with the embeddings of $K into $G_i$.
"""
function intersect(V::T...) where T<:GAPGroup
   L = GapObj([G.X for G in V])
   K = GAP.Globals.Intersection(L)::GapObj
   Embds = [_as_subgroup(G, K)[2] for G in V]
   K = _as_subgroup(V[1], K)[1]
   Arr = Tuple(vcat([K],Embds))
   return Arr
end

function intersect(V::AbstractVector{T}) where T<:GAPGroup
   L = GapObj([G.X for G in V])
   K = GAP.Globals.Intersection(L)::GapObj
   Embds = [_as_subgroup(G, K)[2] for G in V]
   K = _as_subgroup(V[1], K)[1]
   Arr = Tuple(vcat([K],Embds))
   return Arr
end
#T why duplicate this code?
