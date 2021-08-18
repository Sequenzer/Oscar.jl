import Base.intersect
import GAP.@gapattribute

export
    centralizer,
    centre, hascentre, setcentre,
    characteristic_subgroups, hascharacteristic_subgroups, setcharacteristic_subgroups,
    derived_series, hasderived_series, setderived_series,
    derived_subgroup, hasderived_subgroup, setderived_subgroup,
    embedding,
    index,
    ischaracteristic,
    isnilpotent, hasisnilpotent, setisnilpotent,
    issolvable, hasissolvable, setissolvable,
    issupersolvable, hasissupersolvable, setissupersolvable,
    maximal_normal_subgroups, hasmaximal_normal_subgroups, setmaximal_normal_subgroups,
    maximal_subgroups, hasmaximal_subgroups, setmaximal_subgroups,
    minimal_normal_subgroups, hasminimal_normal_subgroups, setminimal_normal_subgroups,
    normal_subgroups, hasnormal_subgroups, setnormal_subgroups,
    quo,
    sub,
    trivial_subgroup, hastrivial_subgroup, settrivial_subgroup

################################################################################
#
#  Subgroup function
#
################################################################################

function _as_subgroup_bare(G::T, H::GapObj) where T
  if T==PermGroup
    H1 = T(H, G.deg)
  elseif T<:MatrixGroup
    H1 = MatrixGroup(G.deg,G.ring)
    H1.ring_iso = G.ring_iso
    H1.mat_iso = G.mat_iso
    H1.X = H
  else
    H1 = T(H)
  end
  return H1
end

function _as_subgroup(G::T, H::GapObj, ::Type{S}) where { T, S }
  H1 = _as_subgroup_bare(G, H)
  return H1, hom(H1, G, x::S -> group_element(G, x.X))
end

function _as_subgroup(G::T, H::GapObj) where T <: GAPGroup
  return _as_subgroup(G, H, elem_type(G))
end

function sub(G::T, elements::Vector{S}) where T <: GAPGroup where S <: GAPGroupElem
  @assert elem_type(G) == S
  elems_in_GAP = GAP.julia_to_gap(GapObj[x.X for x in elements])
  H = GAP.Globals.Subgroup(G.X,elems_in_GAP)
  #H is the group. I need to return the inclusion map too
  return _as_subgroup(G, H)
end

function sub(L::GAPGroupElem...)
   if length(L)==0 throw(ArgumentError("Empty list")) end
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

"""
    issubgroup(G::T, H::T) where T <: GAPGroup

Return (`true`,`f`) if `H` is a subgroup of `G`, where `f` is the embedding homomorphism of `H` into `G`, otherwise return (`false`,`nothing`).
"""
function issubgroup(G::T, H::T) where T <: GAPGroup
   if !all(h -> h in G, gens(H))
      return (false, nothing)
   else
      return (true, _as_subgroup(G, H.X)[2])
   end
end

"""
    embedding(G::T, H::T) where T <: GAPGroup

Return the embedding morphism of `H` into `G`.
An exception is thrown if `H` is not a subgroup of `G`.
"""
function embedding(G::T, H::T) where T <: GAPGroup
   a,f = issubgroup(G,H)
   if !a
      throw(ArgumentError("H is not a subgroup of G"))
   else
      return f
   end
end

@gapattribute trivial_subgroup(G::GAPGroup) =_as_subgroup(G, GAP.Globals.TrivialSubgroup(G.X))
@doc """
    trivial_subgroup(G::GAPGroup)

Return the trivial subgroup of `G`,
together with its embedding morphism into `G`.
""" trivial_subgroup


###############################################################################
#
#  Index
#
###############################################################################

"""
    index(::Type{I} = fmpz, G::T, H::T) where I <: Union{Integer, fmpz} where T <: GAPGroup

Return the index of `H` in `G`, as an instance of `I`.
"""
index(G::T, H::T) where T <: GAPGroup = index(fmpz, G, H)

function index(::Type{I}, G::T, H::T) where I <: Union{Base.Integer, fmpz} where T <: GAPGroup
   i = GAP.Globals.Index(G.X, H.X)
   if i === GAP.Globals.infinity
      error("index() not supported for subgroup of infinite index, use isfinite()")
   end
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
    res[i] = _as_subgroup_bare(G, subs[i])
  end
  return res
end


@gapattribute normal_subgroups(G::GAPGroup) =
  _as_subgroups(G, GAP.Globals.NormalSubgroups(G.X))
"""
    normal_subgroups(G::Group)

Return the vector of normal subgroups of `G` (see [`isnormal`](@ref)).
""" normal_subgroups

"""
    subgroups(G::Group)

Return the vector of all subgroups of `G`.
"""
function subgroups(G::GAPGroup)
  return _as_subgroups(G, GAP.Globals.AllSubgroups(G.X))
end

@gapattribute maximal_subgroups(G::GAPGroup) =
  _as_subgroups(G, GAP.Globals.MaximalSubgroups(G.X))
"""
    maximal_subgroups(G::Group)

Return the vector of maximal subgroups of `G`.
""" maximal_subgroups

@gapattribute maximal_normal_subgroups(G::GAPGroup) =
  _as_subgroups(G, GAP.Globals.MaximalNormalSubgroups(G.X))
"""
    maximal_normal_subgroups(G::Group)

Return the vector of maximal normal subgroups of `G`,
i. e., of those proper normal subgroups of `G` that are maximal
among the proper normal subgroups.
""" maximal_normal_subgroups

@gapattribute minimal_normal_subgroups(G::GAPGroup) =
  _as_subgroups(G, GAP.Globals.MinimalNormalSubgroups(G.X))
"""
    minimal_normal_subgroups(G::Group)

Return the vector of minimal normal subgroups of `G`,
i. e., of those nontrivial normal subgroups of `G` that are minimal
among the nontrivial normal subgroups.
""" minimal_normal_subgroups

@gapattribute characteristic_subgroups(G::GAPGroup) =
  _as_subgroups(G, GAP.Globals.CharacteristicSubgroups(G.X))
"""
    characteristic_subgroups(G::Group)

Return the list of characteristic subgroups of `G`,
i.e., those subgroups that are invariant under all automorphisms of `G`.
""" characteristic_subgroups

@gapattribute centre(G::GAPGroup) = _as_subgroup(G, GAP.Globals.Centre(G.X))
@doc Markdown.doc"""
    centre(G::Group)

Return the centre of `G`, i.e.,
the subgroup of all $x$ in `G` such that $x y$ equals $y x$ for every $y$
in `G`, together with its embedding morphism into `G`.
""" centre

@doc Markdown.doc"""
    centralizer(G::Group, H::Group)

Return the centralizer of `H` in `G`, i.e.,
the subgroup of all $g$ in `G` such that $g h$ equals $h g$ for every $h$
in `H`, together with its embedding morphism into `G`.
"""
function centralizer(G::T, H::T) where T <: GAPGroup
  return _as_subgroup(G, GAP.Globals.Centralizer(G.X, H.X))
end

@doc Markdown.doc"""
    centralizer(G::Group, x::GroupElem) 

Return the centralizer of `x` in `G`, i.e.,
the subgroup of all $g$ in `G` such that $g$ `x` equals `x` $g$,
together with its embedding morphism into `G`.
"""
function centralizer(G::GAPGroup, x::GAPGroupElem)
  return _as_subgroup(G, GAP.Globals.Centralizer(G.X, x.X))
end

centraliser = centralizer

################################################################################
#
#  IsNormal, IsCharacteristic, IsSolvable, IsNilpotent
#
################################################################################

"""
    isnormal(G::T, H::T) where T <: GAPGroup

Return whether the subgroup `H` is normal in `G`,
i. e., `H` is invariant under conjugation with elements of `G`.
"""
isnormal(G::T, H::T) where T <: GAPGroup = GAP.Globals.IsNormal(G.X, H.X)

"""
    ischaracteristic(G::T, H::T) where T <: GAPGroup

Return whether the subgroup `H` is characteristic in `G`,
i. e., `H` is invariant under all automorphisms of `G`.
"""
function ischaracteristic(G::T, H::T) where T <: GAPGroup
  return GAP.Globals.IsCharacteristicSubgroup(G.X, H.X)
end

@gapattribute issolvable(G::GAPGroup) = GAP.Globals.IsSolvableGroup(G.X)::Bool
"""
    issolvable(G::GAPGroup)

Return whether `G` is solvable,
i. e., whether [`derived_series`](@ref)(`G`)
reaches the trivial subgroup in a finite number of steps.
""" issolvable

@gapattribute isnilpotent(G::GAPGroup) = GAP.Globals.IsNilpotentGroup(G.X)::Bool
"""
    isnilpotent(G::GAPGroup)

Return whether `G` is nilpotent,
i. e., whether the lower central series of `G` reaches the trivial subgroup
in a finite number of steps.
""" isnilpotent

@gapattribute issupersolvable(G::GAPGroup) = GAP.Globals.IsSupersolvableGroup(G.X)::Bool
"""
    issupersolvable(G::GAPGroup)

Return whether `G` is supersolvable,
i. e., `G` is finite and has a normal series with cyclic factors.
""" issupersolvable

################################################################################
#
#  Quotient function
#
################################################################################

function quo(G::FPGroup, elements::Vector{S}) where T <: GAPGroup where S <: GAPGroupElem
  @assert elem_type(G) == S
  elems_in_gap = GAP.julia_to_gap(GapObj[x.X for x in elements])
#T better!
  Q=FPGroup((G.X)/elems_in_gap)
  function proj(x::FPGroupElem)
     return group_element(Q,GAP.Globals.MappedWord(x.X,GAP.Globals.GeneratorsOfGroup(G.X), GAP.Globals.GeneratorsOfGroup(Q.X)))
  end
  return Q, hom(G,Q,proj)
end

"""
    quo(G::T, elements::Vector{S})

Return the quotient group `G/H` of type `FPGroup` (if `T`=`FPGroup`), `PcGroup` (if the quotient group is solvable) or `PermGroup` (otherwise), where `H` is the normal closure of `elements` in `G`.
"""
function quo(G::T, elements::Vector{S}) where T <: GAPGroup where S <: GAPGroupElem
  @assert elem_type(G) == S
  elems_in_gap = GAP.julia_to_gap(GapObj[x.X for x in elements])
  H = GAP.Globals.NormalClosure(G.X,GAP.Globals.Group(elems_in_gap))
  @assert GAP.Globals.IsNormal(G.X, H)
  H1 = T(H)
  return quo(G, H1)
end
#T prescribe other type?

"""
    quo(G::T, H::T)

Return the quotient group `G/H` of type `PcGroup` (if the quotient group is solvable) or `PermGroup` (otherwise), together with the projection `G` -> `G/H`.
"""
function quo(G::T, H::T) where T <: GAPGroup
  mp = GAP.Globals.NaturalHomomorphismByNormalSubgroup(G.X, H.X)
  cod = GAP.Globals.ImagesSource(mp)
  S = elem_type(G)
  S1 = _get_type(cod)
  codom = S1(cod)
  mp_julia = __create_fun(mp, codom, S)
  return codom, hom(G, codom, mp_julia)
end

function __create_fun(mp, codom, ::Type{S}) where S
  function mp_julia(x::S)
    el = GAP.Globals.Image(mp, x.X)
    return group_element(codom, el)
  end
  return mp_julia
end

################################################################################
#
#  Derived subgroup and derived series
#  
################################################################################

@gapattribute derived_subgroup(G::GAPGroup) =
  _as_subgroup(G, GAP.Globals.DerivedSubgroup(G.X))
"""
    derived_subgroup(G::GAPGroup)

Return the derived subgroup of `G`, i.e.,
the subgroup generated by all commutators of `G`.
""" derived_subgroup

@gapattribute derived_series(G::GAPGroup) = _as_subgroups(G, GAP.Globals.DerivedSeries(G.X))
@doc Markdown.doc"""
    derived_series(G::GAPGroup)

Return the vector $[ G_1, G_2, \ldots ]$,
where $G_1 =$ `G` and $G_{i+1} =$ `derived_subgroup`$(G_i)$.
""" derived_series


################################################################################
#
#  Intersection
#
################################################################################

@doc Markdown.doc"""
    intersect(V::T...) where T <: Group
    intersect(V::AbstractVector{T}) where T <: Group

If `V` is $[ G_1, G_2, \ldots, G_n ]$,
return the intersection $K$ of the groups $G_1, G_2, \ldots, G_n$,
together with the embeddings of $K into $G_i$.
"""
function intersect(V::T...) where T<:GAPGroup
   L = GAP.julia_to_gap([G.X for G in V])
   K = GAP.Globals.Intersection(L)
   Embds = [_as_subgroup(G, K)[2] for G in V]
   K = _as_subgroup(V[1], K)[1]
   Arr = Tuple(vcat([K],Embds))
   return Arr
end

function intersect(V::AbstractVector{T}) where T<:GAPGroup
   L = GAP.julia_to_gap([G.X for G in V])
   K = GAP.Globals.Intersection(L)
   Embds = [_as_subgroup(G, K)[2] for G in V]
   K = _as_subgroup(V[1], K)[1]
   Arr = Tuple(vcat([K],Embds))
   return Arr
end
#T why duplicate this code?


################################################################################
#
#  Conversions between types
#
################################################################################

_get_iso_function(::Type{PermGroup}) = GAP.Globals.IsomorphismPermGroup
_get_iso_function(::Type{FPGroup}) = GAP.Globals.IsomorphismFpGroup
_get_iso_function(::Type{PcGroup}) = GAP.Globals.IsomorphismPcGroup

function isomorphic_group(::Type{T}, G::GAPGroup) where T <: GAPGroup
  f = _get_iso_function(T)
  mp = f(G.X)
  G1 = T(GAP.Globals.ImagesSource(mp))
  fmap = _hom_from_gap_map(G, G1, mp)
  return G1, fmap
end
