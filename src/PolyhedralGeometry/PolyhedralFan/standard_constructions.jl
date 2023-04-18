#TODO: inward/outward options? via polymake changes?

"""
    normal_fan(P::Polyhedron)

Return the normal fan of `P`. The maximal cones of the normal fan of `P` are
dual to the edge cones at the vertices of `P`.

# Examples
The rays of a normal fan of a cube point in every positive and negative unit
direction.
```jldoctest
julia> C = cube(3);

julia> NF = normal_fan(C)
Polyhedral fan in ambient dimension 3

julia> rays(NF)
6-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0, 0]
 [-1, 0, 0]
 [0, 1, 0]
 [0, -1, 0]
 [0, 0, 1]
 [0, 0, -1]
```
"""
function normal_fan(P::Polyhedron{T}) where T<:scalar_types
   pmp = pm_object(P)
   pmnf = Polymake.fan.normal_fan(pmp)
   return PolyhedralFan{T}(pmnf)
end

"""
    face_fan(P::Polyhedron)

Return the face fan of `P`. The polytope `P` has to contain the origin, then
the maximal cones of the face fan of `P` are the cones over the facets of `P`.

# Examples
By definition, this bounded polyhedron's number of facets equals the amount of
maximal cones of its face fan.
```jldoctest
julia> C = cross_polytope(3);

julia> FF = face_fan(C)
Polyhedral fan in ambient dimension 3

julia> n_maximal_cones(FF) == nfacets(C)
true
```
"""
function face_fan(P::Polyhedron{T}) where T<:scalar_types
   pmp = pm_object(P)
   pmff = Polymake.fan.face_fan(pmp)
   return PolyhedralFan{T}(pmff)
end


###############################################################################
## Star subdivision
###############################################################################

@doc raw"""
    star_subdivision(PF::PolyhedralFan, n::Int)

Return the star subdivision of a polyhedral fan at its n-th torus orbit.
Note that this torus orbit need not be maximal. We follow definition 3.3.17
of [CLS11](@cite).

# Examples
```jldoctest
julia> star = star_subdivision(normal_fan(simplex(3)), 1)
Polyhedral fan in ambient dimension 3

julia> rays(star)
5-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 0, 0]
 [0, 0, 1]
 [1, 1, 1]
 [0, 1, 0]
 [-1, -1, -1]

julia> ray_indices(maximal_cones(star))
6×5 IncidenceMatrix
[1, 2, 3]
[2, 3, 4]
[1, 3, 4]
[2, 4, 5]
[1, 2, 5]
[1, 4, 5]
```
"""
function star_subdivision(PF::_FanLikeType{T}, n::Int) where T<:scalar_types
  
  # check if n-th cone exist
  @req n <= n_cones(PF) "Cannot subdivide cone $n as it does not exist"
  
  # check if n-th cone can be subdivided
  cone_list = cones(PF)
  nthcone = Polymake.row(cone_list, n)
  @req length(nthcone) > 1 "Cannot subdivide cone $n as it is generated by a single ray"
  
  # check if nth cone is smooth
  @req _cone_is_smooth(PF, nthcone) "Cannot subdivide maximal cone $n as it is not smooth"
  
  # compute new ray
  R = Polymake.common.primitive(Oscar.pm_object(PF).RAYS)
  newindex = size(R,1) + 1
  newray = sum([R[i,:] for i in nthcone])
  
  # compute new rays
  newrays = [R; transpose(newray)]
  
  # compute the new maximal cones
  max_cone_indices = ray_indices(maximal_cones(PF))
  cone_indices = Polymake.row(cone_list, n)
  newmaxcones = (Vector{Int})[]
  for i in 1:n_maximal_cones(PF)
    indices_ith_max_cone = Polymake.row(max_cone_indices, i)
    if issubset(cone_indices, indices_ith_max_cone)
      @req _cone_is_smooth(PF, indices_ith_max_cone) "All cones containing sigma need to be smooth"
      for subset in subsets(Vector{Int64}(cone_indices), length(cone_indices)-1)
        tmp = Base.setdiff(indices_ith_max_cone, cone_indices)
        tmp = Base.union(subset, tmp)
        push!(tmp, newindex)
        push!(newmaxcones, tmp)
      end
    else
      push!(newmaxcones, Vector{Int}(indices_ith_max_cone))
    end
  end
  newmaxcones = IncidenceMatrix(newmaxcones)
  
  # return the new fan
  return PolyhedralFan{T}(newrays, newmaxcones)
  
end


function _cone_is_smooth(PF::_FanLikeType, c::AbstractSet{<:Integer})
  R = matrix(ZZ, rays(PF))
  return _is_unimodular(R[Vector{Int}(c),:])
end

function _is_unimodular(M::ZZMatrix)
  nrows(M) <= ncols(M) || return false
  n = nrows(M)
  return abs(det(snf(M)[:,1:n])) == 1
end


###############################################################################
## Cartesian/Direct product
###############################################################################

@doc raw"""
    *(PF1::PolyhedralFan, PF2::PolyhedralFan)

Return the Cartesian/direct product of two polyhedral fans.

# Examples
```jldoctest
julia> normal_fan(simplex(2))*normal_fan(simplex(3))
Polyhedral fan in ambient dimension 5
```
"""
function Base.:*(PF1::PolyhedralFan, PF2::PolyhedralFan)
    prod = Polymake.fan.product(pm_object(PF1), pm_object(PF2))
    return PolyhedralFan{detect_scalar_type(PolyhedralFan, prod)}(prod)
end
