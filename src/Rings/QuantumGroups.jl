
export QuantumGroup, quantum_group, iszero, quantum_symmetric_group


"""
    QuantumGroup
```jldoctest
free, (x,y,z) = free_associative_algebra(QQ, ["x", "y", "z"]);
f1 = x*y + y*z;
I = ideal([f1]);

Q = quantum_group(I);
iszero(Q, f1, 4)

```
"""
mutable struct QuantumGroup 
  base_ring::FreeAssAlgebra
  relations::FreeAssAlgIdeal

  function QuantumGroup(R::FreeAssAlgebra, relations::FreeAssAlgIdeal)
    r = new()
    r.base_ring = R
    r.relations = relations
    return r
  end
end

Base.show(io::IO, Q::QuantumGroup) = print(io, "QuantumGroup")


quantum_group(R::FreeAssAlgebra, relations::FreeAssAlgIdeal) = QuantumGroup(R, relations)
quantum_group(relations::FreeAssAlgIdeal) = QuantumGroup(base_ring(relations), relations)
function quantum_group(I::Vector{<:FreeAssAlgElem})
  R = parent(I[1])
  @req all(x -> parent(x) == R, I) "parent mismatch"
  return QuantumGroup(parent(I[1]), relations(I))
end

base_ring(Q::QuantumGroup) = Q.base_ring
relations(Q::QuantumGroup) = Q.relations

function iszero(Q::QuantumGroup, f::FreeAssAlgElem, deg_bound::Int=-1)
  return ideal_membership(f, relations(Q), deg_bound)
end

@doc raw"""
    getQuantumPermutationGroup(n::Int,interreduce::Bool=true)

Get the relations that define the quantum permutation group on `n` elements. If interreduce is true, it uses the interreduce function from Oscar to reduce the number of generators.

# Examples
```jldoctest
id = quantum_symmetric_group(3)
gens(relations(id))
# output

20
```
"""
function quantum_symmetric_group(n::Int)
    generator_strings = String[]
    for i in 1:n, j in 1:n
            push!(generator_strings, "u[$i,$j]")
    end
    A, g = free_associative_algebra(Oscar.QQ, generator_strings)
    u = Matrix{elem_type(A)}(undef, n, n)
    for i in 1:n, j in 1:n
            u[i, j] = g[(i-1)*n + j]
    end
    relations = elem_type(A)[]
    #Squared relations
    for i in 1:n, j in 1:n
        new_relation = u[i, j] * u[i, j] - u[i, j]
            push!(relations, new_relation)
            for k in 1:n
                    if k != j
                        new_relation = u[i,j] * u[i, k]
                        push!(relations, new_relation)
                        new_relation = u[j, i]*u[k, i]
                        push!(relations, new_relation)
                    end
            end
    end

    #row and column sum relations
    for i in 1:n
        new_relation_row = -1
        new_relation_col = -1
        for k in 1:n
            new_relation_row += u[i,k]
            new_relation_col += u[k,i]
        end
        push!(relations, new_relation_row)
        push!(relations, new_relation_col)

    end

    return quantum_group(ideal(Vector{elem_type(A)}(relations)))
end

