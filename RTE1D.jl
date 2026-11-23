using Ferrite
using Plots
using SparseArrays
using Test

σ(x::Float64) = 1
σₛ(x::Float64) = 0.

ncell_x = 64
ncell_μ = 32


###############################################################################
###############################################################################
# Get the defined facet tags in a grid
function getfacettags(grid::Grid)
    collect( keys( Ferrite.getfacetsets(grid) ) )
end

# Get the facet IDs and the corresponding cell IDs for "facets" in a given facet tag
function getfacetsetids(grid::Grid, tag)
    first( getfacetset(grid, tag) ).idx     # returns(cellid, facetid)
end

# Get the IDs of cells on the boundary 
function boundary_cellids(grid::Grid)
    Set( union(getfacetsetids(grid, "left")[1], getfacetsetids(grid, "right")[1]) )
end

"""
    Get the coordinate of the boundary node in the given 
    boundary facet (specified by `facetidx`) in `grid`
"""
function get_boundary_point(grid::Grid, facetidx::FacetIndex)
    cell_id, facet_id = facetidx.idx
    cell = getcells(grid, cell_id)
    cell_facets = Ferrite.facets(cell)  # facets given in node ids
    boundary_node_id = cell_facets[facet_id][1]
    (boundary_node_id, get_node_coordinate(grid, boundary_node_id)[1])
end


"""
    Get the `angular cells IDs` corresponding to inflow at a given boundary node coordinate
"""
function inflow_angular_dofs1D(boundary_node_coord)
    inflow_angular_dofs = []

    if boundary_node_coord ≈ LEFT
        for ang_cell in CellIterator(dh_μ)
            if all( x-> x >= 0., Iterators.flatten(getcoordinates(ang_cell)))
                # @show celldofs(ang_cell)
                append!(inflow_angular_dofs, celldofs(ang_cell))
            end
        end
    elseif boundary_node_coord ≈ RIGHT
        for ang_cell in CellIterator(dh_μ)
            if all( x-> x <= 0., Iterators.flatten(getcoordinates(ang_cell)))
                append!(inflow_angular_dofs, celldofs(ang_cell))
            end
        end
    end

    return inflow_angular_dofs
end

function outflow_angular_dofs1D(boundary_node_coord)
    outflow_angular_dofs = []

    if boundary_node_coord ≈ LEFT
        for ang_cell in CellIterator(dh_μ)
            if all( x-> x <= 0., Iterators.flatten(getcoordinates(ang_cell)))
                # @show celldofs(ang_cell)
                append!(outflow_angular_dofs, celldofs(ang_cell))
            end
        end
    elseif boundary_node_coord ≈ RIGHT
        for ang_cell in CellIterator(dh_μ)
            if all( x-> x >= 0., Iterators.flatten(getcoordinates(ang_cell)))
                append!(outflow_angular_dofs, celldofs(ang_cell))
            end
        end
    end

    return outflow_angular_dofs
end
###############################################################################
###############################################################################





###############################################################################
###############################################################################

"""
    Assemble the angular flux vectors (∫ ψₘ dμ and ∫ μ²ψₘ dμ) by iterating over angular cells
"""
function AngularFluxAssembly!(ψ_flux, ψ_flux_μ, ψ_flux_μsquared, cellValues_μ::CellValues, dh_μ::DofHandler)
    for ang_cell in CellIterator(dh_μ)              # iterate over cells
        reinit!(cellValues_μ, ang_cell)
        # println(getcoordinates(ang_cell))   # coordiantes of a cell
        # println(celldofs(ang_cell))         # global dofs that the cell contributes 
        
        # For spatial coordinate of a quadrature point in a cell, 
        # use `spatial_coordinate(cellValues, qp, getcoordinates(cell))`
        
        dofs = celldofs(ang_cell);
        for qp in 1:getnquadpoints(cellValues_μ)    #iterate over quadrature points in a cell
            dμ = getdetJdV(cellValues_μ, qp)
            for j in 1:getnbasefunctions(cellValues_μ)
                ψ_flux[dofs[j]] += shape_value(cellValues_μ, qp, j) * dμ;

                μ = spatial_coordinate(cellValues_μ, qp, getcoordinates(ang_cell))[1];
                ψ_flux_μ[dofs[j]] += μ * shape_value(cellValues_μ, qp, j) * dμ;

                μsquared = μ^2;
                ψ_flux_μsquared[dofs[j]] += μsquared * shape_value(cellValues_μ, qp, j) * dμ;
            end
        end
    end

    return nothing

end

"""
    Assemble the angular mass matrices (∫ μ²ψₘψⱼ dμ) and (∫ ψₘψⱼ dμ) 
    by iterating over angular cells
"""
function AngularMassAssembly!(mass1_μsquared, mass2, mass_absμ, angular_mass_μ, cellValues_μ::CellValues, dh_μ::DofHandler)
    nbasefunctions = getnbasefunctions(cellValues_μ)
    assembler1 = start_assemble(mass1_μsquared)
    assembler2 = start_assemble(mass2)
    assembler3 = start_assemble(mass_absμ)
    assembler4 = start_assemble(angular_mass_μ)
    Ke1 = zeros(nbasefunctions, nbasefunctions)
    Ke2 = zeros(nbasefunctions, nbasefunctions)
    Ke3 = zeros(nbasefunctions, nbasefunctions)
    Ke4 = zeros(nbasefunctions, nbasefunctions)

    for ang_cell in CellIterator(dh_μ)
        fill!(Ke1, 0.0)
        fill!(Ke2, 0.0)
        fill!(Ke3, 0.0)
        fill!(Ke4, 0.0)
        reinit!(cellValues_μ, ang_cell)
        cellDofs = celldofs(ang_cell)

        for qp in 1:getnquadpoints(cellValues_μ)
            dμ = getdetJdV(cellValues_μ, qp)
            μ = spatial_coordinate(cellValues_μ, qp, getcoordinates(ang_cell))[1];
            # @show dμ
            for itp_j in 1:nbasefunctions
                for itp_k in 1:nbasefunctions
                    Ke1[itp_j, itp_k] += μ^2 * shape_value(cellValues_μ, qp, itp_j) * shape_value(cellValues_μ, qp, itp_k) * dμ
                    Ke2[itp_j, itp_k] += shape_value(cellValues_μ, qp, itp_j) * shape_value(cellValues_μ, qp, itp_k) * dμ
                    Ke3[itp_j, itp_k] += abs(μ) * shape_value(cellValues_μ, qp, itp_j) * shape_value(cellValues_μ, qp, itp_k) * dμ
                    Ke4[itp_j, itp_k] += μ * shape_value(cellValues_μ, qp, itp_j) * shape_value(cellValues_μ, qp, itp_k) * dμ
                end
            end
        end

        assemble!(assembler1, cellDofs, Ke1)
        assemble!(assembler2, cellDofs, Ke2)
        assemble!(assembler3, cellDofs, Ke3)
        assemble!(assembler4, cellDofs, Ke4)
    end

    return nothing
end


"""
    Assemble the matrices (∫ 1/σ∇Φₘ⋅∇Φⱼ dx) and (∫ σₛ/(2σ(σ-σₛ))∇Φₘ⋅∇Φⱼ dx) 
    by iterating over spatial cells
"""
function SpatialStiffnessAssembly!(Stiffness_1, Stiffness_2, cellValues_x::CellValues, dh_x::DofHandler)
    nbasefunctions = getnbasefunctions(cellValues_x)
    assembler1 = start_assemble(Stiffness_1)
    assembler2 = start_assemble(Stiffness_2)
    Ke1 = zeros(nbasefunctions, nbasefunctions)
    Ke2 = zeros(nbasefunctions, nbasefunctions)

    for cell in CellIterator(dh_x)
        fill!(Ke1, 0.0)
        fill!(Ke2, 0.0)
        reinit!(cellValues_x, cell)
        cellDofs = celldofs(cell)

        for qp in 1:getnquadpoints(cellValues_x)
            dV = getdetJdV(cellValues_x, qp)
            qp_coord = spatial_coordinate(cellValues_x, qp, getcoordinates(cell))[1]
            w1 = 1/σ(qp_coord)
            w2 = σₛ(qp_coord) / ( 2 * σ(qp_coord) * (σ(qp_coord) - σₛ(qp_coord)) )
            for itp_j in 1:nbasefunctions
                for itp_k in 1:nbasefunctions
                    ∇Φ_j = shape_gradient(cellValues_x, qp, itp_j)
                    ∇Φ_k = shape_gradient(cellValues_x, qp, itp_k)
                    Ke1[itp_j, itp_k] += w1 * ∇Φ_j ⋅ ∇Φ_k * dV
                    Ke2[itp_j, itp_k] += w2 * ∇Φ_j ⋅ ∇Φ_k * dV
                end
            end
        end
        
        assemble!(assembler1, cellDofs, Ke1)
        assemble!(assembler2, cellDofs, Ke2)
    end

    return nothing

end

"""
    Assemble the matrices (∫ σΦₘΦⱼ dx) and (1/2∫ σₛ∇Φₘ⋅∇Φⱼ dx) 
    by iterating over spatial cells
"""
function SpatialMassAssembly!(Mass_1, Mass_2, Mass_3, cellValues_x::CellValues, dh_x::DofHandler)
    nbasefunctions = getnbasefunctions(cellValues_x)
    assembler1 = start_assemble(Mass_1)
    assembler2 = start_assemble(Mass_2)
    assembler3 = start_assemble(Mass_3)
    Ke1 = zeros(nbasefunctions, nbasefunctions)
    Ke2 = zeros(nbasefunctions, nbasefunctions)
    Ke3 = zeros(nbasefunctions, nbasefunctions)

    for cell in CellIterator(dh_x)
        fill!(Ke1, 0.0)
        fill!(Ke2, 0.0)
        fill!(Ke3, 0.0)
        reinit!(cellValues_x, cell)
        cellDofs = celldofs(cell)

        for qp in 1:getnquadpoints(cellValues_x)
            dV = getdetJdV(cellValues_x, qp)
            qp_coord = spatial_coordinate(cellValues_x, qp, getcoordinates(cell))[1]
            w1 = σ(qp_coord)
            w2 = σₛ(qp_coord)/2
            for itp_j in 1:nbasefunctions
                for itp_k in 1:nbasefunctions
                    Φ_j = shape_value(cellValues_x, qp, itp_j)
                    Φ_k = shape_value(cellValues_x, qp, itp_k)
                    Ke1[itp_j, itp_k] += w1 * Φ_j * Φ_k * dV
                    Ke2[itp_j, itp_k] += w2 * Φ_j * Φ_k * dV
                    Ke3[itp_j, itp_k] += Φ_j * Φ_k * dV
                end
            end
        end
        
        assemble!(assembler1, cellDofs, Ke1)
        assemble!(assembler2, cellDofs, Ke2)
        assemble!(assembler3, cellDofs, Ke3)
    end

    return nothing

end

function assemble_global_matrix(
    dh_x, dh_μ, 
    outflow_dofs, outflow_map,
    spatial_stiffness_1, spatial_stiffness_2,
    spatial_mass_1, spatial_mass_2,
    angular_mass_μsquared, angular_mass_absμ, angular_mass,
    ψ_flux, ψ_flux_μ,
)
    nx = ndofs(dh_x)
    nμ = ndofs(dh_μ)
    N = nx * nμ

    angular_coupling = sparse((ψ_flux_μ * transpose(ψ_flux_μ)))
    angular_outer = sparse((ψ_flux * transpose(ψ_flux)))

    M = kron(spatial_stiffness_1, angular_mass_μsquared) +
        kron(spatial_stiffness_2, angular_coupling) +
        kron(spatial_mass_1, angular_mass) -
        kron(spatial_mass_2, angular_outer)

    rows = Int[]
    cols = Int[]
    values = Float64[]

    for (m, n) in outflow_dofs
        row = (m - 1) * nμ + n
        # is_free[row] || continue

        for k in outflow_map[m]
            val = angular_mass_absμ[n, k]
            if val != 0.0
                push!(rows, row)
                push!(cols, (m - 1) * nμ + k)
                push!(values, val)
            end
        end
    end

    M += sparse(rows, cols, values, N, N)
    dropzeros!(M)
    return M
end

###############################################################################
###############################################################################
@time begin

# Spatial domain x ∈ (0, 1); 1st order line elements
LEFT = 0.; RIGHT = 1.;
grid_x = generate_grid(Line, (ncell_x, ), Vec(LEFT, ), Vec(RIGHT, ))
itp_x  = Lagrange{RefLine, 1}()
qr_x   = QuadratureRule{RefLine}(2)
cellValues_x = CellValues(qr_x, itp_x)

# Angular domain μ ∈ (-1, 1); piecewise constant line elements
grid_μ = generate_grid(Line, (ncell_μ, ), Vec(-1., ), Vec(1., ))
itp_μ  = DiscontinuousLagrange{RefLine, 0}()
qr_μ   = QuadratureRule{RefLine}(2)
cellValues_μ = CellValues(qr_μ, itp_μ)

dh_x = DofHandler(grid_x)
add!(dh_x, :Φ, itp_x)
close!(dh_x)

dh_μ = DofHandler(grid_μ)
add!(dh_μ, :ψ, itp_μ)
close!(dh_μ)
@assert iseven(ndofs(dh_μ))

spatial_stiffness_1 = allocate_matrix(dh_x); # weighted by 1/σ
spatial_stiffness_2 = allocate_matrix(dh_x); # weighted by 1/(2(σ - σₛ))
spatial_mass_1      = allocate_matrix(dh_x); # weighted by σ
spatial_mass_2      = allocate_matrix(dh_x); # weighted by -σₛ/2
spatial_mass_3      = allocate_matrix(dh_x); # unweighted spatial mass matrix

angular_mass_μsquared = allocate_matrix(dh_μ); # weighted by μ²
angular_mass          = allocate_matrix(dh_μ); # unweighted angular mass matrix
angular_mass_absμ     = allocate_matrix(dh_μ); # weighted by |μ|
angular_mass_μ        = allocate_matrix(dh_μ); # weighted by μ
ψ_flux          = zeros(ndofs(dh_μ));
ψ_flux_μ        = zeros(ndofs(dh_μ));
ψ_flux_μsquared = zeros(ndofs(dh_μ));

SpatialStiffnessAssembly!(spatial_stiffness_1, spatial_stiffness_2, cellValues_x, dh_x)
SpatialMassAssembly!(spatial_mass_1, spatial_mass_2, spatial_mass_3, cellValues_x, dh_x)

AngularFluxAssembly!(ψ_flux, ψ_flux_μ, ψ_flux_μsquared, cellValues_μ, dh_μ)
AngularMassAssembly!(angular_mass_μsquared, angular_mass, angular_mass_absμ, angular_mass_μ, cellValues_μ, dh_μ)




boundary_facetidxs = union(getfacetset(grid_x, "left"), getfacetset(grid_x, "right"));

inflow_dofs::Vector{Tuple{Int, Int}}  = []
outflow_dofs::Vector{Tuple{Int, Int}} = []
spatial_boundary_cell_ids = Int[]
# iteration over boundary facets
for boundary_facetidx in boundary_facetidxs
    boundary_cell_id, facet_nodes_local_id  = boundary_facetidx.idx
    push!(spatial_boundary_cell_ids, boundary_cell_id)
    (boundary_node_id, boundary_node_coord) = get_boundary_point(grid_x, boundary_facetidx)

    # inflow and outflow angular DOFs for the boundary node
    angular_inflows  = inflow_angular_dofs1D(boundary_node_coord)
    angular_outflows = outflow_angular_dofs1D(boundary_node_coord)
    
    boundary_cell_dofs = celldofs(dh_x, boundary_cell_id) # e.g. [16, 17]
    local_dof_id       = Vec(Ferrite.facetdof_indices(itp_x)[facet_nodes_local_id])
    # global DOF of the boundary facet node
    boundary_facet_dof = boundary_cell_dofs[local_dof_id]   # facet_id-th DOF in the cell
    
    inflow_condition_dof = (boundary_facet_dof, angular_inflows)
    outflow_dof          = (boundary_facet_dof, angular_outflows)

    append!(inflow_dofs, 
            [(spatial_dof, ang_dof) for ang_dof in inflow_condition_dof[2] for spatial_dof in inflow_condition_dof[1]]
        )
    append!(outflow_dofs, 
            [(spatial_dof, ang_dof) for ang_dof in outflow_dof[2] for spatial_dof in outflow_dof[1]]
        )
    
end

spatial_boundary_cell_ids = unique(spatial_boundary_cell_ids)

groups = Dict()
paris = [el[1]=>el[2] for el in inflow_dofs]
for (k, v) in paris
    push!(get!(groups, k, []), v)
end
inflow_spatial_dof_to_inflow_ang_dof = Dict(k => v for (k, v) in groups)
inflow_spatial_dof = [k for k in keys(inflow_spatial_dof_to_inflow_ang_dof)]

groups = Dict()
paris = [el[1]=>el[2] for el in outflow_dofs]
for (k, v) in paris
    push!(get!(groups, k, []), v)
end
outflow_spatial_dof_to_outflow_ang_dof = Dict(k => v for (k, v) in groups)
outflow_spatial_dof = [k for k in keys(outflow_spatial_dof_to_outflow_ang_dof)]

inflow_dofs_set = Set(inflow_dofs)

# free_dofs = Iterators.filter(
#     pair -> pair ∉ inflow_dofs_set,
#     Iterators.product(1:ndofs(dh_x), 1:ndofs(dh_μ))
# )

free_dofs = Iterators.product(1:ndofs(dh_x), 1:ndofs(dh_μ))

M = assemble_global_matrix(
    dh_x, dh_μ, 
    outflow_dofs, outflow_spatial_dof_to_outflow_ang_dof,
    spatial_stiffness_1, spatial_stiffness_2,
    spatial_mass_1, spatial_mass_2,
    angular_mass_μsquared, angular_mass_absμ, angular_mass,
    ψ_flux, ψ_flux_μ,
)

# @test M2 ≈ M

load_vector = zeros(ndofs(dh_x) * ndofs(dh_μ))

# Inflow boundary condition
# defined for (x = 0 and μ > 0) and (x = 1 and μ < 0)
inflow_data = Dict()
for j in eachindex(inflow_spatial_dof)
    inflow_data[inflow_spatial_dof[j]] = ones(length(inflow_spatial_dof_to_inflow_ang_dof[inflow_spatial_dof[j]]))
end

# inflow trace
begin
for (m, n) in inflow_dofs
    for k in eachindex(inflow_spatial_dof_to_inflow_ang_dof[m])
        spatial_dof = m
        ang_dof = inflow_spatial_dof_to_inflow_ang_dof[spatial_dof][k]
        load_vector[(m-1)*ndofs(dh_μ) + n] += inflow_data[spatial_dof][k] * angular_mass_absμ[n, ang_dof]
    end
end
end

begin
sol_vec = M \ load_vector
end

sol = Matrix{Float64}(undef, ndofs(dh_x), ndofs(dh_μ))
for m in 1:ndofs(dh_x), n in 1:ndofs(dh_μ)
    sol[m, n] = sol_vec[(m - 1) * ndofs(dh_μ) + n]
end


end

# for j in eachindex(inflow_boundary_dof)
#     spatial_dof = inflow_boundary_dof[j];
#     ang_dofs = inflow_boundary_dof_to_inflow_ang_dof[spatial_dof];
#     for k in eachindex(ang_dofs)
#         ang_dof = ang_dofs[k];
#         sol[spatial_dof, ang_dof] = inflow_data[j, k];
#     end
# end

# heatmap(sol, c=:thermal, dpi = 200)
# plot(sol[1, :], dpi = 200)

fluxes = sol * ψ_flux
plot(fluxes, dpi = 200)

# μ_centers = [spatial_coordinate(cellValues_μ, 1, getcoordinates(cell))[1] for cell in CellIterator(dh_μ)]
# plot(μ_centers, sol[1, :], xlabel="μ", ylabel="φ(0, μ)", label="Weak Robin BC", lw=2, dpi=200)

# using BenchmarkTools
# @btime AngularMassAssembly!(K_μ, cellValues_μ, dh_μ) evals = 10 setup = (fill!(K_μ, 0.))
