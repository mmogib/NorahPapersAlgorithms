using Polyhedra, CDDLib

H = hrep([HalfSpace([-1, -1], -2), HalfSpace([-1, 0], 0), HalfSpace([0, -1], 0)])
P = polyhedron(H, CDDLib.Library());
vr = vrep(P)
vrmat = MixedMatVRep(vr)
lva = [vrmat.V vrmat.R]
find_verex_index(v, A) = begin
    _, n = size(A)
    findfirst([lva[:, col] == v for col in 1:n])
end
# vreByDD = doubledescription(H).
adjdata = Dict{Int,Vector{Int}}()
println("new ..........")
for pi in eachindex(points(vr))
    current_vertix = get(vr, pi)
    incident_hspaces_idxes = incidenthalfspaceindices(P, pi)
    println("Start")
    @show current_vertix

    println("this is vertex number $(find_verex_index(current_vertix,lva))")
    @show incidenthalfspaces(P, pi)
    adjdata[find_verex_index(current_vertix, lva)] = []
    for hs in incident_hspaces_idxes
        @show get(P, hs)
        incident_points = incidentpoints(P, hs)
        incident_rays = incidentrays(P, hs)
        @show incident_rays
        if length(incident_points) > 1
            # @show incident_points
            adjacent_vertices = filter(pnt -> pnt != current_vertix, incident_points)
            @show adjacent_vertices[1]
            println("vertex number $(find_verex_index(adjacent_vertices[1],lva))")
            adjdata[find_verex_index(current_vertix, lva)] = push!(adjdata[find_verex_index(current_vertix, lva)], find_verex_index(adjacent_vertices[1], lva))
            println("----")
        else
            for i in eachindex(incident_rays)
                another_vertex = current_vertix + incident_rays[i].a
                println("vertex number $(find_verex_index(incident_rays[i].a,lva))")
                adjdata[find_verex_index(current_vertix, lva)] = push!(adjdata[find_verex_index(current_vertix, lva)], find_verex_index(incident_rays[i].a, lva))
                @show another_vertex
            end
            println("****")
        end
    end
    println("End")
end

for pi in eachindex(rays(vr))
    @show get(vr, pi)
end