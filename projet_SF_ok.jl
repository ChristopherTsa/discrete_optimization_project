using JuMP, Gurobi, Plots, Random
include("lecture_distances.jl")
file = "instances/instance_80_1.txt"
instance_name = split(basename(file), ".")[1]


# Fonction pour visualiser la solution optimale avec les régions
function visualize_solution(n::Int64, coords::Matrix{Int64}, arcs::Vector{Tuple{Int, Int}}, regions::Dict{Int64, Vector{Int64}})
    Random.seed!(1000)

    # Générer une couleur unique pour chaque région
    region_colors = Dict{Int64, RGB{Float64}}()
    for region in keys(regions)
        region_colors[region] = RGB(rand(), rand(), rand())  # Couleurs aléatoires pour chaque région
    end

    # Tracer les aérodromes en fonction de leur région
    p = plot(legend=:outerright, title="Chemin optimal pour l'$instance_name", xlabel="x", ylabel="y")
    for region in keys(regions)
        region_coords = coords[regions[region], :]
        scatter!(region_coords[:, 1], region_coords[:, 2], label="Région $region", color=region_colors[region], markersize=10)
    end

    # Ajouter les labels aux aérodromes
    for i in 1:n
        annotate!(coords[i, 1], coords[i, 2], text(string(i), :center, 12))
    end

    # Tracer les arcs utilisés
    for (i, j) in arcs
        x = [coords[i, 1], coords[j, 1]]
        y = [coords[i, 2], coords[j, 2]]
        plot!(x, y, arrow=:arrow, label="", linewidth=2, color=:red)
    end
    savefig(p, "results/chemin_optimal_$instance_name.png")
    display(p)
end


function resolution(n::Int64, D::Matrix{Int64}, R::Int64, Amin::Int64, d::Int64, f::Int64, coords::Matrix{Int64}, Nr::Int64, regions::Dict{Int64, Vector{Int64}})
    
    model = Model(Gurobi.Optimizer)

    # Déclaration des variables
    @variable(model, x[i=1:n, j=1:n], Bin)   # x[i,j] = 1 if traveling from i to j
    @variable(model, 0 <= q[i=1:n, j=1:n])  # Auxiliary flow variable for SF
    @variable(model, z[i=1:n; i!=d], Bin)             # Binary variable for nodes in path

    # Déclaration de l'objectif
    @objective(model, Min, sum(D[i, j] * x[i, j] for i in 1:n, j in 1:n))

    # Contraintes de distance et de nombre minimum d'aérodromes visités
    @constraint(model, [i in 1:n, j in 1:n], D[i, j] * x[i, j] <= R)
    @constraint(model, sum(sum(x[i,j] for j in 1:n if j != i) for i in 1:n) >= Amin - 1)

    # Contraintes de depart et arrivee
    @constraint(model, sum(x[d,j] for j in 1:n if j != d) - sum(x[j,d] for j in 1:n if j != d) == 1)
    @constraint(model, sum(x[f,j] for j in 1:n if j != f) - sum(x[j,f] for j in 1:n if j != f) == -1)

    # Contrainte de visite de chaque région au moins une fois
    @constraint(model, [k in keys(regions)], sum(sum(x[i,j] for i in 1:n) for j in regions[k]) >= 1)

    @constraint(model, [i in 1:n], x[i, i] == 0)

    @constraint(model, [i in 1:n], sum(x[i,j] for j in 1:n if j != i) <= 1)
    @constraint(model, [i in 1:n; i != d && i != f], sum(x[i,j] for j in 1:n if j != i) == sum(x[j,i] for j in 1:n if j != i))
    
    # Contraintes spécifiques à la formulation MCF
    # Contrainte de capacité des flows
    @constraint(model, [i in 1:n, j in 1:n; i != j], q[i, j] <= (n - 1) * x[i, j])
    
    # Contraintes d'équilibre des flows
    @constraint(model, sum(q[d, j] for j in 1:n if j != d) == sum(z))
    @constraint(model, [k in 1:n; k != d], sum(q[i, k] for i in 1:n if i != k) - sum(q[k, j] for j in 1:n if j != k) == z[k])
    @constraint(model, [k in 1:n; k != d], sum(x[i, k] for i in 1:n if i != k) == z[k])

    # Résolution
    start = time()
    #set_silent(model)
    optimize!(model)
    finish = time()
    println("Resolution time (s): ", finish - start)

    # Vérifier si le modèle a une solution optimale
    if termination_status(model) == MOI.OPTIMAL
        println("Objective value: ", objective_value(model))
        println("Optimal solution:")
        
        # Trouver les arcs utilisés
        arcs = Tuple{Int64, Int64}[]
        for i in 1:n
            for j in 1:n
                if value(x[i, j]) > 0.5  # Si l'arc est utilisé
                    push!(arcs, (i, j))
                    println("Arc: ", i, " -> ", j)
                end
            end
        end

        # Visualisation
        visualize_solution(n, coords, arcs, regions)
    else
        println("No optimal solution found.")
    end
end

# Fonction main
function main()
    n, d, f, Amin, Nr, R, regions, coords, D = readInstance(file)
    resolution(n, D, R, Amin, d, f, coords, Nr, regions)
end

main()