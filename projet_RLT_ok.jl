using JuMP, Gurobi, Plots, Random
include("lecture_distances.jl")
file = "instances/instance_100_1.txt"
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
    
    #declaration du modele
    model = Model(Gurobi.Optimizer)

    #declaration des variables
    @variable(model, x[i=1:n, j=1:n], Bin)  # x[i,j] = 1 si on va de i à j
    @variable(model, a[i=1:n, j=1:n; i != j] >= 0)  # α_{ij} = t_j * x_{ij}
    @variable(model, b[i=1:n, j=1:n; i != j] >= 0)  # β_{ij} = t_i * x_{ij}
    @variable(model, u[i=1:n] >= 0)

    #declaration de l objectif
    @objective(model, Min, sum(D[i,j] * x[i,j] for i in 1:n, j in 1:n))

    #declaration des contraintes
    #contrainte de distance max
    @constraint(model, [i in 1:n, j in 1:n], D[i,j]*x[i,j] <= R)
    
    #contrainte de nombre min d'aerodrome visite
    @constraint(model, sum(sum(x[i,j] for j in 1:n if j != i) for i in 1:n) >= Amin - 1)

    #contraint de continuite de flow
    @constraint(model, [i in 1:n], sum(x[i,j] for j in 1:n if j != i) <= 1)
    @constraint(model, [i in 1:n; i != d && i != f], sum(x[i,j] for j in 1:n if j != i) == sum(x[j,i] for j in 1:n if j != i))

    #contrainte de depart et arrivee
    @constraint(model, sum(x[d,j] for j in 1:n if j != d) - sum(x[j,d] for j in 1:n if j != d) == 1)
    @constraint(model, sum(x[f,j] for j in 1:n if j != f) - sum(x[j,f] for j in 1:n if j != f) == -1)

    #contrainte de region
    @constraint(model, [k in keys(regions)], sum(sum(x[i,j] for i in 1:n) for j in regions[k]) >= 1)

    @constraint(model, [i in 1:n], x[i,i] == 0)
    
    #contraintes RLT
    @constraint(model, [i in 1:n, j in 1:n; i != j], a[i, j] == b[i, j] + x[i, j])
    @constraint(model, [j in 1:n; j != d && j != f], 
        x[d, j] + sum(a[i, j] for i in 1:n if i != j && i != d) - sum(b[j, i] for i in 1:n if i != j) == 0)
    @constraint(model, [i in 1:n, j in 1:n; i != j && i != d], a[i, j] <= (n - 1) * x[i, j])
    @constraint(model, [i in 1:n, j in 1:n; i != j && i != d], x[i, j] <= b[i, j])

    #contraintes de position
    @constraint(model, u[d] == 0)  # Le nœud de départ est en position 1
    @constraint(model, [i in 1:n, j in 1:n; i != j], u[i] + 1 <= u[j] + n * (1 - x[i,j]))
    @constraint(model, [i in 1:n, j in 1:n; i != j], u[j] + 1 <= u[i] + n * (1 - x[j,i]))

    #resolution
    # Start a chronometer
    start = time()
    #set_silent(model)
    optimize!(model)
    finish = time()
    println("Temps de résolution (s) :", finish-start)

    #affichage des resultats
    if termination_status(model) == MOI.OPTIMAL
        println("Objective value: ", objective_value(model))
        println("Solution optimale: ")
        
        # Affichage des arcs utilisés
        arcs = Tuple{Int64, Int64}[]
        for i in 1:n
            for j in 1:n
                if value(x[i,j]) > 0.5  # Si l'arc est utilisé
                    push!(arcs, (i, j))
                    println("Arc: ", i, " -> ", j)
                end
            end
        end

        # Visualisation
        visualize_solution(n, coords, arcs, regions)
    else
        println("Aucune solution optimale trouvée.")
    end


end



function main()
    # Lecture de l'instance
    n,d,f,Amin,Nr,R,regions,coords,D = readInstance(file)
    
    # Appel de la fonction de résolution
    resolution(n, D, R, Amin, d, f, coords, Nr, regions)
end

main()