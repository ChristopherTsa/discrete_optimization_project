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

# Fonction pour détecter et renvoyer les sous-tours
function find_subtours(x_vals, n)
    # Initialiser la liste de sous-tours
    subtours = []

    # Créer une liste de nœuds non encore associés à un sous-tour
    remaining_nodes = Set(1:n)

    model1 = Model(Gurobi.Optimizer)

    while !isempty(remaining_nodes)

        # Variables binaires z[i] indiquant si i appartient au sous-ensemble S
        z = @variable(model1, [i in remaining_nodes], Bin)

        # Objectif : maximiser la somme des arcs x_vals[i, j] qui font partie d'un sous-tour potentiel
        @objective(model1, Max, sum(x_vals[i, j] * z[i] * z[j] for i in remaining_nodes, j in remaining_nodes) - sum(z[i] for i in remaining_nodes) + 1)

        # Résoudre le modèle de séparation
        set_silent(model1)
        JuMP.optimize!(model1)


        # Vérifier si un sous-tour a été trouvé
        if objective_value(model1) > 0  # Si la valeur est positive, un sous-tour existe
            subtour_nodes = Set{Int}()

            for i in remaining_nodes 
                if value(z[i]) > 0  # Si z[i] = 1, alors i est dans le sous-tour
                    push!(subtour_nodes, i)
                end
            end

            # Ajouter le sous-tour détecté à la liste des sous-tours
            if !isempty(subtour_nodes)
                println("Sous-tour détecté: ", subtour_nodes)
                push!(subtours, subtour_nodes)
                # Retirer les nœuds du sous-tour de remaining_nodes
                remaining_nodes = setdiff(remaining_nodes, subtour_nodes)
            elseif isempty(subtour_nodes)
                println("Fin de la détection de sous-tours")
                break
            end
        else
            # Si aucun sous-tour n'est trouvé, sortir de la boucle
            break
        end
    end

    return subtours
end

# Fonction pour ajouter des contraintes de sous-tours
function add_subtour_constraints!(model, subtours, x)
    for subtour in subtours
        # Créer une contrainte pour interdire ce sous-tour spécifique
        @constraint(model, sum(x[i, j] for i in subtour, j in subtour if i != j) <= length(subtour) - 1)
        println("Ajout de la contrainte de coupure pour le sous-tour: ", subtour)
    end
end

function add_initial_constraints(model, x, n)
    # Contrainte pour éliminer les sous-tours de 2 nœuds
    # Empêcher un cycle de deux nœuds : i -> j -> i
    @constraint(model, [i in 1:n, j in 1:n; i !=j], x[i, j] + x[j, i] <= 1)

    # Contrainte pour éliminer les sous-tours de 3 nœuds (triangles)
    # Empêcher un cycle de trois nœuds : i -> j -> k -> i
    @constraint(model, [i in 1:n, j in 1:n, k in j+1:n], x[i, j] + x[j, k] + x[k, i] <= 2)
    @constraint(model, [i in 1:n, j in 1:n, k in j+1:n], x[i, k] + x[k, j] + x[j, i] <= 2)
end

# Fonction de résolution avec génération de contraintes
function resolution(n, D, R, Amin, d, f, coords, Nr, regions)
    model = Model(Gurobi.Optimizer)

    # Déclaration des variables
    @variable(model, x[i=1:n, j=1:n], Bin)  # x[i,j] = 1 si on va de i à j

    # Déclaration de l'objectif
    @objective(model, Min, sum(D[i, j] * x[i, j] for i in 1:n, j in 1:n))

    # Contraintes de distance et de nombre minimum d'aérodromes visités
    @constraint(model, [i in 1:n, j in 1:n], D[i, j] * x[i, j] <= R)
    @constraint(model, sum(sum(x[i, j] for j in 1:n if j != i) for i in 1:n) >= Amin - 1)

    # Contraintes de continuite de flow
    @constraint(model, [i in 1:n], sum(x[i,j] for j in 1:n if j != i) <= 1)
    @constraint(model, [i in 1:n; i != d && i != f], sum(x[i,j] for j in 1:n if j != i) == sum(x[j,i] for j in 1:n if j != i))

    # Contraintes de depart et arrivee
    @constraint(model, sum(x[d,j] for j in 1:n if j != d) - sum(x[j,d] for j in 1:n if j != d) == 1)
    @constraint(model, sum(x[f,j] for j in 1:n if j != f) - sum(x[j,f] for j in 1:n if j != f) == -1)

    # Contrainte de visite de chaque région au moins une fois
    @constraint(model, [k in keys(regions)], sum(sum(x[i,j] for i in 1:n) for j in regions[k]) >= 1)

    @constraint(model, [i in 1:n], x[i,i] == 0)

    # Ajouter les contraintes initiales pour éliminer les sous-tours simples
    add_initial_constraints(model, x, n)

    # Résolution avec génération de contraintes
    iteration = 0
    start = time()
    #set_silent(model)
    while true
        println("Itération: ", iteration)
        JuMP.optimize!(model)

        # Vérifier si le modèle a une solution optimale
        if termination_status(model) != MOI.OPTIMAL
            println("Aucune solution optimale trouvée.")
            return
        end

        # Récupérer les valeurs des variables x et détecter les sous-tours
        x_vals = value.(x)
        subtours = find_subtours(x_vals, n)

        # Si aucun sous-tour n'est trouvé, la solution est optimale
        if isempty(subtours)
            println("Solution optimale trouvée.")
            break
        end

        # Sinon, ajouter les contraintes de coupure pour chaque sous-tour
        println("Sous-tours détectés: ", subtours)
        add_subtour_constraints!(model, subtours, x)
        iteration += 1
    end

    finish = time()
    println("Temps de résolution (s) :", finish-start)

    # Affichage des résultats
    println("Objective value: ", objective_value(model))
    println("Solution optimale: ")
    
    # Affichage des arcs utilisés
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
end

# Fonction principale
function main()
    # Lecture de l'instance
    n, d, f, Amin, Nr, R, regions, coords, D = readInstance(file)
    
    # Appel de la fonction de résolution avec génération de contraintes
    resolution(n, D, R, Amin, d, f, coords, Nr, regions)
end

main()