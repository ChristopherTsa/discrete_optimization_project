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
    # Étape 1 : Récupérer les arcs trouvés et créer un dictionnaire pour les prolongements
    println("Arcs utilisés dans la solution actuelle:")
    arc_dict = Dict{Int, Vector{Int}}()
    L = Set{Tuple{Int, Int}}()
    for i in 1:n
        for j in 1:n
            if x_vals[i, j] > 0.5
                println("Arc: $i -> $j")
                push!(L, (i, j))
                if haskey(arc_dict, i)
                    push!(arc_dict[i], j)
                else
                    arc_dict[i] = [j]
                end
            end
        end
    end

    # Initialiser la liste des sous-tours trouvés et l'ensemble des sous-tours détectés
    subtours = []
    detected_subtours = Set{Set{Int}}()

    # Boucle principale pour détecter les sous-tours
    while !isempty(L)
        # Prendre un arc de départ
        (start, next_node) = first(L)
        delete!(L, (start, next_node))  # Supprimer l'élément de L
        current_path = [(start, next_node)]
        visited_nodes = Set([start, next_node])

        # Chercher une prolongation jusqu'à ce qu'une boucle soit trouvée ou qu'il n'y ait plus de prolongation
        is_cycle = false
        while true
            # Chercher un arc qui prolonge le chemin
            found_extension = false
            if haskey(arc_dict, next_node)
                for to in arc_dict[next_node]
                    if to ∉ visited_nodes  # Prolonge le chemin sans revenir sur un nœud visité
                        push!(current_path, (next_node, to))
                        push!(visited_nodes, to)
                        next_node = to
                        delete!(L, (next_node, to))  # Enlever l'arc prolongé de L
                        found_extension = true
                        break
                    elseif to == start  # Boucle détectée
                        push!(current_path, (next_node, to))
                        is_cycle = true
                        delete!(L, (next_node, to))  # Enlever l'arc de la boucle de L
                        break
                    end
                end
            end

            # Sortir si aucun prolongement n'a été trouvé
            if !found_extension || is_cycle
                break
            end
        end

        # Ajouter le sous-tour si un cycle a été trouvé, sinon ignorer le chemin
        if is_cycle
            # Créer un ensemble des nœuds du cycle et vérifier s'il a déjà été détecté
            subtour_nodes = Set(vcat([node for (node, _) in current_path], [current_path[end][2]]))
            if subtour_nodes ∉ detected_subtours
                push!(subtours, subtour_nodes)
                push!(detected_subtours, subtour_nodes)  # Marquer le sous-tour comme détecté
                println("Sous-tour détecté: ", subtour_nodes)
            end
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

    #contraint de continuite de flow
    @constraint(model, [i in 1:n], sum(x[i,j] for j in 1:n if j != i) <= 1)
    @constraint(model, [i in 1:n; i != d && i != f], sum(x[i,j] for j in 1:n if j != i) == sum(x[j,i] for j in 1:n if j != i))

    # Contraintes de départ et d'arrivée
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