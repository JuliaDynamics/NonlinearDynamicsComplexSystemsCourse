# %% Exercise 1 (playing around with fox rabbit)
# The exercise uses the code of the main notebook.
# Here we're just pasting a code that leads to stable oscillations!

# Stable regime:
plot_population_timeseries(1000;
    total_rabbits = 100,
    rabbit_energy = 15.0, rabbit_Δenergy = 5.0, rabbit_r_prob = 0.50,
    # Fox properties:
    total_foxes = 5,
    fox_energy = 20.0, fox_Δenergy = 2.0, fox_r_prob = 0.03,
    # Grass properties:
    griddims = (20, 20), regrowth_time = 5,
    # General:
    seed = 1234555
)

# %% Exercise 2 (Wright fischer)
using Agents, Random

@agent Haploid NoSpaceAgent begin
    trait::Float64
end

function initialize_model(; population_size = 100, seed = 23182)
    rng = MersenneTwister(seed)
    model = AgentBasedModel(Haploid; rng, scheduler = Schedulers.randomly)
    #Add agents to the model with a random value for the trait
    for _ in 1:population_size
        add_agent!(model, rand(abmrng(model)))
    end
    return model
end

neutralmodel = initialize_model()

# Neutral sampling
modelstep_neutral!(model::ABM) = sample!(model, nagents(model))
# Collect data
using Statistics: mean
neutralmodel = initialize_model()
steps = 1000
adata = [(:trait, mean)]
adf, _ = run!(neutralmodel, dummystep, modelstep_neutral!, steps; adata)

# model with selection
modelstep_selection!(model::ABM) = sample!(model, nagents(model), :trait)

model_with_selection = initialize_model()
steps = 1000
adata = [(:trait, mean)]
adf, _ = run!(model_with_selection, dummystep, modelstep_selection!, steps; adata)


# %% Exercise 3: Spatial rock paper scissors
using Agents, Random

# The 3 subpopulations are arranged in a two-dimensional square lattice with periodic
# boundary conditions. An agent can only interact with its 4 nearest neighbours,
# therefore we use `metric = :manhattan`. We also use `GridSpaceSingle`, because
# in this model there cannot be two agents in the same position.
# (In fact, a more performant, but less intuitive way to write this ABM would be to not use
# agents at all, and only use a spatial property that represents the strategy, as is
# done in the Forest Fire model: )
# https://juliadynamics.github.io/AgentsExampleZoo.jl/dev/examples/forest_fire/#Forest-fire

dims = (50, 50)
space = GridSpaceSingle(dims, periodic = true, metric = :manhattan)

# Agent type
@agent Strategy GridAgent{2} begin
    #:rock, :paper or :scissors
    type::Symbol
end

Rock(id, pos) = Strategy(id, pos, :rock)
Paper(id, pos) = Strategy(id, pos, :paper)
Scissors(id, pos) = Strategy(id, pos, :scissors)

# When we initialize the model, we will use `GridSpaceSingle`
# as the space of choice.
function initialize_model(;
        dims = (60, 60),
        n_rock = 800,
        n_paper = n_rock,
        n_scissors = n_rock,
        selection_rate = 1.0,
        reproduction_rate = 1.0,
        exchange_rate = 0.7,
        seed = 23182,
    )

    rng = MersenneTwister(seed)
    space = GridSpaceSingle(dims, periodic = true, metric = :manhattan)
    #model properties
    properties = (
        selection_rate = selection_rate,
        reproduction_rate = reproduction_rate,
        exchange_rate = exchange_rate,
    )
    model = AgentBasedModel(Strategy, space;
        properties, rng, scheduler = Schedulers.randomly
    )
    # Add agents to the model to a random but unoccupied position
    # using the `add_agent_single!` function
    for (strategy, N) in zip((:rock, :paper, :scissors), (n_rock, n_paper, n_scissors))
        for _ in 1:N
            add_agent_single!(model, strategy)
        end
    end
    return model
end

model = initialize_model()

# Time evolution rules:
function strategy_step!(strategy, model)
    # propensities
    a1 = model.selection_rate * nagents(model)
    a2 = model.reproduction_rate * nagents(model)
    a3 = model.exchange_rate * nagents(model)
    # total propensities
    a0 = a1 + a2 + a3
    p = rand(abmrng(model))
    if 0 <= p <= a1/a0
        selection!(strategy, model)
    elseif a1/a0 <= p < (a1+a2)/a0
        reproduce!(strategy, model)
    elseif (a1+a2)/a0 <= p < 1
        move_or_swap!(strategy, model)
    end
end

function selection!(strategy, model)
    contender = random_nearby_agent(strategy, model)
    if !isnothing(contender)
        if strategy.type == :rock && contender.type == :scissors
            kill_agent!(contender, model)
        elseif strategy.type == :scissors && contender.type == :paper
            kill_agent!(contender, model)
        elseif strategy.type == :paper && contender.type == :rock
            kill_agent!(contender, model)
        end
    end
    return
end

function reproduce!(strategy, model)
    for pos in nearby_positions(strategy, model)
        # empty position to put an offspring in
        if isempty(pos, model)
            add_agent!(strategy.pos, model, strategy.type)
        end
        break
    end
    return
end

function move_or_swap!(strategy, model)
    # get a random position
    # (in the future, this will be a dedicated function from Agents.jl)
    function random_nearby_position(pos, model)
        return rand(abmrng(model), collect(nearby_positions(pos, model)))
    end
    rand_pos = random_nearby_position(strategy.pos, model)
    if isempty(rand_pos, model)
        move_agent!(strategy, rand_pos, model)
    else
        # To swap two agents we don't need to move or do anything complicated,
        # we simply swap their strategies!
        neighbor = model[id_in_position(rand_pos, model)]
        strategy.type, neighbor.type = neighbor.type, strategy.type
    end
    return
end

step!(model, strategy_step!)
model

step!(model, strategy_step!, 10)
model

# There is a bug here and the model gets more agents than positions.
# Will be fixed soon! Follow:
# https://github.com/JuliaDynamics/Agents.jl/issues/757

# %%
# Visualize
using InteractiveDynamics
# This model is too costly to animate in real time. We make a video
using CairoMakie
CairoMakie.activate!()

function strategycolor(a)
    if a.type == :rock
        :blue
    elseif a.type == :paper
        :orange
    else
        :red
    end
end

plotkwargs = (;
    ac = strategycolor,
    am = :rect,
)

abmvideo(
    "rockpaperscissors.mp4",
    model,
    strategy_step!;
    frames = 100,
    framerate = 8,
    title = "Spatial Rock-Paper-Scissors game",
    plotkwargs...,
)
