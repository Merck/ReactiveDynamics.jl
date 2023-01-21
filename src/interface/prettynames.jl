# parameter of agent classes, and their pretty names
const prettynames = Dict(:transition => Dict(:rate => [:rate],
                                             :post_action => [:post_action, :post],
                                             :name => [:name],
                                             :priority => [:priority],
                                             :probability_of_success => [
                                                 :probability,
                                                 :prob,
                                                 :pos,
                                                 :p,
                                             ],
                                             :capacity => [:cap, :capacity],
                                             :cycle_time => [:ct, :cycle_time],
                                             :max_life_time => [
                                                 :lifetime,
                                                 :maxlifetime,
                                                 :maxtime,
                                                 :timetolive,
                                             ]),
                         :species => Dict(:init_uncertainty => [:uncertainty, :randomness]),
                         :event => Dict(:name => [:name]),
                         :sampleable => Dict(:name => [:name], :every => [:every],
                                             :on => [:on]))
