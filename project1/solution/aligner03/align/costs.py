import aligner03.utils

default_mutation_costs = aligner03.utils.minidict({
    # Deletion
    'D': -2,
    # Insertion
    'I': -2,
    # Mutation
    'X': -1,
    # Match
    '=': 3,
})
