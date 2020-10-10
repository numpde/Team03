import humdum.utils

default_mutation_costs = humdum.utils.minidict({
    # Deletion
    'D': -2,
    # Insertion
    'I': -2,
    # Mutation
    'X': -1,
    # Match
    '=': 3,
})
