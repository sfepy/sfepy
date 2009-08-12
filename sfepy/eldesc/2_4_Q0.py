name = '2_4_Q0'

desc = {
    'family'      : 'LagrangeTensor111',
    'approx_order' : 0
}

geometry = '2_4'

nodes = {
    'v' : {
        'mode'  : 'generate',
        'order' : 'mei'
    },
    's2' : {
        'mode'  : 'generate',
        'order' : 'mi'
    }
}

# geometry key : qp key
key_map = {'v' : 'v', 's2' : 's2'}

base_funs = {
    'v' : {
        'mode' : 'generate',
        'grad' : 1
    },
    's2' : {
        'mode' : 'generate',
        'grad' : 1
    }
} 
