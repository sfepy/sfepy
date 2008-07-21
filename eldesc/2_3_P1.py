name = '2_3_P1'

desc = {
    'family'      : 'Lagrange',
    'approx_order' : 1
}

geometry = '2_3'

nodes = {
    'v' : {
        'mode'  : 'generate',
        'order' : 'veb'
    },
    's2' : {
        'mode'  : 'generate',
        'order' : 'vb'
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
