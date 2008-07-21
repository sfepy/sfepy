name = '3_4_P1B'

desc = {
    'family'      : 'LagrangeP1B',
    'approx_order' : 1
}

geometry = '3_4'

nodes = {
    'v' : {
        'mode'  : 'generate',
        'order' : 'vefb',
    },
    's3' : {
        'mode'  : 'generate',
        'order' : 'veb'
    }
}

# geometry key : qp key
key_map = {'v' : 'v', 's3' : 's3'}

base_funs = {
    'v' : {
        'mode' : 'generate',
        'grad' : 1
    },
    's3' : {
        'mode' : 'generate',
        'grad' : 1
    }
}
