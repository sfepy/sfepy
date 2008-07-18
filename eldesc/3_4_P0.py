name = '3_4_P0'

desc = {
    'family'      : 'Lagrange',
    'approx_order' : 0
}

geometry = '3_4'

nodes = {
    'v' : {
        'mode'  : 'generate',
        'order' : 'mesi'
    },
    's3' : {
        'mode'  : 'generate',
        'order' : 'mei'
    }
}

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
