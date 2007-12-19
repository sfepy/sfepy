# last revision: 19.03.2007
name = '2_3_P2'

desc = {
    'family'      : 'Lagrange',
    'approxOrder' : 2
}

geometry = '2_3'

nodes = {
    'v' : {
        'mode'  : 'generate',
        'order' : 'veb',
    },
    's3' : {
        'mode'  : 'generate',
        'order' : 'vb'
    }
}

# geometry key : qp key
keyMap = {'v' : 'v', 's2' : 's3'}

baseFuns = {
    'v' : {
        'mode' : 'generate',
        'grad' : 1
    },
    's3' : {
        'mode' : 'generate',
        'grad' : 1
    }
}
