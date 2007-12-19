name = '2_3_P1B'

desc = {
    'family'      : 'LagrangeP1B',
    'approxOrder' : 1
}

geometry = '2_3'

nodes = {
    'v' : {
        'mode'  : 'generate',
        'order' : 'veb',
    },
    's2' : {
        'mode'  : 'generate',
        'order' : 'vb'
    }
}

# geometry key : qp key
keyMap = {'v' : 'v', 's2' : 's2'}

baseFuns = {
    'v' : {
        'mode' : 'generate',
        'grad' : 1
    },
    's2' : {
        'mode' : 'generate',
        'grad' : 1
    }
}
