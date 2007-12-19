name = '3_4_P2B'

desc = {
    'family'      : 'LagrangeP2B',
    'approxOrder' : 2
}

geometry = '3_4'

nodes = {
    'v' : {
        'mode'  : 'generate',
        'order' : 'vefb',
    },
    's6' : {
        'mode'  : 'generate',
        'order' : 'veb'
    }
}

# geometry key : qp key
keyMap = {'v' : 'v', 's3' : 's6'}

baseFuns = {
    'v' : {
        'mode' : 'generate',
        'grad' : 1
    },
    's6' : {
        'mode' : 'generate',
        'grad' : 1
    }
}
