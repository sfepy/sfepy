name = '3_4_P1'

desc = {
    'family'      : 'Lagrange',
    'approxOrder' : 1
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

# geometry key : qp key
keyMap = {'v' : 'v', 's3' : 's3'}

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
