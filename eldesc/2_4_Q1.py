name = '2_4_Q1'

desc = {
    'family'      : 'LagrangeTensor111',
    'approxOrder' : 1
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
