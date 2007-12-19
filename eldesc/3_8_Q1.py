name = '3_8_Q1'

desc = {
    'family'      : 'LagrangeTensor111',
    'approxOrder' : 1
}

geometry = '3_8'

nodes = {
    'v' : {
        'mode'  : 'generate',
        'order' : 'mesi'
    },
    's4' : {
        'mode'  : 'generate',
        'order' : 'mei'
    }
}

# geometry key : qp key
keyMap = {'v' : 'v', 's4' : 's4'}

baseFuns = {
    'v' : {
        'mode' : 'generate',
        'grad' : 1
    },
    's4' : {
        'mode' : 'generate',
        'grad' : 1
    }
}
