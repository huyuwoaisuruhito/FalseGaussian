PT = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al':13, 'Si': 14, 'P':15, 'S': 16, 
    'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 
}

IPT = {v:k for k,v in PT.items()}

def get_atomic_number(s):
    return PT[s]

def get_name(n):
    return IPT[n]