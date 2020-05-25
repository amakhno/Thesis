
array_of_interesting_nucldies_mother = [
    "cd106",
    "cd108",
    " kr78",
    " kr80",
    "ba130",
    "ba132",
    "er164",
    "xe124",
    "xe126",
    "sn112",
    "sn114",
    "ce136",
    " sr84",
    "te120"
]

array_of_interesting_nucldies_pramother = [' br78',
' br80',
' rb78',
' rb80',
' rb84',
'  y84',
'ag106',
'ag108',
'in106',
'in108',
'in112',
'in114',
'sb112',
'sb114',
'sb120',
' i120',
' i124',
' i126',
'cs124',
'cs126',
'cs130',
'cs132',
'la130',
'la132',
'la136',
'pr136',
'ho164',
'tm164']

def get_product(lines, arr):
    reac_type = int(lines[0].replace('\n', ''))
    # nuclides = lines[1].split()
    n = 5
    line = lines[1]
    nuclides = [line[i:i+n] for i in range(0, len(line), n)]
    nuclides.remove("     ")
    if (reac_type == 1):
        nuc_input = [nuclides[0]]
        products = [nuclides[1]]
    if (reac_type == 2):        
        nuc_input = [nuclides[0]]
        products = nuclides[1:3]
    if (reac_type == 3):        
        nuc_input = [nuclides[0]]
        products = nuclides[1:4]
    if (reac_type == 4):        
        nuc_input = nuclides[0:2]
        products = [nuclides[2]]
    if (reac_type == 5):
        nuc_input = nuclides[0:2]
        products = nuclides[2:4]
    if (reac_type == 6):
        nuc_input = nuclides[0:2]
        products = nuclides[2:5]
    if (reac_type == 7):
        nuc_input = nuclides[0:2]
        products = nuclides[2:6]
    if (reac_type == 8):
        nuc_input = nuclides[0:3]
        products = [nuclides[3]]
    if (reac_type == 9):
        nuc_input = nuclides[0:3]
        products = nuclides[3:5]
    if (reac_type == 10):
        nuc_input = nuclides[0:4]
        products = nuclides[4:6]
    if (reac_type == 11):        
        nuc_input = nuclides[0]
        products = nuclides[1:5]
    if (len(list(set(products) & set(arr)))):
        print(str(nuc_input) + '->' + str(products))
        return [nuc_input, products]
    return []    

def read_file(file_name):
    if (isinstance(file_name, str)):
        with open(file_name, 'r') as f:
            lines = f.readlines()
        lines_count = len(lines)
        print(lines_count)
        arr_of_mother = []
        arr_of_pramother = []
        if (lines_count < 4):
            raise Exception("There are no any reactions in the file")
        for line_number in range(0, lines_count, 4):
            inp_and_prod = [] # get_product(lines[line_number: line_number + 4], array_of_interesting_nucldies_mother)
            if (len(inp_and_prod) > 0):
                arr_of_mother.append(inp_and_prod)
            inp_and_prod = get_product(lines[line_number: line_number + 4], array_of_interesting_nucldies_pramother)
            if (len(inp_and_prod) > 0):
                arr_of_pramother.append(inp_and_prod)
        print("Eyahl")
    else:
        raise Exception("Input name is not str")

file_name = 'reaclib'
read_file(file_name)