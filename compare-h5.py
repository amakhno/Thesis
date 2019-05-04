import pandas as pd

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

dictionary = {
    "cd": 48,
    "kr": 36,
    "ba": 56,
    "er": 68,
    "xe": 54,
    "sn": 50,
    "ce": 58,
    "sr": 38,
    "te": 43
}

def compare_h5(fn1, fn2):
    def get_index(store, a, z):
        if (isinstance(store, pd.HDFStore)):
            sa = pd.Series(store.root.A)
            sz = pd.Series(store.root.Z)
            a_arr = sa[sa == a].keys()
            z_arr = sz[sz == z].keys()
            index = a_arr.intersection(z_arr)
            if (len(index) == 1):
                return index[0]
            else:
                raise "Not found element"
        else:
            raise "Error"

    store1 = pd.io.pytables.HDFStore(fn1)
    store2 = pd.io.pytables.HDFStore(fn2)
    #y1 = pd.DataFrame(store1.root.Y)
    #y2 = pd.DataFrame(store2.root.Y)
    
    for index in range(0, len(array_of_interesting_nucldies_mother)):
        elem_name = array_of_interesting_nucldies_mother[index]
        z = dictionary[''.join(x for x in elem_name if x.isalpha())]
        a = int(''.join(x for x in elem_name if x.isdigit()))
        elem_index = get_index(store1, a, z)
        size_t1, _ = store1.root.Y.shape
        y1 = store1.root.Y[size_t1 - 1][elem_index]
        size_t2, _ = store2.root.Y.shape
        y2 = store2.root.Y[size_t2 - 1][elem_index]
        print("{0}, WithMy: {1}, WithoutMy: {2}".format(elem_name, y1, y2))
    #y3 = abs(y1 - y2)
    #y3_sum = y3.sum().sum()
    store1.close()
    store2.close()
    #return y3_sum

folder_name = 'out-r'
for index in range(1, 2):
    # print(folder_name + '/WithMy-' + str(index))
    compare_h5(folder_name + '/WithMy-' + str(index) + '.h5', 
        folder_name + '/WithoutMy-' + str(index) + '.h5')
    # print("multyplier: 1e{0}, diff: {1}".format(index, diff))