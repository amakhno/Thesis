import pandas as pd

def compare_h5(fn1, fn2):
    def get_index(store, a, z):
        if (isinstance(store, pd.HDFStore)):
            sa = pd.Series(store.root.A)
            sz = pd.Series(store.root.Z)
            a = sa[sa == a].keys()
            z = sz[sz == z].keys()
            index = a.intersection(z)
            if (len(index) == 1):
                return index[0]
            else:
                raise "Not found element"
        else:
            raise "Error"

    store1 = pd.io.pytables.HDFStore(fn1)
    store2 = pd.io.pytables.HDFStore(fn2)
    y1 = pd.DataFrame(store1.root.Y)
    y2 = pd.DataFrame(store2.root.Y)
    y3 = abs(y1 - y2)
    y3_sum = y3.sum().sum()
    store1.close()
    store2.close()
    return y3_sum

folder_name = 'out-r'
for index in range(1, 2):
    # print(folder_name + '/WithMy-' + str(index))
    diff = compare_h5(folder_name + '/WithMy-' + str(index) + '.h5', 
        folder_name + '/WithoutMy-' + str(index) + '.h5')
    print("multyplier: 1e{0}, diff: {1}".format(index, diff))