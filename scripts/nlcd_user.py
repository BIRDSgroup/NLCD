from nlcd_main import *
def nlcd_group(data,shuffles,algo):
    output=[]
    for i in range(len(data)):
        L=data[i][0]
        A=data[i][1]
        B=data[i][2]
        out=combine_tests(L,A,B,shuffles,algo)
        output.append(out)
    return pd.DataFrame(output,columns=['p_final','p_LassocB','p_LassocA|B','p_AassocB|L','p_LindB|A','OS Test 2','OS Test 4'])

def nlcd_single(L,A,B,shuffles,algo):
    out=combine_tests(L,A,B,shuffles,algo)
    print("The final p value is ",out[0])
    print("Test 1 L assoc B ",out[1])
    print("Test 2 L assoc A | B ",out[2])
    print("Test 3 A assoc B | L ",out[3])
    print("Test 4 L ind B | A ",out[4])
    print("Overlap score from Test 2 ",out[5])
    print("Overlap score from Test 4 ",out[6])




