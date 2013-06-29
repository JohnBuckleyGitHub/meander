

def cont_checker():
    c = 2
    d = 3
    a = 'dog'
    b = 'cat'
    e = 7
    ls = [a, b, c, d, e]
    for i in ls:
        try:
            sump = i + 3
        except:
            print('could not add ', i)
            continue
        print(i)
        print('should not error on', sump)