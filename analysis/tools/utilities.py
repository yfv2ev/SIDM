def print_list(l):
    print('\n'.join(l))


def print_debug(name, var, print_mode=True):
    if print_mode:
        print(f"{name}: {var}")


def partition_list(l, condition):
    passes = []
    fails = []
    for x in l:
        if condition(x):
            passes.append(x)
        else:
            fails.append(x)
    return passes, fails