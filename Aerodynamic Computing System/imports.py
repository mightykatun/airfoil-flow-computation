#read from config.txt

def var(query, form):
    file = open("config.txt", "r", encoding="utf8")
    variables = file.readlines()
    i = 0
    while True:
        temp = variables[i].split(": ", 1)
        if query in temp[0]:
            break
        i += 1
    temp[1] = temp[1].strip("\n")
    try:
        float(temp[1])
    except ValueError:
        return str(temp[1])
    else:
        if form == "int":
            return int(temp[1])
        elif form == "str":
            return str(int(temp[1]))
        elif form == "float":
            return float(temp[1])
