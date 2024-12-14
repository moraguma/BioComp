while True:
    val = input().split(" ")
    s = (float(val[0]) + float(val[1]) + float(val[2])) / 2
    print(f"{s - float(val[1])} {s - float(val[2])} {s - float(val[0])}")
