

f = [0, 1, 2, 3, 4, 5]
print(f)

s = []
s.append(5)
s.append(4)

for id in sorted(s, reverse=True):
    del f[id]

print(f)