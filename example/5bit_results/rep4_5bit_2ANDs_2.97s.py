# [4278251760, 4042260480]
S = []
for X in range(32):
	x=[(X>>i)&1 for i in range(5)]
	g = [0]*2; y = [0]*5

	g[0] = (x[2]^x[3]) & (x[4]) ^ x[2]
	g[1] = (x[2]) & (x[4]) 

	y[0] = x[0] 
	y[1] = x[1] 
	y[2] = g[0] 
	y[3] = x[3] ^ g[1] 
	y[4] = x[4] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4))

print(S)
print([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 24, 25, 26, 27, 28, 29, 30, 31, 20, 21, 22, 23])
if tuple(S)==tuple([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 24, 25, 26, 27, 28, 29, 30, 31, 20, 21, 22, 23]):
	print(True)
else:
	print(False)
